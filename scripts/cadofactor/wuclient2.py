#!/usr/bin/env python

import sys
import os
import fcntl
import errno
import stat
import optparse
import shutil
import time
if sys.version_info[0] == 3:
    # pylint: disable=E0611
    # pylint: disable=F0401
    import urllib.request as urllib_request
    import urllib.error as urllib_error
elif sys.version_info[0] == 2:
    import urllib2 as urllib_request
    import urllib2 as urllib_error
import subprocess
import hashlib
import logging
import socket
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import email.encoders
import email.generator
from string import Template
from io import BytesIO
from workunit import Workunit

if sys.version_info[0] == 3 and sys.version_info[1] < 4:
    # In Python 3.[012], use a fixed BytesGenerator which accepts a bytes 
    # input. The fact that the BytesGenerator in Python 3.[012] doesn't 
    # is a bug, see http://bugs.python.org/issue16564
    print ("Using work-around")
    class FixedBytesGenerator(email.generator.BytesGenerator):
        # pylint: disable=W0232
        # pylint: disable=E1101
        # pylint: disable=E1102
        # pylint: disable=E1002
        def _handle_bytes(self, msg):
            payload = msg.get_payload()
            if payload is None:
                return
            if isinstance(payload, bytes):
                # Payload is bytes, output is bytes - just write them
                self._fp.write(payload)
            elif isinstance(payload, str):
                super(FixedBytesGenerator, self)._handle_text(msg)
            else:
                # Payload is neither bytes nor string - this can't be right
                raise TypeError('bytes payload expected: %s' % type(payload))
        _writeBody = _handle_bytes
elif sys.version_info[0] == 2:
    # In Python 2.x, use the regular email generator
    # pylint: disable=C0103
    FixedBytesGenerator = email.generator.Generator
elif tuple(sys.version_info[:2]) >= (3, 4):
    # The tuple() guarantees a tuple >= tuple comparison; comparing 
    # tuple>=list does not work as desired
    # In Python >=3.3, use the bug-fixed bytes generator
    # pylint: disable=E1101
    # pylint: disable=C0103
    FixedBytesGenerator = email.generator.BytesGenerator


class WuMIMEMultipart(MIMEMultipart):
    ''' Defines convenience functions for attaching files and data to a 
    MIMEMultipart object
    '''
    
    def attach_data(self, name, filename, data):
        ''' Attach the data as a file

        name is the string that is sent to the server as the name of the form 
        input field for the upload; for us it is always 'result'.
        filename is the string that is sent to the server as the source file 
        name, this is the name as given in the RESULT lines, or stdout1 for 
        stdout of the first command that ran, etc.
        data is the content of the file to send.
        '''
        result = MIMEApplication(data, _encoder=email.encoders.encode_noop)
        result.add_header('Content-Disposition', 'form-data', 
                          name=name, filename=filename)
        self.attach(result)
    
    def attach_file(self, name, filename, filepath):
        ''' Attach the file as file

        Parameters as in attach_data(), but filepath is the path to the file 
        whose data should be sent
        '''
        logging.debug ("Adding result file %s to upload", filepath)
        with open(filepath, "rb") as infile:
            filedata = infile.read()
        self.attach_data(name, filename, filedata)

    def attach_key(self, key, value):
        ''' Attach a simple key=value pair '''
        attachment = MIMEText(str(value))
        attachment.add_header('Content-Disposition', 'form-data', 
                              name=key)
        self.attach(attachment)

    def flatten(self, debug = 0):
        ''' Flatten the mimedata with BytesGenerator and return bytes array '''
        if debug >= 2:
            print("Headers of mimedata as a dictionary:")
            print(dict(self.items()))
        bio = BytesIO()
        gen = FixedBytesGenerator(bio)
        gen.flatten(self, unixfrom=False)
        postdata = bio.getvalue() + b"\n"
        if debug >= 2:
            print("Postdata as a bytes array:")
            print(postdata)
        return postdata


class FileLockedException(IOError):
    """ Locking a file for exclusive access failed """
    pass

def open_exclusive(filename):
    """ Open a file and get an exlcusige lock on it """
    fileobj = open(filename)
    try:
        fcntl.flock(fileobj, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except IOError as err:
        if err.errno == errno.EACCES or err.errno == errno.EAGAIN:
            fileobj.close()
            raise FileLockedException(errno.EACCES, "File locked", filename)
        raise
    return fileobj

def close_exclusive(fileobj):
    """ Close a file, releasing any held lock on it """
    fcntl.flock(fileobj, fcntl.LOCK_UN)
    fileobj.close()

class WorkunitProcessor(object):
    def __init__(self, workunit, settings):
        self.settings = settings
        self.workunit = workunit
        self.errorcode = 0 # will get set if any command exits with code != 0
        self.failedcommand = None # If any command exits with code != 0, this
                                  # get set to the index of the failed command
        self.stdio = {"stdout": [], "stderr": []}

    def  __str__(self):
        return "Processor for Workunit:\n%s" % super(WorkunitProcessor, self)

    def renice(self):
        os.nice(int(self.settings["NICENESS"]))

    def run_commands(self):
        if self.result_exists():
            return True
        for (counter, command) in enumerate(self.workunit.get("COMMAND", [])):
            command = Template(command).safe_substitute(self.settings)
            logging.info ("Running command for workunit %s: %s", 
                          self.workunit.get_id(), command)

            # If niceness command line parameter was set, call self.renice() 
            # in child process, before executing command
            if int(self.settings["NICENESS"]) > 0:
                renice_func = self.renice
            else:
                renice_func = None

            # Run the command
            child = subprocess.Popen(command, shell=True, 
                                     stdout = subprocess.PIPE, 
                                     stderr = subprocess.PIPE, 
                                     preexec_fn = renice_func)
            # Wait for command to finish executing, capturing stdout and stderr 
            # in output tuple
            (child_stdout, child_stderr) = child.communicate()

            self.stdio["stdout"].append(child_stdout)
            self.stdio["stderr"].append(child_stderr)

            if child.returncode != 0:
                logging.error ("Command exited with exit code %s", 
                               child.returncode) 
                if self.stdio["stdout"][-1]:
                    logging.error ("Stdout: %s", self.stdio["stdout"][-1]) 
                if self.stdio["stderr"][-1]:
                    logging.error ("Stderr: %s", self.stdio["stderr"][-1]) 
                self.failedcommand = counter
                self.errorcode = child.returncode
                return False
            else:
                logging.debug ("Command exited successfully")
        return True

    def result_exists(self):
        ''' Check whether all result files already exist, returns True of False 
        '''
        for filename in self.workunit.get("RESULT", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            if not os.path.isfile(filepath):
                logging.info ("Result file %s does not exist", filepath)
                return False
            logging.info ("Result file %s already exists", filepath)
        logging.info ("All result files already exist")
        return True

    def cleanup(self):
        ''' Delete uploaded result files and files from DELETE lines '''
        logging.info ("Cleaning up for workunit %s", self.workunit.get_id())
        for filename in self.workunit.get("RESULT", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            logging.info ("Removing result file %s", filepath)
            os.remove(filepath)
        for filename in self.workunit.get("DELETE", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            logging.info ("Removing file %s", filepath)
            os.remove(filepath)

class WorkunitClient(object):
    def __init__(self, settings):
        self.settings = settings
        
        self.wu_filename = os.path.join(self.settings["DLDIR"], 
                                        self.settings["WU_FILENAME"])
        self.download_wu()

        # Get an exclusive lock to avoid two clients working on the same 
        # workunit
        try:
            self.wu_file = open_exclusive(self.wu_filename)
        except FileLockedException:
            logging.error("File '%s' is already locked. This may "
                          "indicate that two clients with clientid '%s' are "
                          "running. Terminating.", 
                          self.wu_filename, self.settings["CLIENTID"])
            raise

        logging.debug ("Parsing workunit from file %s", self.wu_filename)
        wu_text = self.wu_file.read()
        # WU file stays open so we keep the lock

        self.workunit = Workunit(wu_text)
        logging.debug ("Workunit ID is %s", self.workunit.get_id())
    
    def download_wu(self):
        # Download the WU file if none exists
        url = self.settings["GETWUPATH"]
        options = "clientid=" + self.settings["CLIENTID"]
        while not self.get_missing_file(url, self.wu_filename, options=options):
            logging.error("Error downloading workunit file")
            wait = float(self.settings["DOWNLOADRETRY"])
            logging.error("Waiting %s seconds before retrying", wait)
            time.sleep(wait)

    def cleanup(self):
        logging.info ("Removing workunit file %s", self.wu_filename)
        close_exclusive(self.wu_file)
        os.remove(self.wu_filename)

    @staticmethod
    def _urlopen(request, wait, is_upload = False):
        """ Wrapper around urllib2.urlopen() that retries in case of
        connection failure.
        """
        conn = None
        while conn is None:
            try:
                conn = urllib_request.urlopen(request)
            except (NameError, urllib_error.URLError) as error:
                conn = None
            except (socket.error) as error:
                if error.errno == errno.ECONNRESET:
                    conn = None
                else:
                    raise
            if not conn:
                logging.error("%s of result failed, %s", 
                              'Upload' if is_upload else 'Download', 
                              str(error))
                logging.error("Waiting %s seconds before retrying", wait)
                time.sleep(wait)
        return conn

    @staticmethod
    def get_content_charset(conn):
        """ Returns the character set of the server's response. 
        
        Defaults to latin-1 if no charset header was sent.
        The encoding may matter if the path names for the uploaded files 
        contain special characters, like accents.
        """
        if sys.version_info[0] == 3:
            return conn.info().get_content_charset()
        else:
            encoding = "latin-1" # Default value
            for item in conn.info().getplist():
                pair = item.split("=")
                if len(pair) == 2 and pair[0].strip() == "charset":
                    encoding = pair[1].strip()
        return encoding

    def get_file(self, urlpath, dlpath = None, options = None):
        # print('get_file("' + urlpath + '", "' + dlpath + '")')
        if dlpath == None:
            filename = urlpath.split("/")[-1]
            dlpath = os.path.join(self.settings["DLDIR"], filename)
        url = self.settings["SERVER"] + "/" + urlpath
        if options:
            url = url + "?" + options
        logging.info ("Downloading %s to %s", url, dlpath)
        wait = float(self.settings["DOWNLOADRETRY"])
        request = self._urlopen(url, wait)
        # Try to open the file exclusively
        try:
            fd = os.open(dlpath, os.O_CREAT | os.O_WRONLY | os.O_EXCL)
        except OSError as err:
            if err.errno == 17: # File exists error
                # There is a possible race condition here. If process A creates 
                # the file, then process B tries and finds that the file exists
                # and immediately get a shared lock for reading, then process A 
                # can never get an exclusive lock for writing.
                # To avoid this, we let process B wait until the file has 
                # positive size, which implies that process A must have the 
                # lock already. After 60 seconds, assume the file really has 0 
                # bytes and return
                logging.warning("Looks like another process already created "
                                "file %s", dlpath)
                slept = 0
                timeout = 60
                while slept < timeout and os.path.getsize(dlpath) == 0:
                    logging.warning("Sleeping until %s contains data", dlpath)
                    time.sleep(1)
                    slept += 1
                if slept == timeout:
                    logging.warning("Slept %d seconds, %s still has no data", 
                                    timeout, dlpath)
                return
            else:
                raise
        fcntl.flock(fd, fcntl.LOCK_EX)
        outfile = os.fdopen(fd, "wb")
        shutil.copyfileobj (request, outfile)
        fcntl.flock(fd, fcntl.LOCK_UN)
        outfile.close() # This should also close the fd
        request.close()
    
    def get_missing_file(self, urlpath, filename, checksum = None, 
                         options = None):
        # print('get_missing_file(%s, %s, %s)' % (urlpath, filename, checksum))
        if os.path.isfile(filename):
            logging.info ("%s already exists, not downloading", filename)
            if checksum is None:
                return True
            filesum = self.do_checksum(filename)
            if filesum.lower() == checksum.lower():
                return True
            logging.error ("Existing file %s has wrong checksum %s, "
                           "workunit specified %s. Deleting file.", 
                           filename, filesum, checksum)
            os.remove(filename)

        # If checksum is wrong and does not change during two downloads, exit 
        # with failue, as apparently the file on the server and checksum in 
        # workunit do not agree
        last_filesum = None
        while True:
            try:
                self.get_file(urlpath, filename, options)
            except urllib_error.HTTPError as error:
                logging.error (str(error))
                return False
            if checksum is None:
                return True
            filesum = self.do_checksum(filename)
            if filesum.lower() == checksum.lower():
                return True
            if not last_filesum is None and filesum == last_filesum:
                logging.error ("Downloaded file %s has same wrong checksum %s" 
                               " again. Exiting.", filename, filesum)
                return False
            logging.error ("Downloaded file %s has wrong checksum %s, " 
                           "workunit specified %s. Deleting file.", 
                           filename, filesum, checksum)
            os.remove(filename)
            last_filesum = filesum
        return True

    def get_files(self):
        for (filename, checksum) in self.workunit.get("FILE", []) + \
                self.workunit.get("EXECFILE", []):
            templ = Template(filename)
            archname = templ.safe_substitute({"ARCH": self.settings["ARCH"]})
            dlname = templ.safe_substitute({"ARCH": ""})
            dlpath = os.path.join(self.settings["DLDIR"], dlname)
            if not self.get_missing_file (archname, dlpath, checksum):
                return False
            # Try to lock the file once to be sure that download has finished
            # if another wuclient is doing the downloading
            with open(dlpath) as file_to_lock:
                fcntl.flock(file_to_lock, fcntl.LOCK_SH)
                fcntl.flock(file_to_lock, fcntl.LOCK_UN)
            
            if filename in dict(self.workunit.get("EXECFILE", [])):
                mode = os.stat(dlpath).st_mode
                if mode & stat.S_IXUSR == 0:
                    logging.info ("Setting executable flag for %s", dlpath)
                    os.chmod(dlpath, mode | stat.S_IXUSR)
        return True

    def attach_result(self, processor, mimedata):
        # Build a multi-part MIME document containing the WU id and result file
        mimedata.attach_key("WUid", self.workunit.get_id())
        mimedata.attach_key("clientid", self.settings["CLIENTID"])
        if processor.errorcode:
            mimedata.attach_key("errorcode", processor.errorcode)
        if processor.failedcommand:
            mimedata.attach_key("failedcommand", processor.failedcommand)
        for filename in self.workunit.get("RESULT", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            mimedata.attach_file("results", filename, filepath)
        for name in processor.stdio:
            for (counter, data) in enumerate(processor.stdio[name]):
                if data:
                    logging.debug ("Adding %s for command %s to upload", 
                                   name, counter)
                    mimedata.attach_data("results", name + str(counter), data)
        
        return mimedata

    def upload_result(self, processor):
        # Make POST data
        mimedata = WuMIMEMultipart()
        self.attach_result(processor, mimedata)
        postdata = mimedata.flatten(debug=int(self.settings["DEBUG"]))
        logging.debug("POST data: %s", mimedata)

        url = self.settings["SERVER"].rstrip("/") + "/" + \
                self.settings["POSTRESULTPATH"]
        logging.info("Sending result for workunit %s to %s", 
                     self.workunit.get_id(), url)
        request = urllib_request.Request(url, data=postdata, 
                                         headers=dict(mimedata.items()))
        wait = float(self.settings["DOWNLOADRETRY"])
        conn = self._urlopen(request, wait, is_upload=True)
        if not conn:
            return False
        response = conn.read()

        encoding = self.get_content_charset(conn)
        if sys.version_info[0] == 2:
            response_str = unicode(response, encoding=encoding)
        else:
            response_str = response.decode(encoding=encoding)
        logging.debug ("Server response:\n%s", response_str)
        conn.close()
        return True

    @staticmethod
    def do_checksum(filename, checksum = None):
        """ Computes the SHA1 checksum for a file. If checksum is None, returns 
            the computed checksum. If checksum is not None, return whether the
            computed SHA1 sum and checksum agree """
        blocksize = 65536
        sha1hash = hashlib.sha1() # pylint: disable=E1101
        infile = open(filename, "rb")
        fcntl.flock(infile, fcntl.LOCK_SH)
        
        data = infile.read(blocksize)
        while data:
            sha1hash.update(data)
            data = infile.read(blocksize)
        fcntl.flock(infile, fcntl.LOCK_UN)
        infile.close()
        filesum = sha1hash.hexdigest()
        if checksum is None:
            return filesum
        else:
            return filesum.lower() == checksum.lower()

    def have_terminate_request(self):
        return not self.workunit.get("TERMINATE", None) is None

    def process(self):
        # If all output files exist, send them, return WU as done
        # Otherwise, run commands in WU. If no error and all output 
        #   files exist, send them, return WU as done
        # print(wu)
        if self.have_terminate_request():
            logging.info ("Received TERMINATE, exiting")
            return False

        if not self.get_files():
            return False
        processor = WorkunitProcessor(self.workunit, self.settings)
        processor.run_commands()
        upload_ok = self.upload_result(processor)
        processor.cleanup()
        self.cleanup()
        return upload_ok


# Settings which we require on the command line (no defaults)
REQUIRED_SETTINGS = {"CLIENTID" : (None, "Unique ID for this client"), 
                     "SERVER" : (None, "Base URL for WU server")}

# Optional settings with defaults, overrideable on command line, 
# and a help text
OPTIONAL_SETTINGS = {"WU_FILENAME" : 
                     (None, "Filename under which to store WU files"), 
                     "DLDIR" : ('download/', "Directory for downloading files"),
                     "WORKDIR" : (None, "Directory for result files"),
                     "BASEPATH" : (None, "Base directory for download and work "
                                         "directories"),
                     "GETWUPATH" : 
                     ("/cgi-bin/getwu", 
                      "Path segment of URL for requesting WUs from server"), 
                     "POSTRESULTPATH" : 
                     ("/cgi-bin/upload.py", 
                      "Path segment of URL for reporting results to server"), 
                     "DEBUG" : ("0", "Debugging verbosity"),
                     "ARCH" : ("", "Architecture string for this client"),
                     "DOWNLOADRETRY" : 
                     ("300", "Time to wait before download retries"),
                     "NICENESS" : 
                     ("0", "Run subprocesses under this niceness"),
                     "LOGLEVEL" : ("INFO", "Verbosity of logging"),
                     "LOGFILE" : (None, "File to which to write log output")
                     }
# Merge the two, removing help string
SETTINGS = dict([(a, b) for (a, (b, c)) in list(REQUIRED_SETTINGS.items()) + \
                                        list(OPTIONAL_SETTINGS.items())])

if __name__ == '__main__':

    def parse_cmdline():
        # Create command line parser from the keys in SETTINGS
        parser = optparse.OptionParser()
        for (arg, default) in REQUIRED_SETTINGS.items():
            parser.add_option('--' + arg.lower(), help=default[1])
        for (arg, default) in OPTIONAL_SETTINGS.items():
            if not default[0] is None:
                parser.add_option('--' + arg.lower(), default=default[0], 
                    help=default[1] + " (default: " + default[0] + ")")
            else:
                parser.add_option('--' + arg.lower(), help=default[1])
        # Parse command line
        (options, args) = parser.parse_args()
        if args:
            sys.stderr.write("Did not understand command line arguments %s\n" %
                             " ".join(args))
            raise Exception()
        # Copy values to SETTINGS
        for arg in SETTINGS.keys():
            if hasattr(options, arg.lower()):
                SETTINGS[arg] = getattr(options, arg.lower())
        for arg in REQUIRED_SETTINGS.keys():
            if SETTINGS[arg] is None:
                raise Exception("Command line parameter --%s is required" 
                                % arg.lower())

    parse_cmdline()
    # If no working directory is given, we use <clientid>.work/
    if SETTINGS["WORKDIR"] is None:
        SETTINGS["WORKDIR"] = SETTINGS["CLIENTID"] + '.work/'
    if not SETTINGS["BASEPATH"] is None:
        SETTINGS["WORKDIR"] = os.path.join(SETTINGS["BASEPATH"], 
                                           SETTINGS["WORKDIR"])
        SETTINGS["DLDIR"] = os.path.join(SETTINGS["BASEPATH"], 
                                         SETTINGS["DLDIR"])
    # If no WU filename is given, we use "WU." + client id
    if SETTINGS["WU_FILENAME"] is None:
        SETTINGS["WU_FILENAME"] = "WU." + SETTINGS["CLIENTID"]

    # Create download and working directories if they don't exist
    if not os.path.isdir(SETTINGS["DLDIR"]):
        os.mkdir(SETTINGS["DLDIR"])
    if not os.path.isdir(SETTINGS["WORKDIR"]):
        os.mkdir(SETTINGS["WORKDIR"])

    loglevel = getattr(logging, SETTINGS["LOGLEVEL"].upper(), None)
    if not isinstance(loglevel, int):
        raise ValueError('Invalid log level: ' + SETTINGS["LOGLEVEL"])
    logfile = SETTINGS["LOGFILE"]
    logging.basicConfig(level=loglevel, filename=logfile)

    # print (str(SETTINGS))

    client_ok = True
    while client_ok:
        client = WorkunitClient(settings = SETTINGS)
        client_ok = client.process()
