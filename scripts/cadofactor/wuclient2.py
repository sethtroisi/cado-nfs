#!/usr/bin/env python2

import sys
import os
import fcntl
import errno
import stat
import optparse
import shutil
import time
try:
    import urllib2 as urllib_request
    import urllib2 as urllib_error
except ImportError:
    import urllib.request as urllib_request
    import urllib.error as urllib_error
import subprocess
import hashlib
import logging
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import email.encoders
import email.generator
from string import Template
from io import BytesIO
from workunit import Workunit

if sys.version_info[0] == 3 and sys.version_info[1] < 3:
    # In Python 3.[012], use a fixed BytesGenerator which accepts a bytes 
    # input. The fact that the BytesGenerator in Python 3.[012] doesn't 
    # is a bug, see http://bugs.python.org/issue16564
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
else:
    # In Python >=3.3, use the bug-fixed bytes generator
    # pylint: disable=E1101
    # pylint: disable=C0103
    FixedBytesGenerator = email.generator.BytesGenerator


class WuMIMEMultipart(MIMEMultipart):
    ''' Defines convenience functions for attaching files and data '''
    
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


class WorkunitProcessor(Workunit):
    def __init__(self, text, settings):
        self.exitcode = 0 # will get set if any command exits with code != 0
        self.failedcommand = None
        self.settings = settings
        super(WorkunitProcessor, self).__init__(text)
        logging.debug ("Workunit ID is %s", self.get_id())
        self.stdio = {"stdout": [], "stderr": []}

    def  __str__(self):
        return "Processor for Workunit:\n%s" % super(WorkunitProcessor, self)

    def have_terminate_request(self):
        return "TERMINATE" in self.wudata

    def renice(self):
        os.nice(int(self.settings["NICENESS"]))

    @staticmethod
    def paste_path(*arr):
        arr2 = [p.rstrip(os.sep) for p in arr[:-1]] + [arr[-1]]
        return os.sep.join(arr2)
    
    def run_commands(self):
        for (counter, command) in enumerate(self.wudata.get("COMMAND", [])):
            command = Template(command).safe_substitute(self.settings)
            logging.info ("Running command for workunit %s: %s", 
                          self.get_id(), command)

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
                self.failedcommand = counter
                self.exitcode = child.returncode
                return False
            else:
                logging.debug ("Command exited successfully")
        return True

    def result_exists(self):
        ''' Check whether all result files already exist, returns True of False 
        '''
        if "RESULT" in self.wudata:
            for filename in self.wudata["RESULT"]:
                filepath = self.paste_path(self.settings["WORKDIR"], filename)
                if not os.path.isfile(filepath):
                    logging.info ("Result file %s does not exist", filepath)
                    return False
                logging.info ("Result file %s already exists", filepath)
        logging.info ("All result files already exist")
        return True

    def cleanup(self):
        ''' Delete uploaded result files and files from DELETE lines '''
        logging.info ("Cleaning up for workunit %s", self.get_id())
        if "RESULT" in self.wudata:
            for filename in self.wudata["RESULT"]:
                filepath = self.paste_path(self.settings["WORKDIR"], filename)
                logging.info ("Removing result file %s", filepath)
                os.remove(filepath)
        if "DELETE" in self.wudata:
            for filename in self.wudata["DELETE"]:
                filepath = self.paste_path(self.settings["WORKDIR"], filename)
                logging.info ("Removing file %s", filepath)
                os.remove(filepath)

class WorkunitProcessorClient(WorkunitProcessor):
    def __init__(self, settings):
        self.settings = settings
        # Download the WU file if none exists
        url = self.settings["GETWUPATH"]
        self.wu_filename = self.paste_path(self.settings["DLDIR"], self.settings["WU_FILENAME"])
        options = "clientid=" + self.settings["CLIENTID"]
        while not self.get_missing_file(url, self.wu_filename, options=options):
            logging.error("Error downloading workunit file")
            wait = float(self.settings["DOWNLOADRETRY"])
            logging.error("Waiting %s seconds before retrying", wait)
            time.sleep(wait)

        # Parse the contents of the WU file
        logging.debug ("Parsing workunit from file %s", self.wu_filename)

        self.wu_file = open(self.wu_filename)

        # Get an exclusive lock to avoid two clients working on the same 
        # workunit
        try:
            fcntl.flock(self.wu_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except IOError as err:
            if err.errno == errno.EACCES or err.errno == errno.EAGAIN:
                logging.error("Workunit file is already locked. This may "
                              "indicate that two clients with clientid %s are "
                              "running. Terminating.", self.settings["CLIENTID"])
            raise

        wu_text = self.wu_file.read()
        super(WorkunitProcessorClient, self).__init__(
            text = wu_text, settings = settings)
    
    def cleanup(self):
        super(WorkunitProcessorClient, self).cleanup()
        logging.info ("Removing workunit file %s", self.wu_filename)
        
        fcntl.flock(self.wu_file, fcntl.LOCK_UN)
        self.wu_file.close()
        os.remove(self.wu_filename)

    def _urlopen(self, request, operation):
        """ Wrapper around urllib2.urlopen() that retries in case of
        connection failure.
        """
        conn = None
        while conn is None:
            try:
                conn = urllib_request.urlopen(request)
            except (NameError, urllib_error.URLError) as error:
                conn = None
                wait = float(self.settings["DOWNLOADRETRY"])
                logging.error("%s of result failed, %s", operation, 
                              str(error))
                logging.error("Waiting %s seconds before retrying", wait)
                time.sleep(wait)
            except Exception:
                raise
        return conn

    def get_file(self, urlpath, dlpath = None, options = None):
        # print('get_file("' + urlpath + '", "' + dlpath + '")')
        if dlpath == None:
            dlpath = self.paste_path(self.settings["DLDIR"], urlpath.split("/")[-1])
        url = self.settings["SERVER"] + "/" + urlpath
        if options:
            url = url + "?" + options
        logging.info ("Downloading %s to %s", url, dlpath)
        request = self._urlopen(url, "Download")
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
        outfile = os.fdopen(fd, "w")
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
        for (filename, checksum) in self.wudata.get("FILE", []) + \
                self.wudata.get("EXECFILE", []):
            templ = Template(filename)
            archname = templ.safe_substitute({"ARCH": self.settings["ARCH"]})
            dlname = templ.safe_substitute({"ARCH": ""})
            dlpath = self.paste_path(self.settings["DLDIR"], dlname)
            if not self.get_missing_file (archname, dlpath, checksum):
                return False
            # Try to lock the file once to be sure that download has finished
            # if another wuclient is doing the downloading
            with open(dlpath) as f:
                fcntl.flock(f, fcntl.LOCK_SH)
                fcntl.flock(f, fcntl.LOCK_UN)
            
            if filename in dict(self.wudata.get("EXECFILE", [])):
                mode = os.stat(dlpath).st_mode
                if mode & stat.S_IXUSR == 0:
                    logging.info ("Setting executable flag for %s", dlpath)
                    os.chmod(dlpath, mode | stat.S_IXUSR)
        return True

    def _flatten(self, postdata):
        ''' Flatten the postdata with BytesGenerator and return bytes array '''
        if int(self.settings["DEBUG"]) >= 2:
            print("Headers of postdata as a dictionary:")
            print(dict(postdata.items()))
        bio = BytesIO()
        gen = FixedBytesGenerator(bio)
        gen.flatten(postdata, unixfrom=False)
        postdata = bio.getvalue() + b"\n"
        if int(self.settings["DEBUG"]) >= 2:
            print("Postdata as a bytes array:")
            print(postdata)
        return postdata

    def _make_mimedata(self):
        # Build a multi-part MIME document containing the WU id and result file
        mimedata = WuMIMEMultipart()
        mimedata.attach_key("WUid", self.get_id())
        mimedata.attach_key("clientid", self.settings["CLIENTID"])
        if self.exitcode:
            mimedata.attach_key("exitcode", self.exitcode)
        if self.failedcommand:
            mimedata.attach_key("failedcommand", self.failedcommand)
        if "RESULT" in self.wudata:
            for filename in self.wudata["RESULT"]:
                filepath = self.paste_path(self.settings["WORKDIR"], filename)
                mimedata.attach_file("results", filename, filepath)
        for name in self.stdio:
            for (counter, data) in enumerate(self.stdio[name]):
                if data:
                    logging.debug ("Adding %s for command %s to upload", 
                                   name, counter)
                    mimedata.attach_data("results", name + str(counter), data)
        
        return mimedata

    def upload_result(self):
        mimedata = self._make_mimedata()

        url = self.settings["SERVER"] + "/" + self.settings["POSTRESULTPATH"]
        logging.info("Sending result for workunit %s to %s", 
                     self.get_id(), url)
        request = urllib_request.Request(url, data=self._flatten(mimedata), 
                                         headers=dict(mimedata.items()))
        conn = self._urlopen(request, "Upload")
        if not conn:
            return False
        response = conn.read()

        encoding = "latin-1" # Default value
        # Find out the encoding the server response uses. This may matter if 
        # the path names for the uploaded files contain special characters,
        # like accents
        if sys.version_info[0] == 3:
            encoding = conn.info().get_content_charset()
        else:
            for item in conn.info().getplist():
                pair = item.split("=")
                if len(pair) == 2 and pair[0].strip() == "charset":
                    encoding = pair[1].strip()

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
        if not self.result_exists():
            if not self.run_commands():
                return False
        if not self.upload_result():
            return False
        self.cleanup()
        return True


# Settings which we require on the command line (no defaults)
REQUIRED_SETTINGS = {"CLIENTID" : (None, "Unique ID for this client"), 
                     "SERVER" : (None, "Base URL for WU server")}

# Optional settings with defaults, overrideable on command line, 
# and a help text
OPTIONAL_SETTINGS = {"WU_FILENAME" : 
                     (None, "Filename under which to store WU files"), 
                     "DLDIR" : ('download/', "Directory for downloading files"), 
                     "WORKDIR" : (None, "Directory for result files"),
                     "BASEPATH" : (None, "Base directory for download and work directories"),
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
        SETTINGS["WORKDIR"] = os.path.join(SETTINGS["BASEPATH"], SETTINGS["WORKDIR"])
        SETTINGS["DLDIR"] = os.path.join(SETTINGS["BASEPATH"], SETTINGS["DLDIR"])
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

    ok = True
    while ok:
        processor = WorkunitProcessorClient(settings = SETTINGS)
        ok = processor.process()
