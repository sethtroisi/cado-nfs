#!/usr/bin/env python2

# TODO: file locking for downloading files, so several clients can use the same files (factorbase etc.) safely
# TODO: file locking for workunit files, so name collision can be detected

import sys
import os
import stat
import optparse
import shutil
import time
import urllib2
import subprocess
import hashlib
import logging
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import email.encoders
import email.generator
from string import Template
from io import StringIO,BytesIO,open
from workunit import Workunit

if sys.version_info[0] == 2:
    # In Python 2.x, use the regular email generator
    class FixedBytesGenerator(email.generator.Generator):
        pass
elif sys.version_info[0] == 3 and sys.version_info[1] < 3:
    # In Python 3.[012], use a fixed BytesGenerator which accepts a bytes 
    # input. The fact that the BytesGenerator in Python 3.[012] does not is 
    # a bug, see http://bugs.python.org/issue16564
    class FixedBytesGenerator(email.generator.BytesGenerator):
        def _handle_bytes(self, msg):
            payload = msg.get_payload()
            if payload is None:
                return
            if isinstance(payload, bytes):
                # Payload is bytes, output is bytes - just write them
                self._fp.write(payload)
            elif isinstance(payload, str):
                super(FixedBytesGenerator,self)._handle_text(msg)
            else:
                # Payload is neither bytes nor string - this can't be right
                raise TypeError('bytes payload expected: %s' % type(payload))
        _writeBody = _handle_bytes


# Settings which we require on the command line (no defaults)
REQUIRED_SETTINGS = {"CLIENTID" : (None, "Unique ID for this client"), 
                     "DLDIR" : (None, "Directory for downloading files"), 
                     "SERVER" : (None, "Base URL for WU server"), 
                     "WORKDIR" : (None, "Directory for result files")}
# Optional settings with defaults, overrideable on command line, and a help text
OPTIONAL_SETTINGS = {"WU_FILENAME" : (None, "Filename under which to store WU files"), 
                     "GETWUPATH" : ("/cgi-bin/getwu", "Path segment of URL for requesting WUs from server"), 
                     "POSTRESULTPATH" : ("/cgi-bin/upload.py", "Path segment of URL for reporting results to server"), 
                     "DEBUG" : ("0", "Debugging verbosity"),
                     "ARCH" : ("", "Architecture string for this client"),
                     "DOWNLOADRETRY" : ("300", "Time to wait before download retries"),
                     "NICENESS" : ("0", "Run subprocesses under this niceness"),
                     "LOGLEVEL" : ("INFO", "Verbosity of logging"),
                     "LOGFILE" : (None, "File to which to write log output")}
# Merge the two, removing help string
SETTINGS = dict([(a,b) for (a,(b,c)) in list(REQUIRED_SETTINGS.items()) + \
                                        list(OPTIONAL_SETTINGS.items())])

def get_file(urlpath, dlpath = None, options = None):
    # print('get_file("' + urlpath + '", "' + dlpath + '")')
    if dlpath == None:
        dlpath = SETTINGS["DLDIR"] + os.sep + urlpath.split("/")[-1]
    url = SETTINGS["SERVER"] + "/" + urlpath
    if options:
        url = url + "?" + options
    logging.info ("Downloading " + url + " to " + dlpath);
    request = None
    while request == None:
        try:
            request = urllib2.urlopen(url)
        except urllib2.HTTPError, e:
            logging.error (unicode(e))
            return False
        except (NameError, urllib2.URLError), e:
            request = None
            wait = float(SETTINGS["DOWNLOADRETRY"])
            logging.error("Download of " + urlpath + " failed, " + unicode(e))
            logging.error("Waiting " + str(wait) + " seconds before retrying")
            time.sleep(wait)
    file = open(dlpath, "wb")
    shutil.copyfileobj (request, file)
    file.close()
    request.close()
    return True

def do_checksum(filename, checksum = None):
    """ Computes the SHA1 checksum for a file. If checksum is None, returns 
        the computed checksum. If checksum is not None, return whether the
        computed SHA1 sum and checksum agree """
    blocksize = 65536
    m = hashlib.sha1()
    file = open(filename, "rb")
    data = file.read(blocksize)
    while data:
        m.update(data)
        data = file.read(blocksize)
    file.close()
    filesum = m.hexdigest()
    if checksum is None:
        return filesum
    else:
        return filesum.lower() == checksum.lower()

def get_missing_file(urlpath, filename, checksum = None, options = None):
    # print('get_missing_file("' + urlpath + '", "' + filename + '", ' + str(checksum) + ')')
    if os.path.isfile(filename):
        logging.info (filename + " already exists, not downloading")
        if checksum is None:
            return True
        filesum = do_checksum(filename)
        if filesum.lower() == checksum.lower():
            return True
        logging.error ("Existing file " + filename + " has wrong checksum " + filesum + 
             ", workunit specified " + checksum +". Deleting file.")
        os.remove(filename)
    
    # If checksum is wrong and does not change during two downloads, exit with 
    # failue, as apparently the file on the server and checksum in workunit do 
    # not agree
    last_filesum = None
    while True:
        if not get_file(urlpath, filename, options):
            return False
        if checksum is None:
            return True
        filesum = do_checksum(filename)
        if filesum.lower() == checksum.lower():
            return True
        if not last_filesum is None and filesum == last_filesum:
            logging.error ("Downloaded file " + filename + " has same wrong checksum " 
                 + filesum + " again. Exiting.")
            return False
        logging.error ("Downloaded file " + filename + " has wrong checksum " + 
             filesum + ", workunit specified " + checksum + ". Deleting file.")
        os.remove(filename)
        last_filesum = filesum

class WorkunitProcessor(Workunit):
    
    def __init__(self, filepath, debug = 0):
        self.exitcode = 0 # will get set if any command exits with code != 0
        self.failedcommand = None
        self.clientid = SETTINGS["CLIENTID"]
        self.debug = debug # Controls debugging output

        logging.debug ("Parsing workunit from file " + filepath)
        wu_file = open(filepath)
        wu_text = wu_file.read()
        wu_file.close()
        super(WorkunitProcessor,self).__init__(wu_text)
        self.WUid = self.get_id()
        logging.debug (" done, workunit ID is " + self.get_id())
        self.stdout = []
        self.stderr = []

    def  __str__(self):
        return "Processor for Workunit:\n" + super(WorkunitProcessor,self).__str__()

    def have_terminate_request(self):
        return "TERMINATE" in self.wudata

    def get_files(self):
        for (filename, checksum) in self.wudata.get("FILE", []) + self.wudata.get("EXECFILE", []):
            archname = Template(filename).safe_substitute({"ARCH": SETTINGS["ARCH"]})
            dlname = Template(filename).safe_substitute({"ARCH": ""})
            if not get_missing_file (archname, SETTINGS["DLDIR"] + os.sep + dlname, checksum):
                return False
        for (filename, checksum) in self.wudata.get("EXECFILE", []):
            dlname = Template(filename).safe_substitute({"ARCH": ""})
            path = SETTINGS["DLDIR"] + os.sep + dlname
            mode = os.stat(path).st_mode
            if mode & stat.S_IXUSR == 0:
                logging.info ("Setting executable flag for " + path)
                os.chmod(path, mode | stat.S_IXUSR)
        return True

    @staticmethod
    def renice():
        os.nice(int(SETTINGS["NICENESS"]))
    
    def run_commands(self):
        for (counter, command) in enumerate(self.wudata.get("COMMAND", [])):
            command = Template(command).safe_substitute(SETTINGS)
            logging.info ("Running command for workunit " + self.get_id() + ": " + command)

            # If niceness command line parameter was set, call self.renice() in
            # child process, before executing command
            if int(SETTINGS["NICENESS"]) > 0:
                renice_func = self.renice
            else:
                renice_func = None

            # Run the command
            child = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, 
                stderr = subprocess.PIPE, preexec_fn = renice_func)
            # Wait for command to finish executing, capturing stdout and stderr 
            # in output tuple
            (child_stdout, child_stderr) = child.communicate()

            if len(child_stdout) > 0:
                self.stdout.append(child_stdout)
            else:
                self.stdout.append(None)
            if len(child_stderr) > 0:
                self.stderr.append(child_stderr)
            else:
                self.stderr.append(None)

            if child.returncode != 0:
                logging.error ("Command exited with exit code " + str(child.returncode)) 
                self.failedcommand = counter
                self.exitcode = child.returncode
                return False
            else:
                logging.debug ("Command exited successfully")
        return True

    def upload_result(self):
        # Build a multi-part MIME document containing the WU id and result file
        postdata = MIMEMultipart()
        for key in ("WUid", "clientid", "exitcode", "failedcommand"):
            if not getattr(self, key, None) is None:
                attachment = MIMEText(unicode(getattr(self, key)))
                attachment.add_header('Content-Disposition', 'form-data', name=key)
                postdata.attach(attachment)
        if "RESULT" in self.wudata:
            for f in self.wudata["RESULT"]:
                filepath = SETTINGS["WORKDIR"] + os.sep + f
                logging.debug ("Adding result file " + filepath + " to upload")
                file = open(filepath, "rb")
                filedata = file.read()
                file.close()
                # HTTP does not use a Content-Transfer-Encoding, so use noop encoder
                # This fails, probably related to http://bugs.python.org/issue4768
                result = MIMEApplication(filedata, _encoder=email.encoders.encode_noop)
                result.add_header('Content-Disposition', 'form-data', name="results", 
                                  filename=f)
                postdata.attach(result)
        for (name, arr) in (("stdout", self.stdout), ("stderr", self.stderr)):
            for (counter, f) in enumerate(arr):
                if not f is None:
                    logging.debug ("Adding " + name + " for command " + str(counter) + " to upload")
                    result = MIMEApplication(f, _encoder=email.encoders.encode_noop)
                    result.add_header('Content-Disposition', 'form-data', name="results", 
                                      filename=name + str(counter))
                    postdata.attach(result)
        
        if self.debug >= 2:
            print("Headers of postdata as a dictionary:")
            print(dict(postdata.items()))
        fp = BytesIO()
        g = FixedBytesGenerator(fp)
        g.flatten(postdata, unixfrom=False)
        postdata3 = fp.getvalue() + b"\n"
        if self.debug >= 2:
            print("Postdata as a bytes array:")
            print(postdata3)

        url = SETTINGS["SERVER"] + "/" + SETTINGS["POSTRESULTPATH"]
        logging.info("Sending result for workunit " + self.get_id() + " to " + url)
        request = urllib2.Request(url, data=postdata3, headers=dict(postdata.items()))
        conn = None;
        while conn is None:
            try:
                conn = urllib2.urlopen(request)
            except (urllib2.URLError), e:
                conn = None
                wait = float(SETTINGS["DOWNLOADRETRY"])
                logging.error("Upload of result failed, " + unicode(e))
                logging.error("Waiting " + str(wait) + " seconds before retrying")
                time.sleep(wait)
        response = conn.read()
        encoding = None
        # Find out the encoding the server response uses. This may matter if 
        # the path names for the uploaded files contain special characters,
        # like accents
        content_type = conn.info().get("Content-Type", default=None)
        if not content_type is None:
            # If there are multiple header lines with the same key, their 
            # values are joind with "," separators
            for h in content_type.split(","):
                for g in h.split(";"):
                    f = g.split("=")
                    if len(f) == 2 and f[0].strip() == "charset":
                        encoding = f[1].strip()
        if encoding is None:
            encoding = "latin-1"
        logging.debug ("Server response:\n" + unicode(response, encoding=encoding))
        conn.close()
        return True

    def result_exists(self):
        if "RESULT" in self.wudata:
            for f in self.wudata["RESULT"]:
                filepath = SETTINGS["WORKDIR"] + os.sep + f
                if not os.path.isfile(filepath):
                    logging.info ("Result file " + filepath + " does not exist")
                    return False
                logging.info ("Result file " + filepath + " already exists")
        logging.info ("All result files already exist")
        return True

    def cleanup(self):
        logging.info ("Cleaning up for workunit " + self.get_id())
        if "RESULT" in self.wudata:
            for f in self.wudata["RESULT"]:
                filepath = SETTINGS["WORKDIR"] + os.sep + f
                logging.info ("Removing result file " + filepath)
                os.remove(filepath)
        if "DELETE" in self.wudata:
            for f in self.wudata["DELETE"]:
                filepath = SETTINGS["WORKDIR"] + os.sep + f
                logging.info ("Removing file " + filepath)
                os.remove(filepath)

    def process(self):
        # If all output files exist, send them, return WU as done
        # Otherwise, run commands in WU. If no error and all output 
        #   files exist, send them, return WU as done
        # print(wu)
        if not self.get_files():
            return False
        if not self.result_exists():
            if not self.run_commands():
                return False
        if not self.upload_result():
            return False
        self.cleanup()
        return True

def do_work():
    rc = True
    wu_filename = SETTINGS["DLDIR"] + os.sep + SETTINGS["WU_FILENAME"]
    if not get_missing_file(SETTINGS["GETWUPATH"], wu_filename, 
                            options="clientid=" + SETTINGS["CLIENTID"]):
        return False
    wu = WorkunitProcessor(wu_filename, int(SETTINGS["DEBUG"]))
    if wu.have_terminate_request():
        logging.info ("Received TERMINATE, exiting")
        rc = False
    elif not wu.process():
        return False
    logging.info ("Removing workunit file " + wu_filename)
    os.remove(wu_filename)
    return rc

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
                raise Exception("Command line parameter --" + arg.lower() + " is required")

    parse_cmdline()
    # If no WU filename is given, we use "WU." + client id
    if SETTINGS["WU_FILENAME"] is None:
        SETTINGS["WU_FILENAME"] = "WU." + SETTINGS["CLIENTID"]

    logopt = {}
    logopt["level"] = getattr(logging, SETTINGS["LOGLEVEL"].upper(), None)
    if not isinstance(logopt["level"], int):
        raise ValueError('Invalid log level: ' + SETTINGS["LOGLEVEL"])
    if not SETTINGS["LOGFILE"] is None:
        logopt["filename"] = SETTINGS["LOGFILE"]
    logging.basicConfig(**logopt)

    # print (str(SETTINGS))

    while do_work():
        pass

