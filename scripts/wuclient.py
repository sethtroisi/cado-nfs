#!/usr/bin/env python3

import sys
import os
import stat
import argparse
import shutil
import time
import urllib.request
import subprocess
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import email.encoders
import email.generator
from string import Template
from io import StringIO,BytesIO
from workunit import Workunit

class FixedBytesGenerator(email.generator.BytesGenerator):
    def _handle_bytes(self, msg):
        payload = msg.get_payload()
        if payload is None:
            return
        if isinstance(payload, bytes):
            # Payload is bytes, output is bytes - just write them
            self._fp.write(payload)
        else:
            # Payload is neither bytes not string - this can't be right
            raise TypeError('bytes payload expected: %s' % type(payload))
    _writeBody = _handle_bytes


# Settings which we require on the command line (no defaults)
REQUIRED_SETTINGS = {"CLIENTID" : ("", "Unique ID for this client"), 
                     "DLDIR" : ("", "Directory for downloading files"), 
                     "SERVER" : ("", "Base URL for WU server"), 
                     "WORKDIR" : ("", "Directory for result files")}
# Optional settings with defaults, overrideable on command line, and a help text
OPTIONAL_SETTINGS = {"WU_FILENAME" : ("WU", "Filename under which to store WU files"), 
                     "GETWUPATH" : ("/cgi-bin/getwu", "Path segment of URL for requesting WUs from server"), 
                     "POSTRESULTPATH" : ("/cgi-bin/upload.py", "Path segment of URL for reporting results to server"), 
                     "DEBUG" : ("0", "Debugging verbosity")}
# Merge the two, removing help string
SETTINGS = dict([(a,b) for (a,(b,c)) in list(REQUIRED_SETTINGS.items()) + \
                                        list(OPTIONAL_SETTINGS.items())])

def get_file(urlpath, dlpath = None, options = None):
    # print('get_file("' + urlpath + '", "' + dlpath + '")')
    if dlpath == None:
        dlpath = SETTINGS["DLDIR"] + "/" + os.basename(urlpath) # FIXME: should be url base name
    url = SETTINGS["SERVER"] + "/" + urlpath
    if options:
        url = url + "?" + options
    print ("Downloading " + url + " to " + dlpath);
    request = urllib.request.urlopen(url)
    file = open(dlpath, "wb")
    shutil.copyfileobj (request, file)
    file.close()
    request.close()

def get_missing_file(urlpath, filename, checksum = None, options = None):
    # print('get_missing_file("' + urlpath + '", "' + filename + '", ' + str(checksum) + ')')
    if not os.path.isfile(filename):
        get_file(urlpath, filename, options)
        # FIXME CHECKSUM
    else:
        print (filename + " already exists, not downloading")
    if checksum:
        # FIXME CHECKSUM
        pass
    return True

class Workunit_Processor(Workunit):
    
    def __init__(self, filepath, debug = 0):
        self.exitcode = 0 # will get set if any command exits with code != 0
        self.debug = debug # Controls debugging output

        if self.debug >= 1:
            print ("Parsing workunit from file " + filepath)
        wu_file = open(filepath)
        wu_text = wu_file.read()
        wu_file.close()
        self.wu = Workunit(wu_text)
        if self.debug >= 1:
            print (" done, workunit ID is " + self.wu.get_id())

    def  __str__(self):
        return "Processor for Workunit:\n" + self.wu.__str__()

    def get_files(self):
        for (filename, checksum) in self.wu.data["FILE"] + self.wu.data["EXECFILE"]:
            if not get_missing_file (filename, SETTINGS["DLDIR"] + '/' + filename, checksum):
                return False
        for (filename, checksum) in self.wu.data["EXECFILE"]:
            path = SETTINGS["DLDIR"] + '/' + filename
            mode = os.stat(path).st_mode
            if mode & stat.S_IXUSR == 0:
                print ("Setting executable flag for " + path)
                os.chmod(path, mode | stat.S_IXUSR)
        return True
    
    def run_commands(self):
        for command in self.wu.data["COMMAND"]:
            command = Template(command).safe_substitute(SETTINGS)
            if self.debug >= 0:
                print ("Running command for workunit " + self.wu.get_id() + ": " + command)
            rc = subprocess.call(command, shell=True)
            if rc != 0:
                print ("Command exited with exit code " + str(rc)) 
                self.exitcode = rc
                return False
            elif self.debug >= 1:
                print ("Command exited successfully")
        return True

    def upload_result(self):
        # Build a multi-part MIME document containing the WU id and result file
        postdata = MIMEMultipart()
        WUid = MIMEText(self.wu.get_id())
        WUid.add_header('Content-Disposition', 'form-data', name="WUid")
        postdata.attach(WUid)
        clientid = MIMEText(SETTINGS["CLIENTID"])
        clientid.add_header('Content-Disposition', 'form-data', name="clientid")
        postdata.attach(clientid)
        if self.exitcode > 0:
            rc = MIMEText(str(self.exitcode))
            WUid.add_header('Content-Disposition', 'form-data', name="exitcode")
            postdata.attach(rc)
        if "RESULT" in self.wu.data:
            filepath = SETTINGS["WORKDIR"] + "/" + self.wu.data["RESULT"]
            if self.debug >= 1:
                print ("Adding result file " + filepath + " to upload")
            file = open(filepath, "rb")
            filedata = file.read()
            file.close()
            # HTTP does not use a Content-Transfer-Encoding, so use noop encoder
            # This fails, probably related to http://bugs.python.org/issue4768
            result = MIMEApplication(filedata, _encoder=email.encoders.encode_noop)
            result.add_header('Content-Disposition', 'form-data', name="results", 
                              filename=self.wu.data["RESULT"])
            postdata.attach(result)
        if self.debug >= 2:
            print("Headers of postdata as a dictionary:")
            print(dict(postdata.items()))
        # Ugly hack: overwrite method for writing headers to suppress them
        # We pass the MIME headers to the request below via the headers= argument, 
        # and don't want them to occur again as part of the POST data
        if False:
            postdata2 = postdata.as_string(unixfrom=False) + "\n"
            if self.debug >= 2:
                print("Postdata as a string:")
                print(postdata2)
            postdata3 = bytes(postdata2, encoding="utf-8")
        else:
            fp = BytesIO()
            g = FixedBytesGenerator(fp)
            g.flatten(postdata, unixfrom=False)
            postdata3 = fp.getvalue() + b"\n"
        if self.debug >= 2:
            print("Postdata as a bytes array:")
            print(postdata3)
        url = SETTINGS["SERVER"] + "/" + SETTINGS["POSTRESULTPATH"]
        request = urllib.request.Request(url, data=postdata3, headers=dict(postdata.items()))
        conn = urllib.request.urlopen(request)
        print ("Server response:")
        for line in conn:
            print(line)
        conn.close()
        return True

    def result_exists(self):
        if "RESULT" in self.wu.data:
            filepath = SETTINGS["WORKDIR"] + "/" + self.wu.data["RESULT"]
            if not os.path.isfile(filepath):
                return False
            print ("Result file " + filepath + " already exists")
        print ("All result files already exist")
        return True

    def cleanup(self):
        print ("Cleaning up for workunit " + self.wu.get_id())
        if "RESULT" in self.wu.data:
            filepath = SETTINGS["WORKDIR"] + "/" + self.wu.data["RESULT"]
            print ("Removing result file " + filepath)
            os.remove(filepath)

    def process(self):
        # If all output files exist, send them, return WU as done
        # Otherwise, run commands in WU. If no error and all output 
        #   files exist, send them, return WU as done
        # print(str(wu))
        if not self.get_files():
            return False
        if not self.result_exists():
            self.run_commands()
        if not self.upload_result():
            return False
        self.cleanup()
        return True

def do_work():
    wu_filename = SETTINGS["DLDIR"] + "/" + SETTINGS["WU_FILENAME"]
    if not get_missing_file(SETTINGS["GETWUPATH"], wu_filename, options="clientid=" + SETTINGS["CLIENTID"]):
        return False
    wu = Workunit_Processor(wu_filename, int(SETTINGS["DEBUG"]))
    if not wu.process():
        return False
    print ("Removing workunit file " + wu_filename)
    os.remove(wu_filename)
    return True

if __name__ == '__main__':
    # Create command line parser from the keys in SETTINGS
    parser = argparse.ArgumentParser()
    for arg in REQUIRED_SETTINGS.keys():
        parser.add_argument('--' + arg.lower(), required = True,
        help=REQUIRED_SETTINGS[arg][1])
    for arg in OPTIONAL_SETTINGS.keys():
        parser.add_argument('--' + arg.lower(), required = False, 
            default=OPTIONAL_SETTINGS[arg][0], 
            help=OPTIONAL_SETTINGS[arg][1] + " (default: " + OPTIONAL_SETTINGS[arg][0] + ")")
    # Parse command line, store as dictionary
    args = vars(parser.parse_args())
    # Copy values to SETTINGS
    for arg in SETTINGS.keys():
        if arg.lower() in args:
            SETTINGS[arg] = args[arg.lower()].rstrip("/")

    # print (str(SETTINGS))

    do_work()
