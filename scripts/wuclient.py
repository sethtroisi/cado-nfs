#!/usr/bin/env python3

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
from string import Template

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

def get_file(urlpath, dlpath = None):
    # print('get_file("' + urlpath + '", "' + dlpath + '")')
    if dlpath == None:
        dlpath = SETTINGS["DLDIR"] + "/" + os.basename(urlpath) # FIXME: should be url base name
    url = SETTINGS["SERVER"] + "/" + urlpath
    print ("Downloading " + url + " to " + dlpath);
    request = urllib.request.urlopen(url)
    file = open(dlpath, "wb")
    shutil.copyfileobj (request, file)
    file.close()
    request.close()

def get_missing_file(urlpath, filename, checksum = None):
    # print('get_missing_file("' + urlpath + '", "' + filename + '", ' + str(checksum) + ')')
    if not os.path.isfile(filename):
        get_file(urlpath, filename)
        # FIXME CHECKSUM
    else:
        print (filename + " already exists, not downloading")
    if checksum:
        # FIXME CHECKSUM
        pass
    return True

def _do_nothing(msg):
    pass

class Workunit():
    # Keys that must occur only once
    SCALAR_KEYS = ("WORKUNIT", "RESULT")
    # Keys that can occur multiple times
    LIST_KEYS = ("FILE", "EXECFILE", "CHECKSUM", "COMMAND")

    parsed = None
    exitcode = 0 # will get set if any command exits with code != 0
    debug = 1 # Controls debugging output
    
    def parse(self, filepath):
        wu = {}
        for key in self.LIST_KEYS:
          wu[key] = [] # Init to empty list
        if self.debug >= 1:
            print ("Parsing workunit from file " + filepath)
        wu_file = open(filepath)
        for line in wu_file:
            (key, value) = line.split(" ", 1)
            value = value.rstrip() # Drop trailing whitespace, incl. CR/LF
            if key in self.SCALAR_KEYS:
                if key in wu:
                    print("Error: key " + key + " redefined")
                    return False
                wu[key] = value
            elif key in self.LIST_KEYS:
                # (EXEC)FILE/CHECKSUM keys get special treatment, as CHECKSUM 
                # always refers to the previous (EXEC)FILE
                if key == "FILE" or key == "EXECFILE":
                    wu["CHECKSUM"].append(None)
                if key == "CHECKSUM":
                    if len(wu[key]) == 0 or wu[key][-1] != None:
                        print("Error: extraneous " + key)
                        return False
                    wu[key][-1] = value
                else:
                    wu[key].append(value)
            else:
                print("Error: key " + key + " not recognized")
                return False
        wu_file.close()
        self.parsed = wu
        if self.debug >= 1:
            print (" done, workunit ID is " + self.parsed["WORKUNIT"])
        return True

    def get_files(self):
        for filename in self.parsed["FILE"] + self.parsed["EXECFILE"]:
            if not get_missing_file (filename, SETTINGS["DLDIR"] + '/' + filename):
                return False
        for filename in self.parsed["EXECFILE"]:
            path = SETTINGS["DLDIR"] + '/' + filename
            mode = os.stat(path).st_mode
            if mode & stat.S_IXUSR == 0:
                print ("Setting executable flag for " + path)
                os.chmod(path, mode | stat.S_IXUSR)
        return True
    
    def run_commands(self):
        for command in self.parsed["COMMAND"]:
            command = Template(command).safe_substitute(SETTINGS)
            if self.debug >= 1:
                print ("Running command for " + self.parsed["WORKUNIT"] + ": " + command)
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
        WUid = MIMEText(self.parsed["WORKUNIT"])
        WUid.add_header('Content-Disposition', 'form-data', name="WUid")
        postdata.attach(WUid)
        if self.exitcode > 0:
            rc = MIMEText(str(self.exitcode))
            WUid.add_header('Content-Disposition', 'form-data', name="exitcode")
            postdata.attach(rc)
        if "RESULT" in self.parsed:
            filepath = SETTINGS["WORKDIR"] + "/" + self.parsed["RESULT"]
            print ("Adding result file " + filepath + " to upload")
            file = open(filepath, "rb")
            filedata = file.read()
            file.close()
            result = MIMEApplication(filedata)
            result.add_header('Content-Disposition', 'form-data', name="results", 
                              filename=self.parsed["RESULT"]);
            postdata.attach(result)
        if self.debug >= 2:
            print("Headers of postdata as a dictionary:")
            print(dict(postdata.items()))
        # Ugly hack: overwrite method for writing headers to suppress them
        postdata._write_headers = _do_nothing
        postdata2 = postdata.as_string(unixfrom=False) + "\n"
        if self.debug >= 2:
            print("Postdata as a string:")
            print(postdata2)
        postdata3 = bytes(postdata2, encoding="utf-8")
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

    def cleanup(self):
        print ("Cleaning up for workunit " + self.parsed["WORKUNIT"])
        if "RESULT" in self.parsed:
            filepath = SETTINGS["WORKDIR"] + "/" + self.parsed["RESULT"]
            print ("Removing result file " + filepath)
            os.remove(filepath)

    def process_file(self, filepath):
        if int(SETTINGS["DEBUG"]) > 0:
            print ("Processing work unit file " + filepath + ":")
            for line in open(filepath):
                print(line.rstrip())
        # If all output files exist, send them, return WU as done
        # Otherwise, run commands in WU. If no error and all output 
        #   files exist, send them, return WU as done
        if not self.parse(filepath):
            return False
        # print(str(wu))
        if not self.get_files():
            return False
        if not self.run_commands():
            return False
        if not self.upload_result():
            return False
        self.cleanup()

def do_work():
    wu_filename = SETTINGS["DLDIR"] + "/" + SETTINGS["WU_FILENAME"]
    if not get_missing_file(SETTINGS["GETWUPATH"], wu_filename):
        return False
    wu = Workunit()
    if not wu.process_file(wu_filename):
        return False
    print ("Removing workunit file " + wu_filename)
    os.remove(wu_filename)

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
