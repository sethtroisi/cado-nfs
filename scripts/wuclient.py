#!/usr/bin/env python3

import os
import stat
import argparse
import shutil
import time
import urllib.request
import subprocess
from string import Template

# Settings which we require on the command line (no defaults)
REQUIRED_SETTINGS = {"CLIENTID" : "", "DLDIR" : "", "SERVER" : "", "WORKDIR" : ""}
# Optional settings (with defaults, overrideable on command line)
OPTIONAL_SETTINGS = {"WU_FILENAME" : "WU", "GETWUPATH" : "/cgi-bin/getwu", 
                     "POSTRESULTPATH" : "cgi-bin/upload.py"}
# Merge the two
SETTINGS = dict(list(REQUIRED_SETTINGS.items()) + list(OPTIONAL_SETTINGS.items()))

# Keys that must occur only once
WU_SCALAR_KEYS = ("WORKUNIT")
# Keys that can occur multiple times
WU_LIST_KEYS = ("FILE", "EXECFILE", "CHECKSUM", "COMMAND", "RESULT")

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
    if checksum:
        # FIXME CHECKSUM
        pass

def parse_WU(filename):
    wu = {}
    for key in WU_LIST_KEYS:
      wu[key] = []
    wu_file = open(filename)
    for line in wu_file:
        (key, value) = line.split(" ", 1)
        value = value.rstrip()
        if key in WU_SCALAR_KEYS:
            if key in wu:
                print("Error: key " + key + " redefined")
                return None
            wu[key] = value
        elif key in WU_LIST_KEYS:
            # (EXEC)FILE/CHECKSUM keys get special treatment, as CHECKSUM 
            # always refers to the previous (EXEC)FILE
            if key == "FILE" or key == "EXECFILE":
                wu["CHECKSUM"].append(None)
            if key == "CHECKSUM":
                if len(wu[key]) == 0 or wu[key][-1] != None:
                    print("Error: extraneous " + key)
                    return None
                wu[key][-1] = value
            else:
                wu[key].append(value)
        else:
            print("Error: key " + key + " not recognized")
            return None
    return wu

def upload_results(wu):
    pass

def process_WU(filename):
    print ("Processing work unit file " + filename + ":")
    for line in open(filename):
        print(line.rstrip())
    # If all output files exist, send them, return WU as done
    # Otherwise, run commands in WU. If no error and all output 
    #   files exist, send them, return WU as done
    wu = parse_WU(filename)
    if wu == None:
        return False
    # print(str(wu))
    for filename in wu["FILE"] + wu["EXECFILE"]:
        get_missing_file (filename, SETTINGS["DLDIR"] + '/' + filename)
    for filename in wu["EXECFILE"]:
        path = SETTINGS["DLDIR"] + '/' + filename
        mode = os.stat(path).st_mode
        if mode & stat.S_IXUSR == 0:
            print ("Setting executable flag for " + path)
            os.chmod(path, mode | stat.S_IXUSR)
    for command in wu["COMMAND"]:
        print (command)
        command = Template(command).safe_substitute(SETTINGS)
        print (command)
        rc = subprocess.call(command, shell=True)
        if rc != 0:
            return False
    upload_result(wu)


def do_work():
    wu_filename = SETTINGS["DLDIR"] + "/" + SETTINGS["WU_FILENAME"]
    get_missing_file(SETTINGS["GETWUPATH"], wu_filename)
    process_WU(wu_filename)

if __name__ == '__main__':
    # Create command line parser from the keys in SETTINGS
    parser = argparse.ArgumentParser()
    for arg in REQUIRED_SETTINGS.keys():
        parser.add_argument('--' + arg.lower(), required = True)
    for arg in OPTIONAL_SETTINGS.keys():
        parser.add_argument('--' + arg.lower(), required = False, 
            default=OPTIONAL_SETTINGS[arg])
    # Parse command line, store as dictionary
    args = vars(parser.parse_args())
    # Copy values to SETTINGS
    for arg in SETTINGS.keys():
        if arg.lower() in args:
            SETTINGS[arg] = args[arg.lower()].rstrip("/")

    print (str(SETTINGS))

    do_work()
