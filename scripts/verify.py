#!/usr/bin/env python3
import sys
import os
import re
import argparse
import sqlite3
import subprocess
import wudb

def do_check(acc, wu, programs):
    result = None
    if wu["files"] is None:
        return result
    for file in wu["files"]:
        (filename, filepath) = (file["filename"], file["path"])
        for (progname, pattern) in programs:
            m = re.match(pattern, filename)
            if m is None:
                continue
            assert m.start() == 0
            # Check that the whole string matched the pattern,
            # not just a prefix
            if not m.end() == len(filename):
                continue
            command = progname + " " + filepath
            print ("Running: " + command)
            child = subprocess.Popen(command, shell=True, 
                                     stdout = subprocess.PIPE, 
                                     stderr = subprocess.PIPE)
            (child_stdout, child_stderr) = child.communicate()
            if child.returncode == 0:
                # acc.update({'status':wudb.WuStatus.VERIFIED_OK})
                if result is None:
                    result = True
            else:
                # acc.update({'status':wudb.WuStatus.VERIFIED_ERROR})
                result = False
    return result

if __name__ == "__main__":
    wudb_name = "wudb"
    parser = argparse.ArgumentParser()
    parser.add_argument('-file', action = 'store')
    parser.add_argument('-program', nargs = 2, required = True, 
                        action = 'store')
    parser.add_argument('-wudb', action = 'store', default = wudb_name)
    parser.add_argument('-received', action = 'store_true')
    parser.add_argument('-update', action = 'store_true')
    args = parser.parse_args()

    if "wudb" in args:
        wudb_name = args.wudb

    conn = sqlite3.connect(wudb_name)
    acc = wudb.WuAccess(conn)

    if "file" in args:
        # find wu for this file
        # do_check(acc, wu)
        pass

    if args.received:
        # Test all workunits whose status is RECEIVED_OK
        programs = [(args.program[0], args.program[1])]
        received = acc.query(eq={"status": wudb.WuStatus.RECEIVED_OK})
        print ("Found " + str(len(received)) + " workunits")
        for wu in received:
            result = do_check(acc, wu, programs)
            if result is None:
                print ("Warning: none of the supplied programs matched for workunit " 
                       + wu["wuid"])
            else:
                if result == True:
                    print ("Workunit " + wu["wuid"] + " verified ok")
                elif result == False:
                    print ("Workunit " + wu["wuid"] + " verified with error")
            if args.update:
                acc.verification(wu["wuid"], result)
