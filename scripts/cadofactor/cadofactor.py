#!/usr/bin/env python3
import os
import sys
import argparse
import logging

import re
cado_python_libs_path="@CMAKE_INSTALL_PREFIX@/@LIBSUFFIX@"
if re.search("/", cado_python_libs_path):
    sys.path.append(cado_python_libs_path)

import cadologger
import cadotask
import cadoparams
import itertools
from cadocommand import shellquote

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Integer Factorisation with "
                                     "the Number Field Sieve")
    parser.add_argument("--screenlog", help="Screen logging level, e.g., "
                        "INFO/COMMAND/DEBUG", default="INFO", metavar="LEVEL")
    parser.add_argument("--filelog", help="Log file logging level, e.g., "
                        "INFO/COMMAND/DEBUG", default="DEBUG", metavar="LEVEL")
    parser.add_argument("parameters", help="A file with the parameters to use")
    parser.add_argument("options", metavar="OPTION", help="An option as in "
                        "parameter file (format: key=value)", nargs="*")
    args = parser.parse_args()
    paramfile = args.parameters
    screenlvlname = args.screenlog
    filelvlname = args.filelog
    
    screenlvl = getattr(cadologger, screenlvlname.upper())
    logger = logging.getLogger()
    logger.addHandler(cadologger.ScreenHandler(lvl = screenlvl))

    parameters = cadoparams.Parameters()
    parameters.readfile(paramfile)
    parameters.readparams(args.options)
    tasksparams = parameters.myparams({"workdir": str, "name": str}, "tasks")
    
    # Make working directory, if it does not exist
    directory = tasksparams["workdir"]
    if not os.path.isdir(directory):
        os.makedirs(directory)
    
    # Add a logger to capture the command lines of programs we run
    cmdfilename = tasksparams["workdir"] + os.sep + tasksparams["name"] + ".cmd"
    logger.addHandler(cadologger.CmdFileHandler(cmdfilename))
    
    # Add a logger to write debugging information to a log file
    filelvl = getattr(cadologger, filelvlname.upper())
    logfilename = tasksparams["workdir"] + os.sep + tasksparams["name"] + ".log"
    filehandler = cadologger.FileHandler(filename = logfilename, lvl = filelvl)
    logger.addHandler(filehandler)
    
    logger.info("Command line parameters: %s", 
                " ".join(map(shellquote, sys.argv)))

    logger.debug("Root parameter dictionary:\n%s", parameters)

    # Write a snapshot of the parameters to a file
    for counter in itertools.count():
        snapshot_filename = "%s%s%s.parameters_snapshot.%d" % \
                (tasksparams["workdir"], os.sep, tasksparams["name"], counter)
        if not os.path.isfile(snapshot_filename):
            break
    with open(snapshot_filename, "w") as snapshot_file:
        logger.debug("Writing parameter snapshot to %s", snapshot_filename)
        snapshot_file.write(str(parameters))
        snapshot_file.write("\n")
    
    wudb_file = tasksparams["workdir"] + os.sep + tasksparams["name"] + ".db"
    factorjob = cadotask.CompleteFactorization(db=wudb_file,
                                               parameters = parameters,
                                               path_prefix = [])
    factors = factorjob.run()
    
    dlp_param = parameters.myparams(("dlp",), "")
    dlp = dlp_param.get("dlp", False)
    if not dlp:
        if factors is None:
            sys.exit("Error occurred, terminating")
        else:
            print(" ".join(factors))
    else:
        p = int(factors[0])
        ell = int(factors[1])
        log2 = int(factors[2])
        log3 = int(factors[3])
        assert (p-1) % ell == 0
        assert pow(3, log2*((p-1) // ell), p) == pow(2, log3*((p-1) // ell), p)
        print("p = " + str(p))
        print("ell = " + str(ell))
        print("log2 = " + str(log2))
        print("log3 = " + str(log3))
        print("The other logarithms of the factor base elements are in %s" %
                tasksparams["workdir"] + os.sep + tasksparams["name"] +
                ".reconstructlog.dlog")
