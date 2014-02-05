#!/usr/bin/env python3
import os
import sys
import argparse
import logging
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
    parser.add_argument("--old", help="Use old parameter file format",
                        action="store_true")
    parser.add_argument("parameters", help="A file with the parameters to use")
    parser.add_argument("options", metavar="OPTION", help="An option as in "
                        "parameter file (format: key=value)", nargs="*")
    args = parser.parse_args()
    paramfile = args.parameters
    screenlvlname = args.screenlog
    filelvlname = args.filelog
    
    parameters = cadoparams.Parameters()
    if args.old:
        parameters.read_old_defaults()
    parameters.readfile(paramfile, old_format = args.old)
    parameters.readparams(args.options, old_format = args.old)
    tasksparams = parameters.myparams({"workdir": str, "name": str}, "tasks")
    
    screenlvl = getattr(cadologger, screenlvlname.upper())
    logger = logging.getLogger()
    logger.addHandler(cadologger.ScreenHandler(lvl = screenlvl))
    
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
        print("The logarithms of the factor base elements are in %s" %
                tasksparams["workdir"] + os.sep + tasksparams["name"] +
                ".reconstructlog.dlog")
