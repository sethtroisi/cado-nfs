#!/usr/bin/env python3
import os
import argparse
import logging
import cadologger
import cadotask
import cadoparams

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Integer Factorisation with the Number Field Sieve')
    parser.add_argument("--screenlog", help="Screen logging level, e.g., INFO/COMMAND/DEBUG", default="INFO", metavar="LEVEL")
    parser.add_argument("parameters", help="A file with the parameters to use")
    args = parser.parse_args()
    paramfile = args.parameters
    screenlvlname = args.screenlog
    
    parameters = cadoparams.Parameters()
    parameters._readfile(open(paramfile))
    tasksparams = parameters.myparams(("workdir", "name"), "tasks")
    
    screenlvl = getattr(cadologger, screenlvlname.upper())
    logger = logging.getLogger()
    logger.addHandler(cadologger.ScreenHandler(lvl = screenlvl))
    
    logger.debug("Root parameter dictionary:\n%s" % parameters)
    
    # Make working directory, if it does not exist
    directory = tasksparams["workdir"]
    if not os.path.isdir(directory):
        os.makedirs(directory)
    
    # Add a logger to capture the command lines of programs we run
    cmdfilename = tasksparams["workdir"] + os.sep + tasksparams["name"] + ".cmd"
    logger.addHandler(cadologger.CmdFileHandler(cmdfilename))
    
    # Add a logger to write debugging information to a log file
    logfilename = tasksparams["workdir"] + os.sep + tasksparams["name"] + ".log"
    logger.addHandler(cadologger.FileHandler(filename = logfilename, lvl = logging.DEBUG))
    
    wudb_file = tasksparams["workdir"] + os.sep + tasksparams["name"] + ".db"
    factorjob = cadotask.CompleteFactorization(db=wudb_file, parameters = parameters, path_prefix = [])
    factorjob.run()
