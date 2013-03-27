#!/usr/bin/env python3
import sqlite3
import logging
import cadologger
import wudb
import cadotask

# TODO:
# read parameter file and number to be factored
# run the tasks with those parameters

if __name__ == '__main__':
    wudb_file = "wudb"
    logger = cadologger.Logger()
    logger.addHandler(cadologger.ScreenHandler(lvl = logging.DEBUG))
    # logger.addHandler(cadologger.FileHandler(filename = "log", lvl = logging.DEBUG))
    logger.info ('Opening database file "%s"', wudb_file)
    # dbconn = sqlite3.connect("wudb")
    logger.info ("Beginning factorization")
    factorjob = cadotask.CompleteFactorization(wudb_file)
    factorjob.run()
