#!/usr/bin/env python3
import os
import sqlite3
import logging
import cadologger
import wudb
import cadotask
import cadoparams

# FIXME: Is this a Task object? Probably not
# Should this be in cadotask or in cadofactor?
class CompleteFactorization(object):
    """ The complete factorization, aggregate of the individual tasks """
    def __init__ (self, db, *args, **kwargs):
        self.polysel = cadotask.PolyselTask(*args, db=db, **kwargs)
        self.fb = cadotask.FactorBaseTask(self.polysel, *args, db=db, **kwargs)
        self.freerel = cadotask.FreeRelTask(self.polysel, *args, db=db, **kwargs)
        self.sieving = cadotask.SievingTask(self.polysel, self.fb, *args, 
                                            db=db, **kwargs)
        self.dup1 = cadotask.Duplicates1Task(self.sieving, *args, db=db, 
                                             **kwargs)
        self.dup2 = cadotask.Duplicates2Task(self.dup1, *args, db=db, 
                                             **kwargs)
        self.sing = cadotask.PurgeTask(self.polysel, self.freerel, self.dup2, *args, db=db, **kwargs)
        self.merge = cadotask.MergeTask(self.sing, *args, db=db, **kwargs)
        self.linalg = cadotask.LinAlgTask(self.merge, *args, db=db, **kwargs)
        self.sqrt = cadotask.SqrtTask(self.polysel, self.freerel, self.sieving, 
                             self.merge, self.linalg, *args, db=db, **kwargs)
    
    def run(self, *args, **kwargs):
        self.sqrt.run(*args, **kwargs)

if __name__ == '__main__':
    logger = cadologger.Logger()
    logger.addHandler(cadologger.ScreenHandler(lvl = logging.DEBUG))
    # logger.addHandler(cadologger.FileHandler(filename = "log", lvl = logging.DEBUG))
    parameters = cadoparams.Parameters()
    parameters._readfile(open("parameters"))
    tasksparams = parameters.myparams(("workdir", "name"), "tasks")
    wudb_file = tasksparams["workdir"] + os.sep + tasksparams["name"] + ".db"
    logger.info ('Opening database file "%s"', wudb_file)
    # dbconn = sqlite3.connect("wudb")
    logger.info ("Beginning factorization")
    factorjob = CompleteFactorization(wudb_file, parameters = parameters)
    factorjob.run()
