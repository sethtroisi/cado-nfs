import sqlite3
import wudb
import cadoprograms
import cadologger

# TODO:
# Job status. Should probably be in the DB - no need to have another file sitting around with another file format to parse
# Check which tasks need to be done:
#   1. Parameter generation
#   2. Polynomial selection (add WUs to DB until enough polynomials checked/good enough found)
#   3. Factorbase (if hosted on server)
#   4. Free relations
#   5. Sieving until enough relations received (generating WUs on the fly, adding to DB, distributing to clients, checking results)
#   6. Remove duplicates
#   7. Remove singletons, check, go back to Sieving if necessary (simply increase rels_wanted)
#   8. Merge
#   9. Call BWC
#   10. SQRT
#   11. Check complete factorization

# Some parameters are provided by the param file but can change during
# the factorization, like rels_wanted. On one hand, we want automatic 
# updates to parameters to be stored in the DB, otoh, we want to allow
# externally setting new parameters. Need to distinguish between new 
# external parameters that overwrite DB, and old external parameters 
# that don't overwrite. Or maybe two ways to specify external params:
# --defaults which does not overwrite, and --forceparam which does

class Task(object):
    """ An class that represents one task that needs to be processed. """
    
    def __init__(self, dependencies, db = None, defaults = None):
        self.dbconn = sqlite3.connect(db)
        self.dependencies = dependencies
        # Derived class must define name
        self.param = wudb.DictDbAccess()
        self.param.attachdb(db, self.name)
        self.param_path = "tasks." + self.name
        if defaults:
            mydefaults = defaults._myparams(self.params, self.param_path)
            self.param.setdefault(None, mydefaults)
        self.programparam = []
        for prog in self.programs:
            name = self.name + '_' + prog.binary
            progparam = wudb.DictDbAccess()
            progparam.attachdb(db, name)
            self.programparam.append(progparam)
            path = self.param_path + '.' + prog.binary
            if defaults:
                mydefaults = defaults._myparams(prog.params, path)
                progparam.setdefault(None, mydefaults)
    
    def setparam(override):
        self.param.update(override._myparams(self.params, self.param_path))
    
    def timestamp(self):
        """ Returns timestamp of completion, or None if not completed """
        if "time_finished" in self.param:
            return self.param["time_finished"]
        else:
            return None
    
    def prereq_ok(self, other):
        """ Returns whether self finished more recently than the other task 
        did. If self or other did not finish yet, returns False
        """
        if self.timestamp() is None or other.timestamp() is None:
            return False
        else:
            return self.timestamp() >= other.timestamp()

    def mark_done(self, done):
        if not done and self.timestamp() is None:
            # Nothing to do
            return
        if done:
            self.param["time_finished"] = gettime()
        else:
            self.param["time_finished"] = None

    def run(self):
        """ Runs the prerequisites. Sub-classes should call this. """
        logger = cadologger.Logger()
        # Check/run the prerequisites
        if not self.dependencies is None:
            for task in self.dependencies:
                task.run()
                # Check if prereq is newer than our timestamp
                if not self.prereq_ok(task):
                    logger.info("Prerequisite %s not ok, setting %s as not "
                                "done", task.title, self.title)
                    self.mark_done(False)

    def submit(self, commands, inputfiles, outputfiles, tempfiles):
        ''' Submit a command that needs to be run. Returns a handle
        which can be used for status check.

        The inputfiles parameters is a list of input files that program 
        needs; they are not automatically filled into the command line(s),
        but need to be listed on the command line(s) explicitly. The are 
        input files list is used to generate FILE lines in work units, 
        for example.
        
        The outputfiles parameter lists the output files that running the 
        command will produce. They are used to produce RESULT lines in 
        workunits. Note that the final file names of the output files may 
        differ, for example, when the output files are uploaded to the server
        and stored under a unique file name.

        The tempfiles parameter lists temporary files that should be deleted 
        after the commands ran successfully. There are used to generate DELETE 
        lines in work units.
        '''
        
        # We should actually generate the workunits here and feed them
        # to the workunit processor. Avoids duplicate code and tests 
        # more of the code path. Status of the workunits should be kept
        # in a wudb, just like the server would. I.e., we do everything
        # here as in the client/server setup, except we skip the server
        # and handing WUs to the wu processor directly
        pass
    
    def status(self, handle):
        ''' Check status of previously submitted commands 
        
        Returns the execution status and the list of output files produced
        '''
        return (status, outputfiles)

class ClientServerTask(Task):
    def __init__(self, server = None, *args, **kwargs):
        self.server = server
        super().__init__(*args, **kwargs)

    def submit(self, commands, inputfiles, outputfiles, tempfiles):
        ''' Submit a command that needs to be run. Uses client/server.
        Returns an index which can be used for status check
        '''
        pass
    
    def status(self, index):
        ''' Check status of previously submitted commands '''
        pass


# Each task has positional parameters with references to the tasks from 
# which it receives its inputs. This will be used to query the referenced 
# tasks' parameters, list of output files, etc. - all the info the current 
# task needs to run and which is produced by previous tasks. Since the 
# dependencies are an intrinsic property of a given task, they are passed 
# as positional parameters rather than as, say, a list. The abstract 
# superclass Task receives the dependecies as a list; this list is stored 
# and used for graph traversal which is implemented in Task.

class PolyselTask(Task):
    """ Finds a polynomial, uses client/server """
    name = "polyselect"
    title = "Polynomial Selection"
    programs = (cadoprograms.Polyselect2l,)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, dependencies = None, **kwargs)

    def run(self):
        # Make command line for polselect2l, run it. 
        # Whole range in one go for now
        args = ()
        kwargs = self.programparam[0]
        p = self.programs[0](args, kwargs)
        p.run()
        p.wait()

class FactorBaseTask(Task):
    """ Generates the factor base for the polynomial(s) """
    name = "factorbase"
    title = "Generate Factor Base"
    programs = ()
    def __init__(self, polsel, *args, **kwargs):
        super().__init__(*args, dependencies = (polsel,), **kwargs)

class FreeRelTask(Task):
    """ Generates free relations for the polynomial(s) """
    name = "freerel"
    title = "Generate Free Relations"
    programs = ()
    def __init__(self, polsel, *args, **kwargs):
        super().__init__(*args, dependencies = (polsel,), **kwargs)

class SievingTask(Task):
    """ Does the sieving, uses client/server """
    name = "sieving"
    title = "Sieving"
    programs = ()
    def __init__(self, polsel, factorbase, *args, **kwargs):
        super().__init__(*args, dependencies = (polsel, factorbase), **kwargs)

class DuplicatesTask(Task):
    """ Removes duplicate relations """
    name = "duplicates"
    title = "Filtering: Duplicate Removal"
    programs = ()
    def __init__(self, sieving, *args, **kwargs):
        super().__init__(*args, dependencies = (sieving,), **kwargs)

class SingletonsTask(Task):
    """ Removes singletons and computes excess """
    name = "singletons"
    title = "Filtering: Singleton removal"
    programs = ()
    def __init__(self, duplicates, *args, **kwargs):
        super().__init__(*args, dependencies = (duplicates,), **kwargs)

class MergeTask(Task):
    """ Merges relations """
    name = "merging"
    title = "Filtering: Merging"
    programs = ()
    def __init__(self, freerel, unique, *args, **kwargs):
        super().__init__(*args, dependencies = (freerel, unique), **kwargs)

class LinAlgTask(Task):
    """ Runs the linear algebra step """
    name = "linalg"
    title = "Linear Algebra"
    programs = ()
    def __init__(self, merge, *args, **kwargs):
        super().__init__(*args, dependencies = (merge,), **kwargs)

class SqrtTask(Task):
    """ Runs the square root """
    name = "sqrt"
    title = "Square Root"
    programs = ()
    def __init__(self, polsel, freerel, sieving, merge, linalg, 
                 *args, **kwargs):
        dep = (polsel, freerel, sieving, merge)
        super().__init__(*args, dependencies = dep, **kwargs)

# FIXME: Is this a Task object? Probably not
# Should this be in cadotask or in cadofactor?
class CompleteFactorization(object):
    """ The complete factorization, aggregate of the individual tasks """
    def __init__ (self, db, *args, **kwargs):
        self.polysel = PolyselTask(*args, db=db, **kwargs)
        self.fb = FactorBaseTask(self.polysel, *args, db=db, **kwargs)
        self.freerel = FreeRelTask(self.polysel, *args, db=db, **kwargs)
        self.sieving = SievingTask(self.polysel, self.fb, *args, db=db, 
                                   **kwargs)
        self.dup = DuplicatesTask(self.sieving, *args, db=db, **kwargs)
        self.sing = SingletonsTask(self.dup, *args, db=db, **kwargs)
        self.merge = MergeTask(self.freerel, self.sing, *args, db=db, **kwargs)
        self.linalg = LinAlgTask(self.merge, *args, db=db, **kwargs)
        self.sqrt = SqrtTask(self.polysel, self.freerel, self.sieving, 
                             self.merge, self.linalg, *args, db=db, **kwargs)
    def run(self):
        self.sqrt.run()
