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

# Parameters are stored in the DB. Do we want one table for each task?
# Lump them all together in one big table? To avoid name conflicts,
# we probably want a hierarchical parameter structure with different
# namespaces for different tasks, which suggests one table per task.
# Say, one table per task. Columns: parameter, value

# Should task dependecy be modeled with inheritance? Is the graph traversal
# algorithm of super() the one we want? Would need our own namespace mechanism
# to keep variables of different tasks separate - messy

class Task(object):
    """ An class that represents one task that needs to be processed.
    """
    name = None
    title = None
    dependencies = ()
    programs = ()

    def __init__(self, dbconn, defaults = None, override = None):
        # Derived class must define name
        self.dbconn = dbconn
        self.dbparam = wudb.DictDbAccess(dbconn, self.name)
        if defaults:
            for key in self.myparams(defaults):
                    self.dbparam.setdefault(key, defaults[key])
        if override:
            self.dbparam.update(self.myparams(override))

    def timestamp(self):
        """ Returns timestamp of completion, or None if not completed """
        if "time_finished" in self.dbparam:
            return self.dbparam["time_finished"]
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
            self.dbparam["time_finished"] = gettime()
        else:
            self.dbparam["time_finished"] = None

    def run(self):
        """ Runs the prerequisites. Sub-classes should call this.
        """
        logger = cadologger.Logger()
        # Check/run the prerequisites
        if not self.dependencies is None:
            for d in self.dependencies:
                task = d(self.dbconn)
                task.run()
                # Check if prereq is newer than our timestamp
                if not self.prereq_ok(task):
                    logger.info("Prerequisite %s not ok, setting %s as not "
                                "done", task.title, self.title)
                    self.mark_done(False)

    def submit(self, command):
        ''' Submit a command, or list of commands, that need to be run '''
        pass
    
    def status(self):
        ''' Check status of previously submitted commands '''
        pass

class DistributedTask(Task):
    def __init__(self, *args, server = None, **kwargs):
        self.server = server
        super().__init__(*args, **kwargs)

# How do we handle transitivity? Should a task depend on all tasks that need
# to run before, or just on the "previous" one and let transitivity handle
# earlier ones?

#class Parameters(object):
#    """ Chooses parameters for a factorization """
#    name = "parameters"
#    title = "Parameter selection"
#    dependencies = None
#    class Table(wudb.DbTable):
#        name = "parameters"
#        fields = (
#            ("rowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"), 
#            ("N", "INTEGER", ""), 
#        )
#        primarykey = fields[0][0]
#        references = None
#        index = {}
#        pass
#    pass

class PolyselTask(DistributedTask):
    """ Finds a polynomial, uses client/server """
    name = "polyselect"
    title = "Polynomial Selection"
    dependencies = ()
    programs = (cadoprograms.Polyselect2l,)
    pass

class FactorBaseTask(Task):
    """ Generates the factor base for the polynomial(s) """
    name = "factorbase"
    title = "Generate Factor Base"
    dependencies = (PolyselTask,)
    pass

class FreeRelTask(Task):
    """ Generates free relations for the polynomial(s) """
    name = "freerel"
    title = "Generate Free Relations"
    dependencies = (PolyselTask,)
    pass

class SievingTask(DistributedTask):
    """ Does the sieving, uses client/server """
    name = "sieving"
    title = "Sieving"
    dependencies = (PolyselTask, FactorBaseTask)
    programs = (cadoprograms.Las,)
    pass

class DuplicatesTask(Task):
    """ Removes duplicate relations """
    name = "duplicates"
    title = "Filtering: Duplicate Removal"
    dependencies = (SievingTask,)
    pass

class SingletonsTask(Task):
    """ Removes singletons and computes excess """
    name = "singletons"
    title = "Filtering: Singleton removal"
    dependencies = (DuplicatesTask,)
    pass

class MergeTask(Task):
    """ Merges relations """
    name = "merging"
    title = "Filtering: Merging"
    dependencies = (SingletonsTask,)
    pass

class LinAlgTask(Task):
    """ Runs the linear algebra step """
    name = "linalg"
    title = "Linear Algebra"
    dependencies = (MergeTask,)
    pass

class SqrtTask(Task):
    """ Runs the square root """
    name = "sqrt"
    title = "Square Root"
    dependencies = (MergeTask, LinAlgTask)
    pass

class CompleteFactorizationTask(Task):
    """ Finishes the factorization, tests the factors, cleans up, etc. """
    name = "complete"
    title = "Complete Factorization"
    dependencies = (SqrtTask,)
    pass


