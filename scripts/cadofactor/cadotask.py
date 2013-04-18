import sqlite3
from datetime import datetime
import re
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
    """ A base class that represents one task that needs to be processed. """
    
    def __init__(self, dependencies, db = None, parameters = None):
        ''' Sets up a database connection and a DB-backed dictionary for 
        parameters. Reads parameters from DB, and merges with hierarchical
        parameters in the parameters argument. Parameters passed in by 
        parameters argument do not override values in the DB-backed 
        parameter dictionary.
        '''
        self.dbconn = sqlite3.connect(db)
        self.dependencies = dependencies
        # Derived class must define name
        self.param = wudb.DictDbAccess(db, self.name)
        self.param_path = "tasks." + self.name
        self.logger = cadologger.Logger()
        self.logger.debug("Enter Task.__init__(%s)", self.name)
        # self.logger.debug("parameters = %s", parameters)
        self.parampath = "tasks" + parameters.get_sep() + self.name
        # DB-backed dictionary with parameters for this task
        self.taskparams = wudb.DictDbAccess()
        # Derived class must define name
        self.taskparams.attachdb(db, self.tablename())
        # Set default parametes for this task, if any are given
        if parameters:
            defaults = parameters.myparams(self.paramnames, self.parampath)
            self.taskparams.setdefault(None, defaults)
        # Set default parameters for our programs
        self.progparams = []
        for prog in self.programs:
            progparams = wudb.DictDbAccess()
            self.progparams.append(progparams)
            progparams.attachdb(db, self.prog_tablename(prog))
            if parameters:
                defaults = parameters.myparams(
                    prog.params_dict(), [self.parampath, prog.name])
                progparams.setdefault(None, defaults)
        self.logger.debug("Exit Task.__init__(%s)", self.name)
        return
    
    def tablename(self):
        """ Return the table name for the DB-backed dictioanry with parameters
        for the current task """
        # Maybe replace SQL-disallowed characters here, like '.' ?
        return self.name
    
    def prog_tablename(self, prog):
        """ Return the table name for the DB-backed dictionary with parameters 
        for program prog
        """
        return self.tablename() + '_' + prog.name
    
    def timestamp(self):
        """ Returns timestamp of completion, or None if not completed """
        if "time_finished" in self.taskparams:
            return self.taskparams["time_finished"]
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
    
    def is_done(self):
        return not (self.timestamp() is None)
    
    def mark_done(self, done):
        if not done and self.timestamp() is None:
            # Nothing to do
            return
        if done:
            self.taskparams["time_finished"] = str(datetime.now())
        else:
            self.taskparams["time_finished"] = None

    def run(self, parameters = None):
        ''' Runs the prerequisites. Sub-classes should call this first in 
        their run() method.
        
        Parameters passed in by parameters DO override values in the
        DB-backed parameter dictionary. 
        '''
        self.logger.debug("Enter Task.run(%s)", self.name)
        self.logger.debug("Task.run(%s): self.is_done() = %s", 
                          self.name, self.is_done())
        # Check/run the prerequisites
        if not self.dependencies is None:
            for task in self.dependencies:
                self.logger.debug("Task.run(%s): Running prerequisite %s",
                                  self.name, task.name)
                task.run(parameters)
                # Check if prereq is newer than our timestamp
                if self.is_done() and not self.prereq_ok(task):
                    self.logger.info(
                        "Prerequisite %s not ok, setting %s as not done", 
                        task.title, self.title)
                    self.mark_done(False)
        
        # Set parameters for our programs
        self.progparams = []
        for prog in self.programs:
            self.progparams.append(wudb.DictDbAccess(db, self.name + '_' + prog.name))
            if parameters:
                parampath = self.param_path + '.' + prog.name
                mydefaults = defaults._myparams(prog.params, parampath)
                self.progparams[-1].update(None, mydefaults)
        if parameters:
            # Override this task's parameters, if parameters are given
            update = parameters.myparams(self.paramnames, self.parampath)
            self.param.update(None, update)
            # Override the programs' parameters, if parameters are given
            for (index, prog) in enumerate(self.programs):
                update = parameters.myparams(
                    prog.params_dict(), [self.parampath, prog.name])
                self.progparams[index].update(None, update)
        
        self.logger.debug("Exit Task.run(" + self.name + ")")
        return
    
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

class Polynomial(object):
    # Keys that can occur in a polynomial file in their preferred ordering,
    # and whether the key is mandatory or not. The preferred ordering is used
    # for turning a polynomial back into a string. 
    keys = (
        ("n", True),
        ("Y0", True),
        ("Y1", True),
        ("c0", True),
        ("c1", True),
        ("c2", True),
        ("c3", True),
        ("c4", True),
        ("c5", False),
        ("c6", False),
        ("m", True),
        ("skew", True)
        )
    def __init__(self, lines):
        """ Parse a polynomial file in the syntax as produced by polyselect2l 
        """
        self.poly = {}
        self.E = 0.
        for line in lines:
            # print ("Parsing line: >%s<" % line)
            # If there is a "No polynomial found" message anywhere in the
            # file, we assume that there is no polynomial. This assumption
            # will be false if, e.g., files get concatenated
            if re.match("No polynomial found", line):
                self.poly = None
                return 
            # If this is a line telling the Murphy E value, 
            # extract the value and store it
            match = re.match("\s*#\s*MurphyE\s*\(.*\)=(.*)$", line)
            if match:
                self.E = float(match.group(1))
                continue
            # Drop comment, strip whitespace
            l = line.split('#', 1)[0].strip()
            # print ("Without comment: >%s<" % l)
            # If nothing is left, process next line
            if not l:
                continue
            # All remaining lines must be of the form "x: y"
            a = l.split(":")
            if not len(a) == 2:
                raise Exception("Invalid line %s in file %s" % (l, filename))
            (key, value) = a
            if not key in dict(self.keys):
                raise Exception("Invalid key %s in line %s in file %s" %
                                (key, l, filename))
            self.poly[key] = value
        for (key, isrequired) in self.keys:
            if isrequired and not key in self.poly:
                raise Exception("Key %s missing in polynomial in file %s" % 
                                (key, filename))
        return
    
    def __bool__(self):
        # Python 3 calls it __bool__()
        return not self.poly is None
    
    def __nonzero__(self):
        # and Python 2 calls it __nonzero__()
        return self.__bool__()
    
    def __str__(self):
        if self.poly is None:
            return ""
        arr = [key + ": " + self.poly[key] for (key, req) in self.keys 
               if key in self.poly]
        return "\n".join(arr)
    
    def __cmp__(self, other):
        """ Compares the Murphy E values """
        return int(other.E < self.E) - int(self.E < other.E)
    # Oh joy, Python 3 did away with __cmp__()
    def __lt__(self, other):
        return self.__cmp__(other) < 0
    def __le__(self, other):
        return self.__cmp__(other) <= 0
    def __eq__(self, other):
        return self.__cmp__(other) == 0
    def __ne__(self, other):
        return self.__cmp__(other) != 0
    def __ge__(self, other):
        return self.__cmp__(other) >= 0
    def __gt__(self, other):
        return self.__cmp__(other) > 0

    def setE(self, E):
        self.E = float(E)


class PolyselTask(Task):
    """ Finds a polynomial, uses client/server """
    name = "polyselect"
    title = "Polynomial Selection"
    programs = (cadoprograms.Polyselect2l,)
    paramnames = ("adrange", "P", "N", "admin", "admax", "name")
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, dependencies = None, **kwargs)
    
    def is_done(self):
        self.logger.debug ("PolyselTask.is_done(): Task parameters: %s", 
                           self.taskparams)
        return int(self.taskparams["admin"]) >= int(self.taskparams["admax"])
    
    def run(self, parameters = None):
        # Make command line for polselect2l, run it. 
        # Whole range in one go for now
        self.logger.debug("Enter PolyselTask.run(" + self.name + ")")
        super().run(parameters = parameters)
        
        print ("PolyselTask.run(): Task parameters:")
        print(self.taskparams)
        
        if "bestpoly" in self.taskparams:
            bestpoly = Polynomial(self.taskparams["bestpoly"].splitlines())
            bestpoly.setE(self.taskparams["bestE"])
            self.logger.info("Best polynomial previously found has "
                             "Murphy_E = %g", bestpoly.E)
        else:
            bestpoly = None
        
        while not self.is_done():
            for key in self.paramnames:
                if not key in self.taskparams:
                    raise Exception("Required parameter %s not set for %s" %
                                    (key, self.name))
            args = (self.taskparams["P"],)
            
            admin = int(self.taskparams["admin"])
            admax = admin + int(self.taskparams["adrange"])
            kwargs = self.progparams[0]
            kwargs["admin"] = str(admin)
            kwargs["admax"] = str(admax)
            outputfile = "%s.polysel.%d-%d" % \
               (self.taskparams["name"], admin, admax)
            p = cadoprograms.Polyselect2l(self.taskparams["N"], 
                                          self.taskparams["P"],
                                          outputfile, kwargs)
            p.run()
            p.wait()
            self.taskparams["admin"] = str(admax)
            with open(outputfile, "r") as f:
                poly = Polynomial(f)
            if not bool(poly):
                self.logger.info("No polynomial found in file %s", 
                                 outputfile)
            elif not poly.E:
                raise Exception("Polynomial in file %s has no Murphy E value" 
                                % outputfile)
            elif not bestpoly or poly > bestpoly:
                bestpoly = poly
                self.taskparams["bestE"] = str(poly.E)
                self.taskparams["bestpoly"] = str(poly)
                self.logger.info("New best polynomial from file %s:"
                                 " Murphy E = %g" % (outputfile, poly.E))
                self.logger.info("%s", poly)
            else:
                self.logger.info("Best polynomial from file %s with E=%g is "
                                 "no better than current best with E=%g",
                                 outputfile, poly.E, bestpoly.E)
            # print(poly)

        self.logger.debug("Best polynomial has Murphy_E =  %g", bestpoly.E)
        self.logger.debug("Marking " + self.name + " done")
        self.mark_done(True)
        self.logger.debug("Exit PolyselTask.run(" + self.name + ")")
        return
    


class FactorBaseTask(Task):
    """ Generates the factor base for the polynomial(s) """
    name = "factorbase"
    title = "Generate Factor Base"
    programs = (cadoprograms.MakeFB,)
    paramnames = ()
    
    def __init__(self, polsel, *args, **kwargs):
        super().__init__(*args, dependencies = (polsel,), **kwargs)

    def run(self, parameters = None):
        self.logger.debug("Enter FactorBaseTask.run(" + self.name + ")")
        super().run(parameters = parameters)
        if not self.is_done():
            args = ()
            kwargs = self.progparams[0]
            p = self.programs[0](args, kwargs)
            p.run()
            p.wait()
            self.logger.debug("Marking " + self.name + " done")
            self.mark_done(True)
        self.logger.debug("Exit FactorBaseTask.run(" + self.name + ")")
        
class FreeRelTask(Task):
    """ Generates free relations for the polynomial(s) """
    name = "freerel"
    title = "Generate Free Relations"
    programs = (cadoprograms.FreeRel,)
    paramnames = ()
    
    def __init__(self, polsel, *args, **kwargs):
        super().__init__(*args, dependencies = (polsel,), **kwargs)

    def run(self, parameters = None):
        self.logger.debug("Enter FreeRelTask.run(" + self.name + ")")
        super().run(parameters = parameters)
        if not self.is_done():
            args = ()
            kwargs = self.progparams[0]
            p = self.programs[0](args, kwargs)
            p.run()
            p.wait()
            self.logger.debug("Marking " + self.name + " done")
            self.mark_done(True)
        self.logger.debug("Exit FreeRelTask.run(" + self.name + ")")

class SievingTask(Task):
    """ Does the sieving, uses client/server """
    name = "sieving"
    title = "Sieving"
    programs = (cadoprograms.Las,)
    paramnames = ()
    
    def __init__(self, polsel, factorbase, *args, **kwargs):
        super().__init__(*args, dependencies = (polsel, factorbase), **kwargs)

    def run(self, parameters = None):
        self.logger.debug("Enter SievingTask.run(" + self.name + ")")
        super().run(parameters = parameters)
        if not self.is_done():
            args = ()
            kwargs = self.progparams[0]
            p = self.programs[0](args, kwargs)
            p.run()
            p.wait()
            self.logger.debug("Marking " + self.name + " done")
            self.mark_done(True)
        self.logger.debug("Exit SievingTask.run(" + self.name + ")")

class DuplicatesTask(Task):
    """ Removes duplicate relations """
    name = "duplicates"
    title = "Filtering: Duplicate Removal"
    programs = (cadoprograms.Duplicates1, cadoprograms.Duplicates2)
    paramnames = ()
    
    def __init__(self, sieving, *args, **kwargs):
        super().__init__(*args, dependencies = (sieving,), **kwargs)

    def run(self, parameters = None):
        self.logger.debug("Enter DuplicatesTask.run(" + self.name + ")")
        super().run(parameters = parameters)
        if not self.is_done():
            args = ()
            kwargs = self.progparams[0]
            p = self.programs[0](args, kwargs)
            p.run()
            p.wait()
            self.logger.debug("Marking " + self.name + " done")
            self.mark_done(True)
        self.logger.debug("Exit DuplicatesTask.run(" + self.name + ")")

class PurgeTask(Task):
    """ Removes singletons and computes excess """
    name = "singletons"
    title = "Filtering: Singleton removal"
    programs = (cadoprograms.Purge,)
    paramnames = ()
    
    def __init__(self, duplicates, *args, **kwargs):
        super().__init__(*args, dependencies = (duplicates,), **kwargs)

    def run(self, parameters = None):
        self.logger.debug("Enter PurgeTask.run(" + self.name + ")")
        super().run(parameters = parameters)
        if not self.is_done():
            args = ()
            kwargs = self.progparams[0]
            p = self.programs[0](args, kwargs)
            p.run()
            p.wait()
            self.logger.debug("Marking " + self.name + " done")
            self.mark_done(True)
        self.logger.debug("Exit PurgeTask.run(" + self.name + ")")

class MergeTask(Task):
    """ Merges relations """
    name = "merging"
    title = "Filtering: Merging"
    programs = (cadoprograms.Merge,)
    paramnames = ()
    
    def __init__(self, freerel, unique, *args, **kwargs):
        super().__init__(*args, dependencies = (freerel, unique), **kwargs)

    def run(self, parameters = None):
        self.logger.debug("Enter MergeTask.run(" + self.name + ")")
        super().run(parameters = parameters)
        if not self.is_done():
            args = ()
            kwargs = self.progparams[0]
            p = self.programs[0](args, kwargs)
            p.run()
            p.wait()
            self.logger.debug("Marking " + self.name + " done")
            self.mark_done(True)
        self.logger.debug("Exit MergeTask.run(" + self.name + ")")

class LinAlgTask(Task):
    """ Runs the linear algebra step """
    name = "linalg"
    title = "Linear Algebra"
    programs = (cadoprograms.BWC,)
    paramnames = ()
    
    def __init__(self, merge, *args, **kwargs):
        super().__init__(*args, dependencies = (merge,), **kwargs)

    def run(self, parameters = None):
        self.logger.debug("Enter LinAlgTask.run(" + self.name + ")")
        super().run(parameters = parameters)
        if not self.is_done():
            args = ()
            kwargs = self.progparams[0]
            p = self.programs[0](args, kwargs)
            p.run()
            p.wait()
            self.logger.debug("Marking " + self.name + " done")
            self.mark_done(True)
        self.logger.debug("Exit LinAlgTask.run(" + self.name + ")")

class SqrtTask(Task):
    """ Runs the square root """
    name = "sqrt"
    title = "Square Root"
    programs = (cadoprograms.Sqrt,)
    paramnames = ()
    
    def __init__(self, polsel, freerel, sieving, merge, linalg, 
                 *args, **kwargs):
        dep = (polsel, freerel, sieving, merge)
        super().__init__(*args, dependencies = dep, **kwargs)

    def run(self, parameters = None):
        self.logger.debug("Enter SqrtTask.run(" + self.name + ")")
        super().run(parameters = parameters)
        if not self.is_done():
            args = ()
            kwargs = self.progparams[0]
            p = self.programs[0](args, kwargs)
            p.run()
            p.wait()
            self.logger.debug("Marking " + self.name + " done")
            self.mark_done(True)
        self.logger.debug("Exit SqrtTask.run(" + self.name + ")")


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
        self.sing = PurgeTask(self.dup, *args, db=db, **kwargs)
        self.merge = MergeTask(self.freerel, self.sing, *args, db=db, **kwargs)
        self.linalg = LinAlgTask(self.merge, *args, db=db, **kwargs)
        self.sqrt = SqrtTask(self.polysel, self.freerel, self.sieving, 
                             self.merge, self.linalg, *args, db=db, **kwargs)
    
    def run(self, *args, **kwargs):
        self.sqrt.run(*args, **kwargs)
