#!/usr/bin/env python3
import sqlite3
from datetime import datetime
import re
import os.path
from fractions import gcd
import abc
import random
import time
import wudb
import logging
import patterns
import cadoprograms
import wuserver

# Some parameters are provided by the param file but can change during
# the factorization, like rels_wanted. On one hand, we want automatic 
# updates to parameters to be stored in the DB, otoh, we want to allow
# externally setting new parameters. Need to distinguish between new 
# external parameters that overwrite DB, and old external parameters 
# that don't overwrite. Or maybe two ways to specify external params:
# --defaults which does not overwrite, and --forceparam which does


class FilesCreator(wudb.DbAccess, metaclass=abc.ABCMeta):
    """ A base class for classes that produce a list of output files, with
    some auxiliary information stored with each file (e.g., nr. of relations).
    This info is stored in the form of a DB-backed dictionary, with the file
    name as the key and the auxiliary data as the value.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.output_files = self.make_db_dict(self.make_tablename("outputfiles"))
    
    def add_output_files(self, filenames):
        """ Adds a dict of files to the list of existing output files """
        for (filename, value) in filenames.items():
            if filename in self.output_files:
                raise KeyError("%s already in output files table" % filename)
            self.output_files[filename] = value
    
    def get_output_filenames(self, condition = None):
        """ Return output file names, optionally those that match a condition
        
        If a condition is given, it must be callable with 1 parameter and
        boolean return type; then only those filenames are returned where
        for the auxiliary data s (i.e., the value stored in the dictionary
        with the file name as key) satisfies condition(s) == True.
        """
        if condition is None:
            return list(self.output_files.keys())
        else:
            return [f for (f,s) in self.output_files.items() if condition(s)]
    
    def forget_output_filenames(self, filenames):
        for f in filenames:
            del(self.output_files[f])
    
    @abc.abstractclassmethod
    def make_tablename(self, name):
        pass


class Task(wudb.DbAccess, metaclass=abc.ABCMeta):
    """ A base class that represents one task that needs to be processed. 
    
    Sub-classes must define class variables:
    """
    
    # Properties that subclasses need to define
    @abc.abstractproperty
    def name(self):
        # The name of the task in a simple form that can be used as
        # a Python dictionary key, a directory name, part of a file name,
        # part of an SQL table name, etc. That pretty much limits it to
        # alphabetic first letter, and alphanumeric rest. 
        pass
    @abc.abstractproperty
    def parampath(self):
        pass
    @abc.abstractproperty
    def title(self):
        # A pretty name for the task, will be used in screen output
        pass
    @abc.abstractproperty
    def programs(self):
        # A list of classes of Programs which this tasks uses
        pass
    @abc.abstractproperty
    def paramnames(self):
        # A list of parameter keywords which this task uses.
        # This is used for extracting relevant parameters from the parameter
        # hierarchical dictionary.
        # Sub-classes need to define a property 'paramnames' which returns a
        # list of parameters they accept, plus super()'s paramnames list
        # Parameters that all tasks use
        return ("name", "workdir", "remove")
    
    def __init__(self, dependencies, *args, parameters = None, **kwargs):
        ''' Sets up a database connection and a DB-backed dictionary for 
        parameters. Reads parameters from DB, and merges with hierarchical
        parameters in the parameters argument. Parameters passed in by 
        parameters argument do not override values in the DB-backed 
        parameter dictionary.
        '''
        
        super().__init__(*args, **kwargs)
        if dependencies:
            self.dependencies = dependencies
            # for d in dependencies:
            #    d.subscribeObserver(self)
        else:
            self.dependencies = {}
        self.logger = logging.getLogger(self.title)
        self.logger.debug("Enter Task.__init__(%s)", 
                          self.name)
        if False:
            self.logger.debug("Enter Task.__init__(): parameters = %s", 
                              parameters)
        # DB-backed dictionary with the state of this task
        self.state = self.make_db_dict(self.make_tablename())
        self.logger.debug("%s: state = %s", self.title, self.state)
        # Set default parametes for this task, if any are given
        if parameters:
            self.params = parameters.myparams(self.paramnames, self.parampath)
            self.logger.debug("%s: params = %s", self.title, self.params)
            self.params.setdefault("remove", False)
        # Set default parameters for our programs
        self.progparams = []
        for prog in self.programs:
            if parameters:
                progparams = parameters.myparams(
                    prog.get_config_keys(), [self.parampath, prog.name])
            else:
                progparams = {}
            self.progparams.append(progparams)
        self.logger.debug("Exit Task.__init__(%s)", self.name)
        return
    
    @staticmethod
    def check_tablename(name):
        no_ = name.replace("_", "")
        if not no_[0].isalpha() or not no_[1:].isalnum():
            raise Exception("%s is not valid for an SQL table name" % name)
    
    def make_tablename(self, extra = None):
        """ Return the table name for the DB-backed dictionary with the state
        for the current task """
        # Maybe replace SQL-disallowed characters here, like digits and '.' ? 
        # Could be tricky to avoid collisions
        name = self.name
        if extra:
            name = name + '_' + extra
        self.check_tablename(name)
        return name
    
    def is_done(self):
        for d in self.dependencies:
            if not d.is_done():
                return False
        return True
    
    def run(self):
        ''' Runs the prerequisites. Sub-classes should call this first in 
        their run() method.
                '''
        self.logger.debug("Enter Task.run(%s)" % type(self))
        self.logger.debug("Task.run(%s): self.is_done() = %s", 
                          self.name, self.is_done())
        # Check/run the prerequisites
        for task in self.dependencies:
            if not task.is_done():
                self.logger.debug("Task.run(%s): Running prerequisite %s",
                                  self.name, task.name)
                task.run()
        
        self.logger.debug("Exit Task.run(%s)" % type(self))
        self.logger.info("Starting")
        return
    
    def _make_basename(self):
        """
        >>> class C(object):
        ...     pass
        >>> f = C()
        >>> f.params = {"workdir": "/foo/bar", "name": "jobname"}
        >>> f.name = "taskname"
        >>> Task._make_basename(f)
        '/foo/bar/jobname.taskname'
        """
        return "%s%s%s.%s" % (self.params["workdir"].rstrip(os.sep), os.sep,
                              self.params["name"], self.name)
    
    def make_output_filename(self, name, subdir = False, dirextra = None):
        """ Make a filename of the form workdir/jobname.taskname.name """
        if subdir:
            return self.make_output_dirname(dirextra) + name
        else:
            assert dirextra == None
            return "%s.%s" % (self._make_basename(), name)

    def make_output_dirname(self, extra = None):
        """ Make a directory name of the form workdir/jobname.taskname/ if 
        extra is not given, or workdir/jobname.taskname/extra/ if it is
        """
        r = self._make_basename() + os.sep
        if extra:
            r += "%s%s" % (extra, os.sep)
        return r
    
    def translate_input_filename(self, filename):
        return filename

    def test_outputfile_exists(self, filename):
        return os.path.isfile(filename)
    
    @staticmethod
    def check_files_exist(filenames, filedesc, shouldexist):
        """ Check that the output files in "filenames" exist or don't exist, 
        according to shouldexist.
        
        Raise IOError if any check fails, return None
        """
        for f in filenames:
            exists = os.path.isfile(f)
            if shouldexist and not exists:
                raise IOError("%s file %s does not exist" % (filedesc, f))
            elif not shouldexist and exists:
                raise IOError("%s file %s already exists" % (filedesc, f))
        return
    
    @staticmethod
    def make_directories(basedir, extra = None):
        if not os.path.isdir(basedir):
            os.mkdir(basedir)
        if extra:
            for subdir in extra:
                dirname = basedir + os.sep + subdir
                if not os.path.isdir(dirname):
                    os.mkdir(dirname)
        return

    # These two function go together, one produces a workunit name from the 
    # name of the factorization, the task name, and a task-provided identifier,
    # and the other function splits them again
    wu_paste_char = '_'
    def make_wuname(self, identifier):
        assert not self.wu_paste_char in self.params["name"]
        assert not self.wu_paste_char in self.name
        assert not self.wu_paste_char in identifier
        return self.wu_paste_char.join([self.params["name"], self.name, identifier])
    
    def split_wuname(self, wuname):
        return wuname.split(self.wu_paste_char)
    
    class ResultInfo(wudb.WuResultMessage):
        def __init__(self, wuid, rc, stdout, stderr, output_files):
            self.wuid = wuid
            self.rc = rc
            self.stdout = stdout
            self.stderr = stderr
            self.output_files = output_files
        def get_wu_id(self):
            return self.wuid
        def get_output_files(self):
            return self.output_files
        def get_stdout(self, command_nr):
            assert command_nr == 0
            return self.stdout
        def get_stderr(self, command_nr):
            assert command_nr == 0
            return self.stderr
        def get_exitcode(self, command_nr):
            assert command_nr == 0
            return self.rc
    
    def submit_command(self, command, identifier):
        ''' Run a command.
        Return the result tuple. If the caller is an Observer, also send
        result to updateObserver().
        '''
        wuname = self.make_wuname(identifier)
        command.run()
        (rc, stdout, stderr) = command.wait()
        result = [wuname, rc, stdout, stderr]
        result.append(command.get_output_files())
        if isinstance(self, patterns.Observer):
            message = Task.ResultInfo(wuname, rc, stdout, stderr,
                                      command.get_output_files())
            self.updateObserver(message)
        return result
    
    def filter_notification(self, message):
        wuid = message.get_wu_id()
        rc = message.get_exitcode(0)
        stdout = message.get_stdout(0)
        stderr = message.get_stderr(0)
        output_files = message.get_output_files()
        self.logger.debug("%s: Received notification for wuid=%s, rc=%d, output_files=[%s]",
                         self.name, wuid, rc, ", ".join(output_files))
        if rc != 0:
            self.logger.error("Return code is: %d", rc)
        if stdout:
            self.logger.debug("stdout is: %s", stdout)
        if stderr:
            self.logger.error("stderr is: %s", stderr)
        if output_files:
            self.logger.debug("Output files are: %s", ", ".join(output_files))
        (name, task, identifier) = self.split_wuname(wuid)
        if name != self.params["name"] or task != self.name:
            # This notification is not for me
            self.logger.debug("%s: Notification is not for me", self.name)
            return
        self.logger.debug("%s: Notification is for me", self.name)
        return identifier

class ClientServerTask(Task, patterns.Observer):
    @abc.abstractproperty
    def paramnames(self):
        return super().paramnames + ("maxwu",)
    
    def __init__(self, db_listener, registered_filenames, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.wuar = self.make_wu_access()
        self.wuar.create_tables()
        self.registered_filenames = registered_filenames
        self.db_listener = db_listener
        self.db_listener.subscribeObserver(self)
        self.state.setdefault("wu_submitted", 0)
        self.state.setdefault("wu_received", 0)
        self.params.setdefault("maxwu", 10)
        assert self.state["wu_received"] <= self.state["wu_submitted"]
    
    def submit_command(self, command, identifier):
        ''' Submit a workunit to the database.
        Return the result tuple. If the caller is an Observer, also send
        result to updateObserver().
        '''
        
        while self.state["wu_submitted"] - self.state["wu_received"] >= self.params["maxwu"]:
            self.wait()
        wuid = self.make_wuname(identifier)
        wutext = command.make_wu(wuid)
        for filename in command.get_exec_files() + command.get_input_files():
            basename = os.path.basename(filename)
            if not basename in self.registered_filenames:
                self.logger.debug("Registering file name %s with target %s", basename, filename)
                self.registered_filenames.setdefault(basename, filename)
            # If it was already defined with a different target, error
            if not self.registered_filenames[basename] == filename:
                raise Exception("Filename %s, to be registered for target %s, already registered for target %s" %
                                basename, filename, self.registered_filenames[basename])
        
        self.logger.info("Adding workunit %s to database", wuid)
        # print ("WU:\n%s" % wutext)
        self.wuar.create(wutext)
        self.state["wu_submitted"] += 1
    
    def verification(self, message, ok):
        not_str = "ok"
        wuid = message.get_wu_id()
        if not ok:
            not_str = "not ok"
        self.logger.info("Marking workunit %s as %s", wuid, not_str)
        self.wuar.verification(message.get_wu_id(), ok)
    
    def test_outputfile_exists(self, filename):
        # Can't test
        return False
    
    def wait(self):
        # If we get notification on new results reliably, we might not
        # need this poll. But they probably won't be totally reliable
        if not self.db_listener.send_result():
            time.sleep(1)


class Polynomial(object):
    # Keys that can occur in a polynomial file in their preferred ordering,
    # and whether the key is mandatory or not. The preferred ordering is used
    # when turning a polynomial back into a string. 
    
    paramnames = ("rlim", "alim", "lpbr", "lpba", "mfbr", "mfba", "rlambda", 
        "alambda")
    keys = ( ("n", True), ("Y0", True), ("Y1", True), ("c0", True), 
            ("c1", True), ("c2", True), ("c3", True), ("c4", True),
            ("c5", False), ("c6", False), ("m", True), ("skew", True) )
    
    def __init__(self, lines):
        """ Parse a polynomial file in the syntax as produced by polyselect2l 
        """
        self.poly = None
        self.E = 0.
        poly = {}
        for line in lines:
            # print ("Parsing line: >%s<" % line)
            # If there is a "No polynomial found" message anywhere in the
            # file, we assume that there is no polynomial. This assumption
            # will be false if, e.g., files get concatenated
            if re.match("No polynomial found", line):
                return
            # If this is a comment line telling the Murphy E value, 
            # extract the value and store it
            match = re.match("\s*#\s*MurphyE\s*\(.*\)=(.*)$", line)
            if match:
                self.E = float(match.group(1))
                continue
            # Drop comment, strip whitespace
            l = line.split('#', 1)[0].strip()
            # If nothing is left, process next line
            if not l:
                continue
            # All remaining lines must be of the form "x: y"
            a = l.split(":")
            if not len(a) == 2:
                raise Exception("Invalid line %s" % l)
            key = a[0].strip()
            value = a[1].strip()
            if not key in dict(self.keys):
                raise Exception("Invalid key %s in line %s" %
                                (key, l))
            poly[key] = value
        for (key, isrequired) in self.keys:
            if isrequired and not key in poly:
                raise Exception("Key %s missing" % key)
        self.poly = poly
        return
    
    def __str__(self):
        arr = [(key + ": " + self.poly[key] + '\n')
               for (key,req) in self.keys if key in self.poly]
        return "".join(arr)

    def __eq__(self, other):
        return self.poly == other.poly
    def __ne__(self, other):
        return self.poly != other.poly

    def is_valid(self):
        return not self.poly is None
    
    def setE(self, E):
        self.E = float(E)
    
    def create_file(self, filename, params):
        # Write polynomial to a file, and add lines with parameters such as 
        # "alim" if supplied in params 
        with open(filename, "w") as f:
            f.write(str(self))
            for key in self.paramnames:
                if key in params:
                    f.write(key + ": %s\n" % params[key])


class PolyselTask(ClientServerTask):
    """ Finds a polynomial, uses client/server """
    @property
    def name(self):
        return "polyselect"
    @property
    def parampath(self):
        return "tasks." + self.name
    @property
    def title(self):
        return "Polynomial Selection"
    @property
    def programs(self):
        return (cadoprograms.Polyselect2l,)
    @property
    def paramnames(self):
        return super().paramnames + \
            ("adrange", "P", "N", "admin", "admax")
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, dependencies = None, **kwargs)
        self.state["adnext"] = \
            max(self.state.get("adnext", 0), int(self.params.get("admin", 0)))
    
    def is_done(self):
        # self.logger.debug ("PolyselTask.is_done(): Task parameters: %s", 
        #                    self.params)
        return super().is_done() and \
               "bestpoly" in self.state and \
               self.state["adnext"] >= int(self.params["admax"]) and \
               self.state["wu_received"] == self.state["wu_submitted"]
    
    def run(self):
        # Make command line for polselect2l, run it. 
        self.logger.debug("Enter PolyselTask.run(" + self.name + ")")
        super().run()
        
        self.logger.debug("PolyselTask.run(): Task state: %s", 
                          self.state)
        self.logger.debug("PolyselTask.run(): Task parameters: %s", 
                          self.params)
        
        if "bestpoly" in self.state:
            self.bestpoly = Polynomial(self.state["bestpoly"].splitlines())
            self.bestpoly.setE(self.state["bestE"])
            self.logger.info("Best polynomial previously found in %s has "
                             "Murphy_E = %g", 
                             self.state["bestfile"], self.bestpoly.E)
        else:
            self.bestpoly = None
            self.logger.info("No polynomial was previously found")
        
        while self.state["adnext"] <= int(self.params["admax"]):
            adstart = self.state["adnext"]
            adend = adstart + int(self.params["adrange"])
            polyselect_params = self.progparams[0].copy()
            polyselect_params["admin"] = str(adstart)
            polyselect_params["admax"] = str(adend)
            outputfile = self.make_output_filename("%d-%d" % (adstart, adend))
            if self.test_outputfile_exists(outputfile):
                self.logger.info("%s already exists, won't generate again",
                                 outputfile)
            else:
                p = self.programs[0](None, polyselect_params, stdout = outputfile)
                self.submit_command(p, "%d-%d" % (adstart, adend))
            self.state["adnext"] = adend
        
        while self.state["wu_received"] < self.state["wu_submitted"]:
            self.wait()
        
        self.logger.info("finished")
        if not self.bestpoly:
            self.logger.error ("No polynomial found")
            return
        self.logger.info("Best polynomial from %s has Murphy_E = %g", 
                          self.state["bestfile"] , self.bestpoly.E)
        self.logger.debug("Exit PolyselTask.run(" + self.name + ")")
        return
    
    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return
        output_files = message.get_output_files()
        assert len(output_files) == 1
        self.state["wu_received"] += 1
        ok = self.parse_poly(output_files[0])
        self.verification(message, ok)
    
    def parse_poly(self, outputfile):
        poly = None
        if outputfile:
            with open(outputfile, "r") as f:
                try:
                    poly = Polynomial(f)
                except Exception as e:
                    self.logger.error("Invalid polyselect file %s: %s", 
                                      outputfile, e)
                    return False
        
        if not poly or not poly.is_valid():
            self.logger.info('No polynomial found in %s', outputfile)
            return False
        if not poly.E:
            self.logger.error("Polynomial in file %s has no Murphy E value" 
                              % outputfile)
            return False
        if not self.bestpoly or poly.E > self.bestpoly.E:
            self.bestpoly = poly
            self.state["bestE"] = poly.E
            self.state["bestpoly"] = str(poly)
            self.state["bestfile"] = outputfile
            self.logger.info("New best polynomial from file %s:"
                             " Murphy E = %g" % (outputfile, poly.E))
            self.logger.debug("New best polynomial is:\n%s", poly)
        else:
            self.logger.info("Best polynomial from file %s with E=%g is "
                             "no better than current best with E=%g",
                             outputfile, poly.E, self.bestpoly.E)
        # print(poly)
        return True
    
    def get_poly(self):
        if "bestpoly" in self.state:
            return Polynomial(self.state["bestpoly"].splitlines())
        else:
            return None


class FactorBaseOrFreerelTask(Task):
    """ Common base class for programs that produce one output file from 
    the polynomial, i.e., factorbase and freerel 
    """
    
    def __init__(self, polyselect, *args, **kwargs):
        super().__init__(*args, dependencies = (polyselect,), **kwargs)
        self.polyselect = polyselect
        # Invariant: if we have a result (in self.state["outputfile"]) then we
        # must also have a polynomial (in self.state["poly"])
        if "outputfile" in self.state:
            assert "poly" in self.state
            # The target file must correspond to the polynomial "poly"
    
    def run(self):
        self.logger.debug("Enter FactorBaseOrFreerelTask.run(%s)", self.name)
        super().run()
        
        # Get best polynomial found by polyselect
        poly = self.polyselect.get_poly()
        if not poly:
            raise Exception("FactorBaseOrFreerelTask(): no polynomial "
                            "received from PolyselTask")
        
        # Check if we have already computed the target file for this polynomial
        if "poly" in self.state:
            prevpoly = Polynomial(self.state["poly"].splitlines())
            if poly != prevpoly:
                if "outputfile" in self.state:
                    self.logger.info("Received different polynomial from %s, discarding old one",
                                     self.polyselect.title)
                    del(self.state["outputfile"])
                self.state["poly"] = str(poly)
        else:
            self.state["poly"] = str(poly)
        
        if not self.is_done():
            # Write polynomial to a file
            polyfile = self.make_output_filename("poly")
            poly.create_file(polyfile, self.params)
            
            # Make file name for factor base file
            outputfile = self.make_output_filename(self.target)
            args = ()
            kwargs = self.progparams[0].copy()
            kwargs["poly"] = polyfile
            if "pmin" in self.programs[0].get_config_keys():
                kwargs.setdefault("pmin", "1")
            if "pmax" in self.programs[0].get_config_keys():
                kwargs.setdefault("pmax", str(2**int(self.params["lpba"])))
            p = self.programs[0](None, kwargs, stdout = outputfile)
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            self.parse_stderr(stderr)
            
            self.state["outputfile"] = outputfile
        self.logger.debug("Exit FactorBaseOrFreerelTask.run(%s)", self.name)
    
    def is_done(self):
        if "outputfile" in self.state\
                and not os.path.isfile(self.state["outputfile"]):
            raise Exception("FactorBaseOrFreerelTask.is_done(%s): marked "
                            "as done but target file %s does not exist" % 
                            (self.name, self.state["outputfile"]))
        return super().is_done() and "outputfile" in self.state
    
    def updateObserver(self, message):
        if isinstance(message, Polynomial):
            # We might start factorbase computation while polynomial 
            # selection is still running
            pass
    
    def get_filename(self):
        if "outputfile" in self.state:
            return self.state["outputfile"]
        else:
            return None

    def parse_stderr(self, stderr):
        pass


class FactorBaseTask(FactorBaseOrFreerelTask):
    """ Generates the factor base for the polynomial(s) """
    @property
    def name(self):
        return "factorbase"
    @property
    def title(self):
        return "Generate Factor Base"
    @property
    def programs(self):
        return (cadoprograms.MakeFB,)
    @property
    def parampath(self):
        return "tasks.sieve." + self.name
    @property
    def paramnames(self):
        return super().paramnames + \
            ("alim", )
    target = "roots"


class FreeRelTask(FactorBaseOrFreerelTask):
    """ Generates free relations for the polynomial(s) """
    @property
    def name(self):
        return "freerel"
    @property
    def title(self):
        return "Generate Free Relations"
    @property
    def programs(self):
        return (cadoprograms.FreeRel,)
    @property
    def parampath(self):
        return "tasks.sieve." + self.name
    @property
    def paramnames(self):
        return super().paramnames + \
            ("lpba", )
    target = "freerel"
    
    def run(self):
        super().run()
        self.logger.info("Found %d free relations" % self.state["nfree"])
    
    def parse_stderr(self, stderr):
        if "nfree" in self.state:
            del(self.state["nfree"])
        for line in stderr.decode("ascii").splitlines():
            match = re.match('# Free relations: (\d+)', line)
            if match:
                if "nfree" in self.state:
                    raise Exception("Received two values for number of free relations")
                self.state["nfree"] = int(match.group(1))
        if not "nfree" in self.state:
            raise Exception("Received no value for number of free relations")
        return
    
    def get_nrels(self):
        return self.state["nfree"]
    

class SievingTask(ClientServerTask, FilesCreator):
    """ Does the sieving, uses client/server """
    @property
    def name(self):
        return "sieving"
    @property
    def title(self):
        return "Lattice Sieving"
    @property
    def programs(self):
        return (cadoprograms.Las,)
    @property
    def parampath(self):
        return "tasks.sieve." + self.name
    @property
    def paramnames(self):
        return super().paramnames + \
            ("qmin", "qrange", "rels_wanted") + Polynomial.paramnames
    # We seek to this many bytes before the EOF to look for the "Total xxx reports" message
    file_end_offset = 1000
    
    def __init__(self, polyselect, factorbase, *args, **kwargs):
        super().__init__(*args, dependencies = (polyselect, factorbase), **kwargs)
        self.polyselect = polyselect
        self.factorbase = factorbase
        # qmin is optional, but if it exists, should be use in preference to alim
        if "qmin" in self.params:
            self.state.setdefault("qnext", int(self.params["qmin"]))
        self.state.setdefault("qnext", int(self.params["alim"]))
        self.state.setdefault("rels_found", 0)
        self.state.setdefault("rels_wanted", 0)
        self.params.setdefault("max_wus", 10)
        if "rels_wanted" in self.params:
            self.state["rels_wanted"] = max(self.state["rels_wanted"], int(self.params["rels_wanted"]))
        else:
            # Choose sensible default value
            pass
    
    def run(self):
        self.logger.debug("Enter SievingTask.run(" + self.name + ")")
        super().run()
        
        # Get best polynomial found by polyselect
        poly = self.polyselect.get_poly()
        if not poly:
            raise Exception("SievingTask(): no polynomial received")
        # Write polynomial to a file
        polyfile = self.make_output_filename("poly")
        poly.create_file(polyfile, self.params)
        
        while not self.is_done():
            args = ()
            kwargs = self.progparams[0].copy()
            q0 = self.state["qnext"]
            q1 = q0 + int(self.params["qrange"])
            outputfile = self.make_output_filename("%d-%d" % (q0, q1))
            self.check_files_exist([outputfile], "output", shouldexist=False)
            kwargs["q0"] = str(q0)
            kwargs["q1"] = str(q1)
            kwargs["poly"] = polyfile
            kwargs["factorbase"] = self.factorbase.get_filename()
            kwargs["out"] = outputfile
            p = self.programs[0](args, kwargs)
            self.submit_command(p, "%d-%d" % (q0, q1))
            self.state["qnext"] = q1
        self.logger.debug("Exit SievingTask.run(" + self.name + ")")
        return
    
    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return
        self.state["wu_received"] += 1
        output_files = message.get_output_files()
        assert len(output_files) == 1
        ok = self.parse_output_file(output_files[0])
        self.verification(message, ok)
    
    def parse_output_file(self, filename):
        size = os.path.getsize(filename)
        with open(filename, "r") as f:
            f.seek(max(size - self.file_end_offset, 0))
            for line in f:
                match = re.match("# Total (\d+) reports ", line)
                if match:
                    rels = int(match.group(1))
                    self.add_output_files({filename: rels})
                    self.state["rels_found"] += rels
                    self.logger.info("Found %d relations in %s, total is now %d",
                                     rels, filename, self.state["rels_found"])
                    return True
        self.logger.error("Number of relations message not found in file %s", filename)
        return False
    
    def get_nrels(self, filename = None):
        """ Return the number of relations found, either the total so far or
        for a given file
        """
        if filename == None:
            return self.state["rels_found"]
        else:
            return self.output_files[filename]

    def request_more_relations(self, additional):
        self.state["rels_wanted"] += additional
        self.logger.info("%s: New goal for number of relations is %d",
                         self.title, self.state["rels_wanted"])
    
    def is_done(self):
        return super().is_done() and \
               self.state["rels_found"] >= int(self.state["rels_wanted"])
    

class Duplicates1Task(Task, FilesCreator):
    """ Removes duplicate relations """
    @property
    def name(self):
        return "duplicates1"
    @property
    def title(self):
        return "Filtering: Duplicate Removal, splitting pass"
    @property
    def programs(self):
        return (cadoprograms.Duplicates1,)
    @property
    def parampath(self):
        return "tasks.filter." + self.name
    @property
    def paramnames(self):
        return super().paramnames + \
            ("nslices_log",)
    
    def __init__(self, sieving, *args, **kwargs):
        super().__init__(*args, dependencies = (sieving,), **kwargs)
        self.sieving = sieving
        self.nr_slices = 2**int(self.params["nslices_log"])
        self.already_split_input = self.make_db_dict(self.make_tablename("infiles"))
        self.slice_relcounts = self.make_db_dict(self.make_tablename("counts"))
        # Default slice counts to 0, in single DB commit
        self.slice_relcounts.setdefault(
            None, {str(i): 0 for i in range(0, self.nr_slices)})
    
    def run(self):
        self.logger.debug("Enter Duplicates1Task.run(" + self.name + ")")
        super().run()
        # Check that previously split files were split into the same number
        # of pieces that we want now
        for (infile, parts) in self.already_split_input.items():
            if int(parts) != self.nr_slices:
                # TODO: ask interactively (or by -recover) whether to delete 
                # old output files and generate new ones, if input file is 
                # still available
                # If input file is not available but the previously split
                # parts are, we could join them again... not sure if want
                raise Exception("%s was previously split into %d parts, "
                                "now %d parts requested",
                                infile, int(parts), self.nr_slices)
            for outfile in self.make_output_filenames(infile):
                if not outfile in self.output_files:
                    # TODO: How to recover from this error? Just re-split again?
                    raise Exception("Output file %s missing in database for "
                                    "supposedly split input file %s" %
                                    (outfile, infile))
        
        # Check that previously split files do, in fact, exist.
        # FIXME: Do we want this? It may be slow when there are many files.
        # Reading the directory and comparing the lists would probably be
        # faster than individual lookups.
        self.check_files_exist(self.get_output_filenames(), "output",
                               shouldexist=True)
        
        files = self.sieving.get_output_filenames()
        newfiles = [f for f in files if not f in self.already_split_input]
        self.logger.debug ("%s: new files to split are: %s", 
                           self.title, newfiles)
        
        if newfiles:
            basedir = self.make_output_dirname()
            self.make_directories(basedir, map(str, range(0, self.nr_slices)))
            # TODO: can we recover from missing input files? Ask Sieving to
            # generate them again? Just ignore the missing ones?
            self.check_files_exist(newfiles, "input", shouldexist = True)
            # Split the new files
            for f in newfiles:
                outfilenames = self.make_output_filenames(f)
                if self.nr_slices == 1:
                    # If we should split into only 1 part, we don't actually
                    # split at all. We simply write the input file name
                    # to the table of output files, so the next stages will 
                    # read the original siever output file, thus avoiding 
                    # having another copy of the data on disk. Since we don't
                    # process the file at all, we need to ask the Siever task
                    # for the relation count in this file
                    current_counts = [self.sieving.get_nrels(f)]
                else:
                    # TODO: how to recover from existing output files?
                    # Simply re-split? Check whether they all exist and assume 
                    # they are correct if they do?
                    self.check_files_exist(outfilenames.keys(), "output", 
                                            shouldexist=False)
                    kwargs = self.progparams[0].copy()
                    kwargs["out"] = self.make_output_dirname()
                    p = self.programs[0]((f,), kwargs)
                    (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
                    # Check that the output files exist now
                    # TODO: How to recover from error? Presumably a dup1
                    # process failed, but that should raise a return code
                    # exception
                    self.check_files_exist(outfilenames.keys(), "output", 
                                           shouldexist=True)
                    current_counts = self.parse_slice_counts(stderr)
                for (idx, count) in enumerate(current_counts):
                    self.slice_relcounts[str(idx)] += count
                self.already_split_input[f] = self.nr_slices
                self.add_output_files(outfilenames)
        totals = ["%d: %d" % (i, self.slice_relcounts[str(i)])
                  for i in range(0, self.nr_slices)]
        self.logger.info("Relations per slice: %s", ", ".join(totals))
        self.logger.debug("Exit Duplicates1Task.run(" + self.name + ")")
        return
    
    def make_output_filenames(self, name):
        """ Make a dictionary of the output file names corresponding to the
        input file named "name" as keys, and the slice number as a string
        as value. If nr_slices == 1, return the input file name, as in that
        case we do not split at all - we just pass the original file to later
        stages.
        """
        if self.nr_slices == 1:
            return {name: 0}
        else:
            r = {}
            basename = os.path.basename(name)
            for i in range(0, self.nr_slices):
                r[self.make_output_filename(basename, True, str(i))] = i
            return r;
    
    def parse_slice_counts(self, stderr):
        """ Takes line of text and looks for slice counts as printed by dup1
        """
        counts = [None] * self.nr_slices
        for line in stderr.decode("ascii").splitlines():
            match = re.match('# slice (\d+) received (\d+) relations', line)
            if match:
                (slicenr, nrrels) = map(int, match.groups())
                if not counts[slicenr] is None:
                    raise Exception("Received two values for relation count "
                                    "in slice %d" % slicenr)
                counts[slicenr] = nrrels
        for (slicenr, nrrels) in enumerate(counts):
            if nrrels is None:
                raise Exception("Received no value for relation count in "
                                "slice %d" % slicenr)
        return counts
    
    def get_nr_slices(self):
        return self.nr_slices
    
    def get_slice_relcount(self, idx):
        return self.slice_relcounts[str(idx)]
    

class Duplicates2Task(Task, FilesCreator):
    """ Removes duplicate relations """
    @property
    def name(self):
        return "duplicates2"
    @property
    def title(self):
        return "Filtering: Duplicate Removal, removal pass"
    @property
    def programs(self):
        return (cadoprograms.Duplicates2,)
    @property
    def parampath(self):
        return "tasks.filter." + self.name
    @property
    def paramnames(self):
        return super().paramnames + \
            ("nslices_log",)
    
    def __init__(self, duplicates1, *args, **kwargs):
        super().__init__(*args, dependencies = (duplicates1,), **kwargs)
        self.duplicates1 = duplicates1
        self.nr_slices = 2**int(self.params["nslices_log"])
        self.already_done_input = self.make_db_dict(self.make_tablename("infiles"))
        self.slice_relcounts = self.make_db_dict(self.make_tablename("counts"))
        self.slice_relcounts.setdefault(
            None, {str(i): 0 for i in range(0, self.nr_slices)})
    
    def run(self):
        self.logger.debug("Enter Duplicates2Task.run(" + self.name + ")")
        super().run()
        
        for i in range(0, self.nr_slices):
            files = self.duplicates1.get_output_filenames(i.__eq__)
            newfiles = [f for f in files if not f in self.already_done_input]
            if not newfiles:
                self.logger.info("No new files for slice %d, nothing to do", i)
                continue
            # If there are any new files in a slice, we remove duplicates on
            # the whole file set of the slice, as we currently cannot store
            # the duplicate removal state to be able to add more relations
            # in another pass
            # Forget about the previous output filenames of this slice
            # FIXME: Should we delete the files, too?
            self.forget_output_filenames(self.get_output_filenames(i.__eq__))
            del(self.slice_relcounts[str(i)])
            rel_count = self.duplicates1.get_slice_relcount(i)
            basedir = self.make_output_dirname()
            self.make_directories(basedir, str(i))
            kwargs = self.progparams[0].copy()
            kwargs["rel_count"] = str(rel_count * 12 // 10)
            kwargs["output_directory"] = self.make_output_dirname(str(i))
            p = self.programs[0](files, kwargs)
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            nr_rels = self.parse_remaining(stderr.decode("ascii").splitlines())
            # Mark input file names and output file names
            for f in files:
                self.already_done_input[f] = i
            outfilenames = {self.make_output_filename(f, i):i for f in files}
            self.add_output_files(outfilenames)
            self.logger.info("%d unique relations remain on slice %d", nr_rels, i)
            self.slice_relcounts[str(i)] = nr_rels
        self.logger.info("%d unique relations remain in total", self.get_relcount())
        self.logger.debug("Exit Duplicates2Task.run(" + self.name + ")")
    
    def make_output_filename(self, f, i):
        basename = os.path.basename(f)
        return super().make_output_filename(basename, True, str(i))
    
    def parse_remaining(self, text):
        # "     112889 remaining relations"
        for line in text:
            match = re.match('\s*(\d+) remaining relations', line)
            if match:
                return int(match.group(1))
        raise Exception("Received no value for remaining relation count")

    def get_relcount(self):
        nrels = 0
        for i in range(0, self.nr_slices):
            nrels += self.slice_relcounts[str(i)]
        return nrels


class PurgeTask(Task):
    """ Removes singletons and computes excess """
    @property
    def name(self):
        return "singletons"
    @property
    def title(self):
        return "Filtering: Singleton removal"
    @property
    def programs(self):
        return (cadoprograms.Purge,)
    @property
    def parampath(self):
        return "tasks.filter." + self.name
    @property
    def paramnames(self):
        return super().paramnames + \
            ("keep", ) + Polynomial.paramnames
    
    # purge -poly c59.poly -keep 160 -nrels 226167 -out c59.purged.gz -npthr 1x1 -basepath /tmp/cado.wNhDsCM6BS -subdirlist c59.subdirlist -filelist c59.filelist
    
    def __init__(self, polyselect, freerel, sieving, duplicates2, *args, **kwargs):
        super().__init__(*args, dependencies = (duplicates2,), **kwargs)
        self.polyselect = polyselect
        self.freerel = freerel
        self.sieving = sieving
        self.duplicates2 = duplicates2
    
    def run(self):
        self.logger.debug("Enter PurgeTask.run(" + self.name + ")")
        if not self.is_done():
            super().run()
            poly = self.polyselect.get_poly()
            polyfile = self.make_output_filename("poly")
            poly.create_file(polyfile, self.params)
            nfree = self.freerel.get_nrels()
            nunique = self.duplicates2.get_relcount()
            nrels = nfree + nunique
            self.logger.info("Reading %d unique and %d free relations, total %d"
                             % (nunique, nfree, nrels))
            if not nrels:
                raise Exception("No relation count received from %s", self.duplicates2.title)
            purgedfile = self.make_output_filename("purged.gz")
            args = self.duplicates2.get_output_filenames() + [self.freerel.get_filename()]
            kwargs = self.progparams[0].copy()
            kwargs["poly"] = polyfile
            kwargs["keep"] = self.params["keep"]
            kwargs["nrels"] = str(nrels)
            kwargs["out"] = purgedfile
            p = self.programs[0](args, kwargs)
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            if self.parse_stderr(stderr):
                [nrows, weight, excess] = self.parse_stdout(stdout)
                self.logger.info("After purge, %d relations remain with weight %s and excess %s"
                                 % (nrows, weight, excess))
                self.state["purgedfile"] = purgedfile
            else:
                self.sieving.request_more_relations(int(self.sieving.get_nrels() * 0.1))
        self.logger.debug("Exit PurgeTask.run(" + self.name + ")")
    
    def get_purged_filename(self):
        return self.state["purgedfile"]
    
    def is_done(self):
        return super().is_done() and "purgedfile" in self.state
    
    def parse_stderr(self, stderr):
        # If stderr ends with 
        # b'excess < 0.10 * #primes. See -required_excess argument.'
        # then we need more relations from filtering and return False
        for line in stderr.decode("ascii").splitlines():
            if re.match("excess < \d+.\d+ \* #primes", line) or \
                    re.match("number of relations <= number of ideals", line):
                self.logger.info("%s: not enough relations" % self.title)
                return False
        return True
    
    def parse_stdout(self, stdout):
        # Program output is expected in the form:
        # b'NROWS:27372 WEIGHT:406777 WEIGHT*NROWS=1.11e+10\nEXCESS: 160\n'
        # but we allow some extra whitespace
        # If output ends with 
        # b'excess < 0.10 * #primes. See -required_excess argument.'
        # then we need more relations from filtering and return False
        r = {}
        keys = ("NROWS", "WEIGHT", "EXCESS")
        for line in stdout.decode("ascii").splitlines():
            for key in keys:
                match = re.search("%s\s*:\s*(\d+)" % key, line)
                if match:
                    r[key] = int(match.group(1))
        for key in keys:
            if not key in r:
                raise Exception("%s: output of %s did not contain value for %s: %s"
                                % (self.title, self.programs[0].name, key, stdout))
        return [r[key] for key in keys]

class MergeTask(Task):
    """ Merges relations """
    @property
    def name(self):
        return "merging"
    @property
    def title(self):
        return "Filtering: Merging"
    @property
    def programs(self):
        return (cadoprograms.Merge, cadoprograms.Replay)
    @property
    def parampath(self):
        return "tasks.filter." + self.name
    @property
    def paramnames(self):
        return super().paramnames + \
            ("skip", "forbw", "coverNmax", "keep", "maxlevel", "ratio")
    
    def __init__(self, purged, *args, **kwargs):
        super().__init__(*args, dependencies = (purged, ), **kwargs)
        self.purged = purged

    def run(self):
        # merge -out c59.merge.his -mat c59.purged.gz -skip 32 -forbw 3 -coverNmax 100 -keep 160 -maxlevel 15 -ratio 1.1
        # replay --binary -skip 32 -purged c59.purged.gz -his c59.merge.his -index c59.index -out c59.small.bin
        self.logger.debug("Enter MergeTask.run(" + self.name + ")")
        if not self.is_done():
            super().run()
            if "indexfile" in self.state:
                del(self.state["indexfile"])
            if "mergedfile" in self.state:
                del(self.state["mergedfile"])
            if "densefile" in self.state:
                del(self.state["densefile"])
            
            historyfile = self.make_output_filename("history")
            args = ()
            kwargs = self.progparams[0].copy()
            kwargs["mat"] = self.purged.get_purged_filename()
            kwargs["out"] = historyfile
            p = self.programs[0](args, kwargs)
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            
            indexfile = self.make_output_filename("index")
            mergedfile = self.make_output_filename("small.bin")
            args = ()
            kwargs = self.progparams[1].copy()
            kwargs["binary"] = "1"
            kwargs["purged"] = self.purged.get_purged_filename()
            kwargs["history"] = historyfile
            kwargs["index"] = indexfile
            kwargs["out"] = mergedfile
            p = self.programs[1](args, kwargs)
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            
            if not os.path.isfile(indexfile):
                raise Exception("Output file %s does not exist" % indexfile)
            if not os.path.isfile(mergedfile):
                raise Exception("Output file %s does not exist" % mergedfile)
            self.state["indexfile"] = indexfile
            self.state["mergedfile"] = mergedfile
            densefilename = self.make_output_filename("small.dense.bin")
            if os.path.isfile(densefilename):
                self.state["densefile"] = densefilename
            
        self.logger.debug("Exit MergeTask.run(" + self.name + ")")
    
    def get_index_filename(self):
        return self.state.get("indexfile", None)
    
    def get_merged_filename(self):
        return self.state.get("mergedfile", None)
    
    def get_dense_filename(self):
        return self.state.get("densefile", None)
    
    def is_done(self):
        return super().is_done() and "mergedfile" in self.state


class LinAlgTask(Task):
    """ Runs the linear algebra step """
    @property
    def name(self):
        return "bwc"
    @property
    def title(self):
        return "Linear Algebra"
    @property
    def programs(self):
        return (cadoprograms.BWC,)
    @property
    def parampath(self):
        return "tasks.linalg." + self.name
    @property
    def paramnames(self):
        return super().paramnames + \
            ()
    # bwc.pl :complete seed=1 thr=1x1 mpi=1x1 matrix=c59.small.bin nullspace=left mm_impl=bucket interleaving=0 interval=100 mn=64 wdir=c59.bwc shuffled_product=1 bwc_bindir=/localdisk/kruppaal/build/cado-nfs/normal/linalg/bwc
    def __init__(self, merge, *args, **kwargs):
        super().__init__(*args, dependencies = (merge,), **kwargs)
        self.merge = merge
    
    def run(self):
        self.logger.debug("Enter LinAlgTask.run(" + self.name + ")")
        if not self.is_done():
            super().run()
            workdir = self.make_output_dirname()
            self.make_directories(workdir)
            mergedfile = self.merge.get_merged_filename()
            args = ()
            kwargs = self.progparams[0].copy()
            kwargs["complete"] = "1"
            kwargs["matrix"] = os.path.realpath(mergedfile)
            kwargs["wdir"] = os.path.realpath(workdir)
            kwargs.setdefault("nullspace", "left")
            p = self.programs[0](args, kwargs)
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            dependencyfilename = self.make_output_filename("W", subdir = True)
            if not os.path.isfile(dependencyfilename):
                raise Exception("Kernel file %s does not exist" % dependencyfilename)
            self.state["dependency"] =  dependencyfilename
        self.logger.debug("Exit LinAlgTask.run(" + self.name + ")")

    def is_done(self):
        return super().is_done() and "dependency" in self.state
    
    def get_dependency_filename(self):
        return self.state.get("dependency", None)
    
    def get_prefix(self):
        return "%s%s%s.%s" % (self.params["workdir"].rstrip(os.sep), os.sep,
                              self.params["name"], "dep")
    
    
class CharactersTask(Task):
    """ Computes Quadratic Characters """
    @property
    def name(self):
        return "characters"
    @property
    def title(self):
        return "Quadratic Characters"
    @property
    def programs(self):
        return (cadoprograms.Characters,)
    @property
    def parampath(self):
        return "tasks.linalg." + self.name
    @property
    def paramnames(self):
        return super().paramnames + \
            Polynomial.paramnames

    def __init__(self, polyselect, purge, merge, linalg, *args, **kwargs):
        self.polyselect = polyselect
        self.purge = purge
        self.merge = merge
        self.linalg = linalg
        super().__init__(*args, dependencies = (polyselect, purge, merge, linalg), **kwargs)
    
    def run(self):
        self.logger.debug("Enter CharactersTask.run(" + self.name + ")")
        if not self.is_done():
            super().run()
            polyfilename = self.make_output_filename("poly")
            poly = self.polyselect.get_poly()
            poly.create_file(polyfilename, self.params)
            kernelfilename = self.make_output_filename("kernel")
            
            purgedfilename = self.purge.get_purged_filename()
            indexfilename = self.merge.get_index_filename()
            densefilename = self.merge.get_dense_filename()
            dependencyfilename = self.linalg.get_dependency_filename()
            
            args = ()
            kwargs = self.progparams[0].copy()
            kwargs["poly"] = polyfilename
            kwargs["purged"] = purgedfilename
            kwargs["index"] = indexfilename
            kwargs["wfile"] = dependencyfilename
            kwargs["out"] = kernelfilename
            if not densefilename is None:
                kwargs["heavyblock"] = densefilename
            p = self.programs[0](args, kwargs)
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            if not os.path.isfile(kernelfilename):
                raise Exception("Output file %s does not exist" % kernelfilename)
            self.state["kernel"] = kernelfilename
        self.logger.debug("Exit CharactersTask.run(" + self.name + ")")
        return
    
    def is_done(self):
        return super().is_done() and "kernel" in self.state
    
    def get_kernel_filename(self):
        return self.state.get("kernel")


class SqrtTask(Task):
    """ Runs the square root """
    @property
    def name(self):
        return "sqrt"
    @property
    def title(self):
        return "Square Root"
    @property
    def programs(self):
        return (cadoprograms.Sqrt,)
    @property
    def parampath(self):
        return "tasks." + self.name
    @property
    def paramnames(self):
        return super().paramnames + \
            ("N",)
    
    def __init__(self, polyselect, freerel, purge, merge, linalg, characters, *args, **kwargs):
        self.polyselect = polyselect
        self.freerel = freerel
        self.purge = purge
        self.merge = merge
        self.linalg = linalg
        self.characters = characters
        dep = (polyselect, freerel, purge, merge, linalg, characters)
        super().__init__(*args, dependencies = dep, **kwargs)
        self.factors = self.make_db_dict(self.make_tablename("factors"))
        self.add_factor(int(self.params["N"]))
    
    def run(self):
        self.logger.debug("Enter SqrtTask.run(" + self.name + ")")
        if not self.is_done():
            super().run()
            polyfilename = self.make_output_filename("poly")
            poly = self.polyselect.get_poly()
            poly.create_file(polyfilename, self.params)
            purgedfilename = self.purge.get_purged_filename()
            indexfilename = self.merge.get_index_filename()
            kernelfilename = self.characters.get_kernel_filename()
            prefix = self.linalg.get_prefix()
            args = ()
            kwargs = self.progparams[0].copy()
            kwargs["ab"] = "1"
            kwargs["poly"] = polyfilename
            kwargs["purged"] = purgedfilename
            kwargs["index"] = indexfilename
            kwargs["kernel"] = kernelfilename
            kwargs["prefix"] = prefix
            p = self.programs[0](args, kwargs)
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            
            while not self.is_done():
                self.state.setdefault("next_dep", 0)
                kwargs["ab"] = "0"
                kwargs["rat"] = "1"
                kwargs["alg"] = "1"
                kwargs["gcd"] = "1"
                kwargs["dep"] = str(self.state["next_dep"])
                p = self.programs[0](args, kwargs)
                p.run()
                (rc, stdout, stderr) = p.wait()
                if not stdout.decode("ascii").strip() == "Failed":
                    factorlist = list(map(int,stdout.decode("ascii").split()))
                    # FIXME: Can sqrt print more/less than 2 factors?
                    assert len(factorlist) == 2
                    self.add_factor(factorlist[0])
                self.state["next_dep"] += 1
            self.logger.info("%s has finished" % self.title)
        self.logger.info("Factors: %s" % " ".join(self.get_factors()))
        self.logger.debug("Exit SqrtTask.run(" + self.name + ")")
    
    def is_done(self):
        for (factor, isprime) in self.factors.items():
            if not isprime:
                return False
        return super().is_done()
    
    def add_factor(self, factor):
        assert factor > 0
        if str(factor) in self.factors:
            return
        for oldfac in list(map(int, self.factors.keys())):
            g = gcd(factor, oldfac)
            if 1 < g and g < factor:
                self.add_factor(g)
                self.add_factor(factor // g)
                break
            if 1 < g and g < oldfac:
                # We get here only if newfac is a proper factor of oldfac
                assert factor == g
                del(self.factors[str(oldfac)])
                self.add_factor(g)
                self.add_factor(oldfac // g)
                break
        else:
            # We get here if the new factor is coprime to all previously
            # known factors
            isprime = SqrtTask.miller_rabin_tests(factor, 10)
            self.factors[str(factor)] = isprime
    
    def get_factors(self):
        return self.factors.keys()
    
    
    @staticmethod
    def miller_rabin_pass(number, base):
        """
        >>> SqrtTask.miller_rabin_pass(3, 2)
        True
        >>> SqrtTask.miller_rabin_pass(9, 2)
        False
        >>> SqrtTask.miller_rabin_pass(91, 2)
        False
        >>> SqrtTask.miller_rabin_pass(1009, 2)
        True
        >>> SqrtTask.miller_rabin_pass(10000000019, 2)
        True
        >>> SqrtTask.miller_rabin_pass(10000000019*10000000021, 2)
        False
        
        # Check some pseudoprimes. First a few Fermat pseudoprimes which
        # Miller-Rabin should recognize as composite
        >>> SqrtTask.miller_rabin_pass(341, 2)
        False
        >>> SqrtTask.miller_rabin_pass(561, 2)
        False
        >>> SqrtTask.miller_rabin_pass(645, 2)
        False
        
        # Now some strong pseudo-primes
        >>> SqrtTask.miller_rabin_pass(2047, 2)
        True
        >>> SqrtTask.miller_rabin_pass(703, 3)
        True
        >>> SqrtTask.miller_rabin_pass(781, 5)
        True
        """
        if number <= 3:
            return number >= 2
        if number % 2 == 0:
            return False
        # random.randrange(n) produces random integer in [0, n-1]. We want [2, n-2]
        po2 = 0
        exponent = number - 1
        while exponent % 2 == 0:
            exponent >>= 1
            po2 += 1
        
        result = pow(base, exponent, number)
        if result == 1:
            return True
        for i in range(0, po2 - 1):
            if result == number - 1:
                return True
            result = pow(result, 2, number)
        return result == number - 1
    
    @staticmethod
    def miller_rabin_tests(number, passes):
        for i in range(0, passes):
            base = random.randrange(number - 3) + 2
            if not SqrtTask.miller_rabin_pass(number, base):
                return False
        return True

# FIXME: Is this a Task object? Probably not
# Should this be in cadotask or in cadofactor? Probably cadotask
class CompleteFactorization(wudb.DbAccess):
    """ The complete factorization, aggregate of the individual tasks """
    """ Runs the square root """
    @property
    def name(self):
        return "factorization"
    @property
    def title(self):
        return "Factorization"
    @property
    def programs(self):
        return ()
    @property
    def parampath(self):
        return "tasks"
    @property
    def paramnames(self):
        return super().paramnames
    
    def __init__ (self, parameters, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.params = parameters.myparams(("name", "workdir"), "tasks")
        self.db_listener = self.make_db_listener()
        registered_filenames = self.make_db_dict(self.params["name"] + '_server_registered_filenames')
        
        # Set up WU server
        uploaddir = self.params["workdir"].rstrip(os.sep) + os.sep + self.params["name"] + ".upload/"
        self.server = wuserver.ServerLauncher("localhost", 8001, False, self.get_db_filename(), 
            registered_filenames, uploaddir, bg = True)
        self.server.serve()
        
        self.polysel = PolyselTask(*args, parameters = parameters, db_listener = self.db_listener,
                                   registered_filenames = registered_filenames, **kwargs)
        self.fb = FactorBaseTask(self.polysel, *args, parameters = parameters, **kwargs)
        self.freerel = FreeRelTask(self.polysel, *args, parameters = parameters, **kwargs)
        self.sieving = SievingTask(self.polysel, self.fb, *args, parameters = parameters,
                                   db_listener = self.db_listener, registered_filenames = registered_filenames, **kwargs)
        self.dup1 = Duplicates1Task(self.sieving, *args, parameters = parameters, **kwargs)
        self.dup2 = Duplicates2Task(self.dup1, *args, parameters = parameters, **kwargs)
        self.purge = PurgeTask(self.polysel, self.freerel, self.sieving, self.dup2,
                               *args, parameters = parameters, **kwargs)
        self.merge = MergeTask(self.purge, *args, parameters = parameters, **kwargs)
        self.linalg = LinAlgTask(self.merge, *args, parameters = parameters, **kwargs)
        self.characters = CharactersTask(self.polysel, self.purge, self.merge,
                                         self.linalg, *args, parameters = parameters, **kwargs)
        self.sqrt = SqrtTask(self.polysel, self.freerel, self.purge,
                             self.merge, self.linalg, self.characters,
                             *args, parameters = parameters, **kwargs)
    
    def run(self, *args, **kwargs):
        self.sqrt.run(*args, **kwargs)
        self.server.shutdown()
