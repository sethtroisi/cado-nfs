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
import cadoparams
import cadocommand
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


class Task(patterns.Colleague, wudb.DbAccess, cadoparams.UseParameters, metaclass=abc.ABCMeta):
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
        return ("name", "workdir")
    
    def __init__(self, *args, **kwargs):
        ''' Sets up a database connection and a DB-backed dictionary for 
        parameters. Reads parameters from DB, and merges with hierarchical
        parameters in the parameters argument. Parameters passed in by 
        parameters argument do not override values in the DB-backed 
        parameter dictionary.
        '''
        
        super().__init__(*args, **kwargs)
        self.logger = logging.getLogger(self.title)
        self.logger.debug("Enter Task.__init__(%s)", 
                          self.name)
        # DB-backed dictionary with the state of this task
        self.state = self.make_db_dict(self.make_tablename())
        self.logger.debug("state = %s", self.state)
        # Set default parametes for this task, if any are given
        self.params = self.myparams(self.paramnames)
        self.logger.debug("params = %s", self.params)
        # Set default parameters for our programs
        self.progparams = []
        for prog in self.programs:
            progparams = self.myparams(prog.get_config_keys(), prog.name)
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
            assert dirextra is None
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
        process = cadocommand.Command(command)
        (rc, stdout, stderr) = process.wait()
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
    
    def send_notification(self, key, value):
        """ Wrapper around Colleague.send_notification() that instantiates a
        Notification with self as the sender
        """
        notification = Notification(self, key, value)
        super().send_notification(notification)
    
    def send_request(self, key, *args):
        """ Wrapper around Colleague.send_request() that instantiates a
        Request with self as the sender
        """
        request = Request(self, key, *args)
        return super().send_request(request)


class ClientServerTask(Task, patterns.Observer):
    @abc.abstractproperty
    def paramnames(self):
        return super().paramnames + ("maxwu",)
    
    def __init__(self, db_listener, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.db_listener = db_listener
        self.db_listener.subscribeObserver(self)
        self.state.setdefault("wu_submitted", 0)
        self.state.setdefault("wu_received", 0)
        self.params.setdefault("maxwu", 10)
        assert self.get_number_outstanding_wus() >= 0
    
    def submit_command(self, command, identifier):
        ''' Submit a workunit to the database.
        Return the result tuple. If the caller is an Observer, also send
        result to updateObserver().
        '''
        
        while self.get_number_outstanding_wus() >= self.params["maxwu"]:
            self.wait()
        wuid = self.make_wuname(identifier)
        wutext = command.make_wu(wuid)
        for filename in command.get_exec_files() + command.get_input_files():
            basename = os.path.basename(filename)
            self.send_notification(Notification.REGISTER_FILENAME, {basename:filename})
        
        self.logger.info("Adding workunit %s to database", wuid)
        self.send_notification(Notification.SUBMIT_WU, wutext)
        self.state["wu_submitted"] += 1
    
    def verification(self, message, ok):
        wuid = message.get_wu_id()
        not_str = "ok"
        if not ok:
            not_str = "not ok"
        self.logger.info("Marking workunit %s as %s", wuid, not_str)
        assert self.get_number_outstanding_wus() >= 1
        # FIXME: these two should be updated atomically
        self.state["wu_received"] += 1
        self.send_notification(Notification.VERIFY_WU, {message.get_wu_id(): ok})
    
    def get_number_outstanding_wus(self):
        return self.state["wu_submitted"] - self.state["wu_received"]

    def test_outputfile_exists(self, filename):
        # Can't test
        return False
    
    def wait(self):
        # If we get notification on new results reliably, we might not
        # need this poll. But they probably won't be totally reliable
        if not self.db_listener.send_result():
            time.sleep(1)
 

class PolyselTask(ClientServerTask):
    """ Finds a polynomial, uses client/server """
    @property
    def name(self):
        return "polyselect"
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
        super().__init__(*args, **kwargs)
        self.state["adnext"] = \
            max(self.state.get("adnext", 0), int(self.params.get("admin", 0)))
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if self.is_done():
            self.logger.info("Polynomial selection already finished - nothing to do")
            return
        
        if "bestpoly" in self.state:
            self.bestpoly = Polynomial(self.state["bestpoly"].splitlines())
            self.bestpoly.setE(self.state["bestE"])
            self.logger.info("Best polynomial previously found in %s has "
                             "Murphy_E = %g", 
                             self.state["bestfile"], self.bestpoly.E)
        else:
            self.bestpoly = None
            self.logger.info("No polynomial was previously found")
        
        # Submit all the WUs we need to reach admax
        while self.need_more_wus():
            self.submit_one_wu()
        
        # Wait for all the WUs to finish
        while self.get_number_outstanding_wus() > 0:
            self.wait()
        
        if not self.bestpoly:
            self.logger.error ("No polynomial found. Consider increasing the "
                               "search range bound admax, or maxnorm")
            return
        self.logger.info("Finished, best polynomial from file %s has Murphy_E "
                         "= %g", self.state["bestfile"] , self.bestpoly.E)
        return
    
    def is_done(self):
        return "bestpoly" in self.state and not self.need_more_wus() and \
            self.get_number_outstanding_wus() == 0
    
    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return
        output_files = message.get_output_files()
        assert len(output_files) == 1
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
    
    def need_more_wus(self):
        return self.state["adnext"] < int(self.params["admax"])
    
    def submit_one_wu(self):
        adstart = self.state["adnext"]
        adend = adstart + int(self.params["adrange"])
        adend = min(adend, int(self.params["admax"]))
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


class FactorBaseOrFreerelTask(Task, metaclass=abc.ABCMeta):
    """ Common base class for programs that produce one output file from 
    the polynomial, i.e., factorbase and freerel 
    """
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Invariant: if we have a result (in self.state["outputfile"]) then we
        # must also have a polynomial (in self.state["poly"])
        if "outputfile" in self.state:
            assert "poly" in self.state
            # The target file must correspond to the polynomial "poly"
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        # Get best polynomial found by polyselect
        poly = self.send_request(Request.GET_POLYNOMIAL)
        if not poly:
            raise Exception("FactorBaseOrFreerelTask(): no polynomial "
                            "received from PolyselTask")
        
        # Check if we have already computed the target file for this polynomial
        if "poly" in self.state:
            prevpoly = Polynomial(self.state["poly"].splitlines())
            if poly != prevpoly:
                if "outputfile" in self.state:
                    self.logger.info("Received different polynomial, discarding old one")
                    del(self.state["outputfile"])
                self.state["poly"] = str(poly)
        else:
            self.state["poly"] = str(poly)
        
        if not "outputfile" in self.state:
            self.logger.info("Starting")
            # Write polynomial to a file
            polyfile = self.make_output_filename("poly")
            poly.create_file(polyfile, self.params)
            
            # Make file name for factor base/free relations file
            outputfile = self.make_output_filename(self.target)

            # Run command to generate factor base/free relations file
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
            self.logger.info("Finished")

        self.check_files_exist([self.state["outputfile"]], "output", 
                               shouldexist=True)
    
    def get_filename(self):
        if "outputfile" in self.state:
            return self.state["outputfile"]
        else:
            return None
    
    @abc.abstractmethod
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
    def paramnames(self):
        return super().paramnames + \
            ("alim", )
    target = "roots"
    def parse_stderr(self, stderr):
        pass


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
    def paramnames(self):
        return super().paramnames + \
            ("lpba", )
    target = "freerel"
    
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
        self.logger.info("Found %d free relations" % self.state["nfree"])
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
    def paramnames(self):
        return super().paramnames + \
            ("qmin", "qrange", "rels_wanted") + Polynomial.paramnames
    # We seek to this many bytes before the EOF to look for the "Total xxx reports" message
    file_end_offset = 1000
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        qmin = int(self.params.get("qmin", 0))
        if "qnext" in self.state:
            self.state["qnext"] = max(self.state["qnext"], qmin)
        else:
            self.state["qnext"] = int(self.params.get("qmin", self.params["alim"]))
            
        self.state.setdefault("rels_found", 0)
        self.state.setdefault("rels_wanted", 0)
        self.params.setdefault("max_wus", 10)
        self.state["rels_wanted"] = max(self.state.get("rels_wanted", 0), 
                                        int(self.params.get("rels_wanted", 0)))
        if self.state["rels_wanted"] == 0:
            # TODO: Choose sensible default value
            pass
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        # Get best polynomial found by polyselect
        poly = self.send_request(Request.GET_POLYNOMIAL)
        if not poly:
            raise Exception("SievingTask(): no polynomial received")
        # Write polynomial to a file
        polyfile = self.make_output_filename("poly")
        poly.create_file(polyfile, self.params)
        
        while self.state["rels_found"] < int(self.state["rels_wanted"]):
            kwargs = self.progparams[0].copy()
            q0 = self.state["qnext"]
            q1 = q0 + int(self.params["qrange"])
            outputfile = self.make_output_filename("%d-%d" % (q0, q1))
            self.check_files_exist([outputfile], "output", shouldexist=False)
            kwargs["q0"] = str(q0)
            kwargs["q1"] = str(q1)
            kwargs["poly"] = polyfile
            kwargs["factorbase"] = self.send_request(Request.GET_FACTORBASE_FILENAME)
            kwargs["out"] = outputfile
            p = self.programs[0](None, kwargs)
            self.submit_command(p, "%d-%d" % (q0, q1))
            self.state["qnext"] = q1
        self.logger.info("Reached target of %s relations, now have %d",
                         int(self.state["rels_wanted"]), self.state["rels_found"])
        self.logger.debug("Exit SievingTask.run(" + self.name + ")")
        return
    
    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return
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
        if filename is None:
            return self.state["rels_found"]
        else:
            return self.output_files[filename]
    
    def request_more_relations(self, additional):
        if additional > 0:
            self.state["rels_wanted"] += additional
        if self.state["rels_wanted"] > self.state["rels_found"]:
            self.send_notification(Notification.WANT_TO_RUN, None)
            self.logger.info("New goal for number of relations is %d, "
                             "need to sieve more",
                             self.state["rels_wanted"])
        else:
            self.logger.info("New goal for number of relations is %d, but "
                             "already have %d. No need to sieve more",
                             self.state["rels_wanted"], self.state["rels_found"])

class Duplicates1Task(Task, FilesCreator):
    """ Removes duplicate relations """
    @property
    def name(self):
        return "duplicates1"
    @property
    def title(self):
        return "Filtering - Duplicate Removal, splitting pass"
    @property
    def programs(self):
        return (cadoprograms.Duplicates1,)
    @property
    def paramnames(self):
        return super().paramnames + \
            ("nslices_log",)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nr_slices = 2**int(self.params["nslices_log"])
        self.already_split_input = self.make_db_dict(self.make_tablename("infiles"))
        self.slice_relcounts = self.make_db_dict(self.make_tablename("counts"))
        # Default slice counts to 0, in single DB commit
        self.slice_relcounts.setdefault(
            None, {str(i): 0 for i in range(0, self.nr_slices)})
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
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
        
        siever_files = self.send_request(Request.GET_SIEVER_FILENAMES)
        newfiles = [f for f in siever_files
                    if not f in self.already_split_input]
        self.logger.debug ("new files to split are: %s", newfiles)
        
        if not newfiles:
            self.logger.info("No new files to split")
        else:
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
                    current_counts = [self.send_request(Request.GET_SIEVER_RELNUMBER, f)]
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
        """ Takes lines of text and looks for slice counts as printed by dup1
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
    
    def get_nrels(self, idx):
        return self.slice_relcounts[str(idx)]
    
    def request_more_relations(self, additional):
        # TODO: Estimate how many more raw relations we need
        self.send_notification(Notification.WANT_MORE_RELATIONS, additional)
        self.send_notification(Notification.WANT_TO_RUN, None)

class Duplicates2Task(Task, FilesCreator):
    """ Removes duplicate relations """
    @property
    def name(self):
        return "duplicates2"
    @property
    def title(self):
        return "Filtering - Duplicate Removal, removal pass"
    @property
    def programs(self):
        return (cadoprograms.Duplicates2,)
    @property
    def paramnames(self):
        return super().paramnames + \
            ("nslices_log",)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nr_slices = 2**int(self.params["nslices_log"])
        self.already_done_input = self.make_db_dict(self.make_tablename("infiles"))
        self.slice_relcounts = self.make_db_dict(self.make_tablename("counts"))
        self.slice_relcounts.setdefault(
            None, {str(i): 0 for i in range(0, self.nr_slices)})
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        self.logger.info("Starting")
        for i in range(0, self.nr_slices):
            files = self.send_request(Request.GET_DUP1_FILENAMES, i.__eq__)
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
            rel_count = self.send_request(Request.GET_DUP1_RELCOUNT, i)
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
                self.already_done_input[f] = True
            outfilenames = {self.make_output_filename(f, i):i for f in files}
            self.add_output_files(outfilenames)
            self.logger.info("%d unique relations remain on slice %d", nr_rels, i)
            self.slice_relcounts[str(i)] = nr_rels
        self.logger.info("%d unique relations remain in total", self.get_nrels())
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

    def get_nrels(self):
        nrels = 0
        for i in range(0, self.nr_slices):
            nrels += self.slice_relcounts[str(i)]
        return nrels

    def request_more_relations(self, additional):
        # TODO: Estimate how many more raw relations we need
        self.send_notification(Notification.WANT_MORE_RELATIONS, additional)
        self.send_notification(Notification.WANT_TO_RUN, None)


class PurgeTask(Task):
    """ Removes singletons and computes excess """
    @property
    def name(self):
        return "singletons"
    @property
    def title(self):
        return "Filtering - Singleton removal"
    @property
    def programs(self):
        return (cadoprograms.Purge,)
    @property
    def paramnames(self):
        return super().paramnames + \
            ("keep", ) + Polynomial.paramnames
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        poly = self.send_request(Request.GET_POLYNOMIAL)
        polyfile = self.make_output_filename("poly")
        poly.create_file(polyfile, self.params)
        nfree = self.send_request(Request.GET_FREEREL_RELCOUNT)
        nunique = self.send_request(Request.GET_UNIQUE_RELCOUNT)
        if not nunique:
            raise Exception("No unique relation count received")
        nrels = nfree + nunique
        
        if "purgedfile" in self.state and nrels == self.state["input_nrels"]:
            self.logger.info("Already have a purged file, and no new input "
                             "relations available. Nothing to do")
            return
        
        self.state.pop("purgedfile", None)
        self.state.pop("input_nrels", None)
        
        self.logger.info("Reading %d unique and %d free relations, total %d"
                         % (nunique, nfree, nrels))
        purgedfile = self.make_output_filename("purged.gz")
        freerel_filename = self.send_request(Request.GET_FREEREL_FILENAME)
        unique_filenames = self.send_request(Request.GET_UNIQUE_FILENAMES)
        args = unique_filenames + [freerel_filename]
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
            # Update both atomically
            self.state.update({"purgedfile": purgedfile, "input_nrels": nrels})
            self.logger.info("Have enough relations")
            self.send_notification(Notification.HAVE_ENOUGH_RELATIONS, None)
        else:
            self.logger.info("Not enough relations, requesting more")
            self.request_more_relations(int(nunique * 0.1))
        self.logger.debug("Exit PurgeTask.run(" + self.name + ")")
    
    def request_more_relations(self, additional):
        # TODO: Estimate how many more unique relations we need
        self.send_notification(Notification.WANT_MORE_RELATIONS, additional)
        self.send_notification(Notification.WANT_TO_RUN, None)
    
    def get_purged_filename(self):
        return self.state["purgedfile"]
    
    def parse_stderr(self, stderr):
        # If stderr ends with 
        # b'excess < 0.10 * #primes. See -required_excess argument.'
        # then we need more relations from filtering and return False
        for line in stderr.decode("ascii").splitlines():
            if re.match("excess < \d+.\d+ \* #primes", line) or \
                    re.match("number of relations <= number of ideals", line):
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
        return "Filtering - Merging"
    @property
    def programs(self):
        return (cadoprograms.Merge, cadoprograms.Replay)
    @property
    def paramnames(self):
        return super().paramnames + \
            ("skip", "forbw", "coverNmax", "keep", "maxlevel", "ratio")
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if not "mergedfile" in self.state:
            self.logger.info("Starting")
            if "indexfile" in self.state:
                del(self.state["indexfile"])
            if "mergedfile" in self.state:
                del(self.state["mergedfile"])
            if "densefile" in self.state:
                del(self.state["densefile"])
            
            purged_filename = self.send_request(Request.GET_PURGED_FILENAME)
            historyfile = self.make_output_filename("history")
            args = ()
            kwargs = self.progparams[0].copy()
            kwargs["mat"] = purged_filename
            kwargs["out"] = historyfile
            p = self.programs[0](args, kwargs)
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            
            indexfile = self.make_output_filename("index")
            mergedfile = self.make_output_filename("small.bin")
            args = ()
            kwargs = self.progparams[1].copy()
            kwargs["binary"] = "1"
            kwargs["purged"] = purged_filename
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
    def paramnames(self):
        return super().paramnames + \
            ()
    # bwc.pl :complete seed=1 thr=1x1 mpi=1x1 matrix=c59.small.bin nullspace=left mm_impl=bucket interleaving=0 interval=100 mn=64 wdir=c59.bwc shuffled_product=1 bwc_bindir=/localdisk/kruppaal/build/cado-nfs/normal/linalg/bwc
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if not "dependency" in self.state:
            self.logger.info("Starting")
            workdir = self.make_output_dirname()
            self.make_directories(workdir)
            mergedfile = self.send_request(Request.GET_MERGED_FILENAME)
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
    def paramnames(self):
        return super().paramnames + \
            Polynomial.paramnames

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if not "kernel" in self.state:
            self.logger.info("Starting")
            polyfilename = self.make_output_filename("poly")
            poly = self.send_request(Request.GET_POLYNOMIAL)
            poly.create_file(polyfilename, self.params)
            kernelfilename = self.make_output_filename("kernel")
            
            purgedfilename = self.send_request(Request.GET_PURGED_FILENAME)
            indexfilename = self.send_request(Request.GET_INDEX_FILENAME)
            densefilename = self.send_request(Request.GET_DENSE_FILENAME)
            dependencyfilename = self.send_request(Request.GET_DEPENDENCY_FILENAME)
            
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
    def paramnames(self):
        return super().paramnames + \
            ("N",)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.factors = self.make_db_dict(self.make_tablename("factors"))
        self.add_factor(int(self.params["N"]))
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        if not self.is_done():
            self.logger.info("Starting")
            polyfilename = self.make_output_filename("poly")
            poly = self.send_request(Request.GET_POLYNOMIAL)
            poly.create_file(polyfilename, self.params)
            
            purgedfilename = self.send_request(Request.GET_PURGED_FILENAME)
            indexfilename = self.send_request(Request.GET_INDEX_FILENAME)
            kernelfilename = self.send_request(Request.GET_KERNEL_FILENAME)
            prefix = self.send_request(Request.GET_LINALG_PREFIX)
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
                (identifier, rc, stdout, stderr, output_files) = \
                    self.submit_command(p, "dep%d" % self.state["next_dep"])
                if not stdout.decode("ascii").strip() == "Failed":
                    factorlist = list(map(int,stdout.decode("ascii").split()))
                    # FIXME: Can sqrt print more/less than 2 factors?
                    assert len(factorlist) == 2
                    self.add_factor(factorlist[0])
                self.state["next_dep"] += 1
            self.logger.info("finished")
        self.logger.info("Factors: %s" % " ".join(self.get_factors()))
        self.logger.debug("Exit SqrtTask.run(" + self.name + ")")
    
    def is_done(self):
        for (factor, isprime) in self.factors.items():
            if not isprime:
                return False
        return True
    
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

class StartClientsTask(Task):
    """ Starts clients on slave machines """
    @property
    def name(self):
        return "slaves"
    @property
    def title(self):
        return "Client Launcher"
    @property
    def programs(self):
        return (cadoprograms.WuClient,)
    @property
    def paramnames(self):
        return super().paramnames + \
            ('hostnames', 'scriptpath', "nrclients")
    
    def __init__(self, address, port, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.used_ids = {}
        self.pids = self.make_db_dict(self.make_tablename("client_pids"))
        self.hosts = self.make_db_dict(self.make_tablename("client_hosts"))
        
        if 'scriptpath' in self.params:
            self.progparams[0]['execpath'] = self.params['scriptpath']
        self.progparams[0].setdefault('server', "http://%s:%d" % (address, port))
        
        # If hostnames are of the form @file, read host names from file,
        # one host name per line
        match = re.match("@(.*)", self.params["hostnames"])
        if match:
            with open(match.group(1)) as f:
                self.hosts_to_launch = [line for line in f]
        else:
            self.hosts_to_launch = self.params["hostnames"].split(",")
    
    def is_alive(self, clientid):
        # Simplistic: just test if process with that pid exists and accepts
        # signals from us. TODO: better testing here, probably with ps|grep
        # or some such
        (rc, stdout, stderr) = self.kill_client(clientid, signal=0)
        return (rc == 0)
    
    def launch_clients(self):
        for host in self.hosts_to_launch:
            self.launch_host(host.strip())
        s = ", ".join(["%s (Host %s, PID %d)" % (cid, self.hosts[cid], self.pids[cid]) for cid in self.pids])
        self.logger.info("Launched clients: %s" % s)
    
    def make_unique_id(self, host):
        # Make a unique client id for host
        clientid = host
        i = 1
        while clientid in self.used_ids:
            assert clientid in self.pids
            assert clientid in self.hosts
            i += 1
            clientid = "%s%d" % (host, i)
        self.used_ids[clientid] = True
        return clientid
    
    # Cases:
    # Client was never started. Start it, add to state
    # Client was started, but does not exist any more. Remove from state, then start and add again
    # Client was started, and does still exists. Nothing to do.
    
    def launch_host(self, host):
        clientid = self.make_unique_id(host)
        # Check if client is already running
        if clientid in self.pids:
            assert self.hosts[clientid] == host
            if self.is_alive(clientid):
                self.logger.info("Client %s on host %s with PID %d already running",
                                 clientid, host, self.pids[clientid])
                return
            else:
                del(self.pids[clientid])
                del(self.hosts[clientid])
        
        self.logger.info("Starting client id %s on host %s", clientid, host)
        self.progparams[0]['clientid'] = clientid
        wuclient = cadoprograms.WuClient([], self.progparams[0], stdout = "/dev/null", stderr = "/dev/null", bg=True)
        process = cadocommand.RemoteCommand(wuclient, host, self.parameters, self.path_prefix)
        (rc, stdout, stderr) = process.wait()
        if rc != 0:
            self.logger.warning("Starting client on %(host)s failed. "
                                "Consult log file for details." % host)
            return
        self.pids[clientid] = int(stdout)
        self.hosts[clientid] = host

    def kill_all_clients(self):
        # Need the keys() to make a copy as dict will change in loop body
        for clientid in list(self.pids.keys()):
            (rc, stdout, stderr) = self.kill_client(clientid)
            if rc == 0:
                self.logger.info("Stopped client %s (Host %s, PID %d)",
                                 clientid, self.hosts[clientid], self.pids[clientid])
                del(self.pids[clientid])
                del(self.hosts[clientid])
            else:
                self.logger.warning("Stopping client %s (Host %s, PID %d) failed",
                                    clientid, self.hosts[clientid], self.pids[clientid])
    
    def kill_client(self, clientid, signal = None):
        pid = self.pids[clientid]
        host = self.hosts[clientid]
        params = {}
        if not signal is None:
            params["signal"] = str(signal)
        kill = cadoprograms.Kill((pid,), params)
        process = cadocommand.RemoteCommand(kill, host, self.parameters, self.path_prefix)
        return process.wait()

class Message(object, metaclass = abc.ABCMeta):
    def __init__(self, sender, key, value):
        self.sender = sender
        self.key = key
        self.value = value
    def get_sender(self):
        return self.sender
    def get_key(self):
        return self.key
    def get_value(self):
        return self.value

class Notification(Message):
    FINISHED_POLYNOMIAL_SELECTION = 0
    WANT_MORE_RELATIONS = 1
    HAVE_ENOUGH_RELATIONS = 2
    REGISTER_FILENAME = 3
    SUBMIT_WU = 4
    VERIFY_WU = 5
    WANT_TO_RUN = 6

class Request(Message):
    GET_POLYNOMIAL = 0
    GET_FACTORBASE_FILENAME = 1
    GET_FREEREL_FILENAME = 2
    GET_FREEREL_RELCOUNT = 3
    GET_SIEVER_FILENAMES = 4
    GET_SIEVER_RELCOUNT = 5
    GET_DUP1_FILENAMES = 6
    GET_DUP1_RELCOUNT = 7
    GET_UNIQUE_RELCOUNT = 8
    GET_UNIQUE_FILENAMES = 9
    GET_PURGED_FILENAME = 10
    GET_MERGED_FILENAME = 11
    GET_INDEX_FILENAME = 12
    GET_DENSE_FILENAME = 13
    GET_DEPENDENCY_FILENAME = 14
    GET_LINALG_PREFIX = 15
    GET_KERNEL_FILENAME = 16
    def __init__(self, sender, key, value = None):
        super().__init__(sender, key, value)


class CompleteFactorization(wudb.DbAccess, cadoparams.UseParameters, patterns.Mediator):
    """ The complete factorization, aggregate of the individual tasks """
    """ Runs the square root """
    @property
    def name(self):
        return "tasks"
    
    CAN_CANCEL_WUS = 0
    
    def __init__ (self, path_prefix, *args, **kwargs):
        super().__init__(*args, path_prefix = path_prefix, **kwargs)
        self.logger = logging.getLogger("Complete Factorization")
        self.params = self.myparams(("name", "workdir"))
        self.db_listener = self.make_db_listener()
        self.registered_filenames = self.make_db_dict(self.params["name"] + '_server_registered_filenames')
        self.chores = []
        
        # Init WU DB
        self.wuar = self.make_wu_access()
        self.wuar.create_tables()
        
        # Set up WU server
        serveraddress = "localhost"
        serverport = 8001
        uploaddir = self.params["workdir"].rstrip(os.sep) + os.sep + self.params["name"] + ".upload/"
        self.server = wuserver.ServerLauncher(serveraddress, serverport, False, self.get_db_filename(), 
            self.registered_filenames, uploaddir, bg = True)
        
        # Init client lists
        self.clients = []
        for (path, key) in self.get_parameters().find(['slaves'], 'hostnames'):
            self.clients.append(StartClientsTask(serveraddress, serverport, *args, 
                                                 mediator = self,
                                                 path_prefix = path, **kwargs))
        
        parampath = self.get_param_path()
        sievepath = parampath + ['sieve']
        filterpath = parampath + ['filter']
        linalgpath = parampath + ['linalg']
        
        self.polysel = PolyselTask(*args,
                                   mediator = self,
                                   path_prefix = parampath,
                                   db_listener = self.db_listener,
                                   **kwargs)
        self.fb = FactorBaseTask(*args,
                                 mediator = self,
                                 path_prefix = sievepath, 
                                 **kwargs)
        self.freerel = FreeRelTask(*args,
                                   mediator = self,
                                   path_prefix = sievepath,
                                   **kwargs)
        self.sieving = SievingTask(*args,
                                   mediator = self,
                                   path_prefix = sievepath,
                                   db_listener = self.db_listener,
                                   **kwargs)
        self.dup1 = Duplicates1Task(*args,
                                    mediator = self,
                                    path_prefix = filterpath,
                                    **kwargs)
        self.dup2 = Duplicates2Task(*args,
                                    mediator = self,
                                    path_prefix = filterpath,
                                    **kwargs)
        self.purge = PurgeTask(*args,
                               mediator = self,
                               path_prefix = filterpath,
                               **kwargs)
        self.merge = MergeTask(*args,
                               mediator = self,
                               path_prefix = filterpath, 
                               **kwargs)
        self.linalg = LinAlgTask(*args,
                                 mediator = self,
                                 path_prefix = linalgpath, 
                                 **kwargs)
        self.characters = CharactersTask(*args,
                                         mediator = self,
                                         path_prefix = linalgpath,
                                         **kwargs)
        self.sqrt = SqrtTask(*args,
                             mediator = self,
                             path_prefix = parampath,
                             **kwargs)
        
        # Defines an order on tasks in which tasks that want to run should be run
        self.tasks = (self.polysel, self.fb, self.freerel, self.sieving,
                      self.dup1, self.dup2, self.purge, self.merge,
                      self.linalg, self.characters, self.sqrt)
        
        # Assume that all tasks want to run. Ff they are finished already, 
        # they will just return immediately
        self.tasks_that_want_to_run = list(self.tasks)
        
        self.request_map = {
            Request.GET_POLYNOMIAL: self.polysel.get_poly,
            Request.GET_FACTORBASE_FILENAME: self.fb.get_filename,
            Request.GET_FREEREL_FILENAME: self.freerel.get_filename,
            Request.GET_FREEREL_RELCOUNT: self.freerel.get_nrels,
            Request.GET_SIEVER_FILENAMES: self.sieving.get_output_filenames,
            Request.GET_SIEVER_RELCOUNT: self.sieving.get_nrels,
            Request.GET_DUP1_FILENAMES: self.dup1.get_output_filenames,
            Request.GET_DUP1_RELCOUNT: self.dup1.get_nrels,
            Request.GET_UNIQUE_RELCOUNT: self.dup2.get_nrels,
            Request.GET_UNIQUE_FILENAMES: self.dup2.get_output_filenames,
            Request.GET_PURGED_FILENAME: self.purge.get_purged_filename,
            Request.GET_MERGED_FILENAME: self.merge.get_merged_filename,
            Request.GET_INDEX_FILENAME: self.merge.get_index_filename,
            Request.GET_DENSE_FILENAME: self.merge.get_dense_filename,
            Request.GET_DEPENDENCY_FILENAME: self.linalg.get_dependency_filename,
            Request.GET_LINALG_PREFIX: self.linalg.get_prefix,
            Request.GET_KERNEL_FILENAME: self.characters.get_kernel_filename
        }
    
    def run(self):
        self.server.serve()
        
        for clients in self.clients:
            clients.launch_clients()
        
        while self.run_next_task():
            self.do_chores()
        
        self.do_chores()
        
        for c in self.clients:
            c.kill_all_clients()
        
        self.server.shutdown()
    
    def run_next_task(self):
        for task in self.tasks:
            if task in self.tasks_that_want_to_run:
                self.logger.info("Next task that wants to run: %s", task.title)
                self.tasks_that_want_to_run.remove(task)
                task.run()
                return True
        return False
    
    def do_chores(self):
        if self.chores:
            chore = self.chores.pop()
            if chore == self.CAN_CANCEL_WUS:
                self.logger.info("Cancelling remaining workunits")
                self.cancel_available_wus()
            else:
                raise Exception("Unknown chore %s" % chore)
    
    def add_wu(self, wutext):
        # print ("WU:\n%s" % wutext)
        self.wuar.create(wutext)
    
    def verify_wu(self, wuinfo):
        for wuid in wuinfo:
            self.wuar.verification(wuid, wuinfo[wuid])
    
    def cancel_available_wus(self):
        self.wuar.cancel_all_available()
    
    def register_filename(self, d):
        for key in d:
            if not key in self.registered_filenames:
                self.logger.debug("Registering file name %s with target %s",
                                  key, d[key])
                self.registered_filenames[key] = d[key]
            elif d[key] != self.registered_filenames[key]:
                # It was already defined with a different target, error
                raise Exception("Filename %s, to be registered for target %s, "
                                "already registered for target %s" %
                                (key, d[key], self.registered_filenames[key]))
            else:
                # Was already registered with the same target. Nothing to do
                pass
    
    def relay_notification(self, message):
        assert isinstance(message, Notification)
        sender = message.get_sender()
        key = message.get_key()
        value = message.get_value()
        self.logger.debug("Received notification from %s, key = %s, values = %s",
                           sender, key, value)
        """ The relay for letting Tasks talk to us and each other """
        if key == Notification.WANT_MORE_RELATIONS:
            if sender is self.purge:
                self.dup2.request_more_relations(value)
            elif sender is self.dup2:
                self.dup1.request_more_relations(value)
            elif sender is self.dup1:
                self.sieving.request_more_relations(value)
            else:
                raise Exception("Got WANT_MORE_RELATIONS from unknown sender")
        elif key == Notification.HAVE_ENOUGH_RELATIONS:
            if sender is self.purge:
                # TODO: cancel only sieving WUs?
                self.chores.append(self.CAN_CANCEL_WUS)
            else:
                raise Exception("Got HAVE_ENOUGH_RELATIONS from unknown sender")
        elif key == Notification.REGISTER_FILENAME:
            if isinstance(sender, ClientServerTask):
                self.register_filename(value)
            else:
              raise Exception("Got REGISTER_FILENAME, but not from a ClientServerTask")
        elif key == Notification.SUBMIT_WU:
            if isinstance(sender, ClientServerTask):
                self.add_wu(value)
            else:
              raise Exception("Got SUBMIT_WU, but not from a ClientServerTask")
        elif key == Notification.VERIFY_WU:
            if isinstance(sender, ClientServerTask):
                self.verify_wu(value)
            else:
              raise Exception("Got VERIFY_WU, but not from a ClientServerTask")
        elif key == Notification.WANT_TO_RUN:
            if sender in self.tasks_that_want_to_run:
                raise Exception("Got request from %s to run, but it was in run queue already",
                                sender)
            else:
                self.tasks_that_want_to_run.append(sender)
        else:
            raise KeyError("Notification from %s has unknown key %s" % (sender, key))
    
    def answer_request(self, request):
        assert isinstance(request, Request)
        sender = request.get_sender()
        key = request.get_key()
        value = request.get_value()
        self.logger.debug("Received request from %s, key = %s, values = %s",
                           sender, key, value)
        if key in self.request_map:
            if value is None:
                return self.request_map[key]()
            else:
                return self.request_map[key](value)
        else:
            raise KeyError("Unknown Request key %s from sender %s" %
                           (str(key), str(sender)))
    
    def handle_message(self, message):
        if isinstance(message, Notification):
            self.relay_notification(Notification)
        elif isinstance(message, Request):
            return self.answer_request(message)
        else:
            raise TypeError("Message is neither Notification nor Request")
