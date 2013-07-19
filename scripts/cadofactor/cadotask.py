#!/usr/bin/env python3
import re
import os.path
from fractions import gcd
import abc
import random
import time
import wudb
import logging
import socket
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
            return [f for (f, s) in self.output_files.items() if condition(s)]
    
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
            if re.match(r"No polynomial found", line):
                return
            # If this is a comment line telling the Murphy E value, 
            # extract the value and store it
            match = re.match(r"\s*#\s*MurphyE\s*\(.*\)=(.*)$", line)
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
        with open(str(filename), "w") as f:
            f.write(str(self))
            for key in self.paramnames:
                if key in params:
                    f.write(key + ": %s\n" % params[key])


class FilePath(object):
    """ A class that represents a path to a file, where the path should be
    somewhat relocateable.
    
    In particular, we separate the path to the working directory, and the file
    path relative to the working directory. For persistent storage in the DB,
    the path relative to the workdir should be used, whereas for any file
    accesses, the full path needs to be used.
    """
    def __init__(self, workdir, filepath):
        self.workdir = workdir.rstrip(os.sep)
        self.filepath = filepath
    def __str__(self):
        return "%s%s%s" % (self.workdir, os.sep, self.filepath)
    def get_relative(self):
        return self.filepath


class WorkDir(object):
    """ A class that allows generating file and directory names under a 
    working directory.
    
    The directory layout is as follows:
    The current project (i.e., the factorization) has a jobname, e.g.,
    "RSA512". Each task has a name, e.g., "sieving".
    A task can create various files under 
    workdir/jobname.taskname.file
    or put them in a subdirectory
    workdir/jobname.taskname/file
    or, for multiple subdirectories,
    workdir/jobname.taskname/subdir/file
    
    >>> f = WorkDir("/foo/bar", "jobname", "taskname")
    >>> str(f.make_dirname())
    '/foo/bar/jobname.taskname/'
    >>> str(f.make_dirname('subdir'))
    '/foo/bar/jobname.taskname/subdir/'
    >>> str(f.make_filename('file'))
    '/foo/bar/jobname.taskname.file'
    >>> str(f.make_filename('file', use_subdir=True))
    '/foo/bar/jobname.taskname/file'
    >>> str(f.make_filename('file', use_subdir=True, subdir='subdir'))
    '/foo/bar/jobname.taskname/subdir/file'
    """
    def __init__(self, workdir, jobname, taskname):
        self.workdir = workdir.rstrip(os.sep)
        self.jobname = jobname
        self.taskname = taskname
    
    def path_in_workdir(self, filename):
        return FilePath(self.workdir, filename)
    
    def _make_path(self, extra):
        """ Make a path of the form: "workdir/jobname.taskname""extra" """
        return self.path_in_workdir("%s.%s%s" % (self.jobname, self.taskname,
                                                 extra))
    
    def make_dirname(self, subdir = None):
        """ Make a directory name of the form workdir/jobname.taskname/ if 
        subdir is not given, or workdir/jobname.taskname/subdir/ if it is
        """
        if subdir:
            return self._make_path("%s%s%s" % (os.sep, subdir, os.sep))
        else:
            return self._make_path(os.sep)
    
    def make_filename(self, name, use_subdir = False, subdir = None):
        """ If use_subdir is False, make a filename of the form 
        workdir/jobname.taskname.name 
        If use_subdir is True and subdir is None, make a filename of the form 
        workdir/jobname.taskname/name
        If use_subdir is True and subdir is a string, make a filename of the form 
        workdir/jobname.taskname/subdir/name
        """
        if use_subdir:
            if subdir:
                return self._make_path("%s%s%s%s" % (os.sep, subdir, os.sep,
                                                     name))
            else:
                return self._make_path("%s%s" % (os.sep, name))
        else:
            assert subdir is None
            return self._make_path(".%s" % name)
    
    def make_directories(self, subdirs = None):
        dirname = str(self.make_dirname())
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        if subdirs:
            for subdir in subdirs:
                dirname = str(self.make_dirname(subdir))
                if not os.path.isdir(dirname):
                    os.mkdir(dirname)
        return


class Task(patterns.Colleague, wudb.DbAccess, cadoparams.UseParameters,
           metaclass=abc.ABCMeta):
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
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        ''' Sets up a database connection and a DB-backed dictionary for 
        parameters. Reads parameters from DB, and merges with hierarchical
        parameters in the parameters argument. Parameters passed in by 
        parameters argument do not override values in the DB-backed 
        parameter dictionary.
        '''
        
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.logger = logging.getLogger(self.title)
        self.logger.debug("Enter Task.__init__(%s)", 
                          self.name)
        # DB-backed dictionary with the state of this task
        self.state = self.make_db_dict(self.make_tablename())
        self.logger.debug("state = %s", self.state)
        # Set default parametes for this task, if any are given
        self.params = self.myparams(self.paramnames)
        self.logger.debug("param_prefix = %s, get_param_path = %s", 
                          self.get_param_prefix(), self.get_param_path())
        self.logger.debug("params = %s", self.params)
        # Set default parameters for our programs
        self.progparams = []
        for prog in self.programs:
            progparams = self.myparams(prog.get_accepted_keys(), prog.name)
            self.progparams.append(progparams)
        # FIXME: whether to init workdir or not should not be controlled via
        # presence of a "workdir" parameter, but by class definition
        if "workdir" in self.params:
            self.workdir = WorkDir(self.params["workdir"], self.params["name"], 
                               self.name)
        self.logger.debug("Exit Task.__init__(%s)", self.name)
        return
    
    @staticmethod
    def check_tablename(name):
        no_ = name.replace("_", "")
        if not no_[0].isalpha() or not no_[1:].isalnum():
            raise Exception("%s is not valid for an SQL table name" % name)
    
    def make_tablename(self, extra = None):
        """ Return a table name for the DB-backed dictionary """
        # Maybe replace SQL-disallowed characters here, like digits and '.' ? 
        # Could be tricky to avoid collisions
        name = self.name
        if extra:
            name = name + '_' + extra
        self.check_tablename(name)
        return name

    def translate_input_filename(self, filename):
        return filename

    def test_outputfile_exists(self, filename):
        return os.path.isfile(str(filename))
    
    @staticmethod
    def check_files_exist(filenames, filedesc, shouldexist):
        """ Check that the output files in "filenames" exist or don't exist, 
        according to shouldexist.
        
        Raise IOError if any check fails, return None
        """
        for f in filenames:
            exists = os.path.isfile(str(f))
            if shouldexist and not exists:
                raise IOError("%s file %s does not exist" % (filedesc, f))
            elif not shouldexist and exists:
                raise IOError("%s file %s already exists" % (filedesc, f))
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
            # pylint: disable=E1101
            self.updateObserver(message)
        return result
    
    def filter_notification(self, message):
        wuid = message.get_wu_id()
        rc = message.get_exitcode(0)
        stdout = message.get_stdout(0)
        stderr = message.get_stderr(0)
        output_files = message.get_output_files()
        self.logger.message("%s: Received notification for wuid=%s, rc=%d, output_files=[%s]",
                            self.name, wuid, rc, ", ".join(output_files))
        if rc != 0:
            self.logger.error("Return code is: %d", rc)
        if stdout:
            self.logger.debug("stdout is: %s", stdout)
        if stderr:
            self.logger.debug("stderr is: %s", stderr)
        if output_files:
            self.logger.message("Output files are: %s", ", ".join(output_files))
        (name, task, identifier) = self.split_wuname(wuid)
        if name != self.params["name"] or task != self.name:
            # This notification is not for me
            self.logger.debug("%s: Notification is not for me", self.name)
            return
        self.logger.message("%s: Notification is for me", self.name)
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
    
    def get_number_outstanding_wus(self):
        return 0
    
    def verification(self, message, ok):
        pass

    def get_state_filename(self, key):
        if not key in self.state:
            return None
        return self.workdir.path_in_workdir(self.state[key])


class ClientServerTask(Task, patterns.Observer):
    @abc.abstractproperty
    def paramnames(self):
        return super().paramnames + ("maxwu",)
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.state.setdefault("wu_submitted", 0)
        self.state.setdefault("wu_received", 0)
        self.params.setdefault("maxwu", "10")
        assert self.get_number_outstanding_wus() >= 0
        self.send_notification(Notification.SUBSCRIBE_WU_NOTIFICATIONS, None)
    
    def submit_command(self, command, identifier):
        ''' Submit a workunit to the database. '''
        
        while self.get_number_available_wus() >= int(self.params["maxwu"]):
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

    def get_number_available_wus(self):
        return self.send_request(Request.GET_NR_AVAILABLE_WU, None)

    def test_outputfile_exists(self, filename):
        # Can't test
        return False
    
    def wait(self):
        # Ask the mediator to check for workunits of status Received, 
        # and if there are any, to send WU result notifications to the 
        # subscribed listeners.
        # If we get notification on new results reliably from the HTTP server,
        # we might not need this poll. But they probably won't be totally
        # reliable
        if not self.send_request(Request.GET_WU_RESULT):
            time.sleep(1)


class PolyselTask(ClientServerTask, patterns.Observer):
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
            ("adrange", "admin", "admax") + \
            Polynomial.paramnames
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.state["adnext"] = \
            max(self.state.get("adnext", 0), int(self.params.get("admin", 0)))
        self.bestpoly = None
        if "bestpoly" in self.state:
            self.bestpoly = Polynomial(self.state["bestpoly"].splitlines())
            self.bestpoly.setE(self.state["bestE"])
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.name, self.state)
        
        if self.is_done():
            self.logger.info("Polynomial selection already finished - nothing to do")
            return
        
        if not self.bestpoly is None:
            self.logger.info("Best polynomial previously found in %s has "
                             "Murphy_E = %g", 
                             self.state["bestfile"], self.bestpoly.E)
        else:
            self.logger.info("No polynomial was previously found")
        
        # Submit all the WUs we need to reach admax
        while self.need_more_wus():
            self.submit_one_wu()
        
        # Wait for all the WUs to finish
        while self.get_number_outstanding_wus() > 0:
            self.wait()
        
        if self.bestpoly is None:
            self.logger.error ("No polynomial found. Consider increasing the "
                               "search range bound admax, or maxnorm")
            return
        self.logger.info("Finished, best polynomial from file %s has Murphy_E "
                         "= %g", self.state["bestfile"] , self.bestpoly.E)
        self.write_poly_file()
        return
    
    def is_done(self):
        return not self.bestpoly is None and not self.need_more_wus() and \
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
        if self.bestpoly is None or poly.E > self.bestpoly.E:
            self.bestpoly = poly
            update = {"bestE": poly.E, "bestpoly": str(poly),
                      "bestfile": outputfile}
            self.state.update(update)
            self.logger.info("New best polynomial from file %s:"
                             " Murphy E = %g" % (outputfile, poly.E))
            self.logger.debug("New best polynomial is:\n%s", poly)
        else:
            self.logger.info("Best polynomial from file %s with E=%g is "
                             "no better than current best with E=%g",
                             outputfile, poly.E, self.bestpoly.E)
        return True
    
    def write_poly_file(self):
        filename = self.workdir.make_filename("poly")
        self.bestpoly.create_file(filename, self.params)
        self.state["polyfilename"] = filename.get_relative()
    
    def get_poly(self):
        if not "bestpoly" in self.state:
            return None
        return Polynomial(self.state["bestpoly"].splitlines())
    
    def get_poly_filename(self):
        return self.get_state_filename("polyfilename")

    def need_more_wus(self):
        return self.state["adnext"] < int(self.params["admax"])
    
    def submit_one_wu(self):
        adstart = self.state["adnext"]
        adend = adstart + int(self.params["adrange"])
        adend = min(adend, int(self.params["admax"]))
        outputfile = self.workdir.make_filename("%d-%d" % (adstart, adend))
        if self.test_outputfile_exists(outputfile):
            self.logger.info("%s already exists, won't generate again",
                             outputfile)
        else:
            # Remove admin and admax from the parameter-file-supplied
            # parameters as those would conflict with the computed values
            self.progparams[0].pop("admin", None)
            self.progparams[0].pop("admax", None)
            p = cadoprograms.Polyselect2l(admin=adstart, admax=adend,
                                          stdout=str(outputfile),
                                          **self.progparams[0])
            self.submit_command(p, "%d-%d" % (adstart, adend))
        self.state["adnext"] = adend


class FactorBaseTask(Task):
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


    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        # Invariant: if we have a result (in self.state["outputfile"]) then we
        # must also have a polynomial (in self.state["poly"])
        if "outputfile" in self.state:
            assert "poly" in self.state
            # The target file must correspond to the polynomial "poly"
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.name, self.state)
        
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
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            
            # Make file name for factor base/free relations file
            outputfilename = self.workdir.make_filename("roots")

            # Run command to generate factor base/free relations file
            p = cadoprograms.MakeFB(poly=polyfilename,
                                    stdout = str(outputfilename),
                                    **self.progparams[0])
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            if rc:
                raise Exception("Program failed")
            
            self.state["outputfile"] = outputfilename.get_relative()
            self.logger.info("Finished")

        self.check_files_exist([self.get_filename()], "output", 
                               shouldexist=True)
    
    def get_filename(self):
        return self.get_state_filename("outputfile")


class FreeRelTask(Task):
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
            ("lpba", "lpbr")
    wanted_regex = {
        'nfree': (r'# Free relations: (\d+)', int),
        'nprimes': (r'Renumbering struct: nprimes=(\d+)', int),
        'minindex': (r'Renumbering struct: min_index=(\d+)', int)
    }
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        # Invariant: if we have a result (in self.state["freerelfilename"])
        # then we must also have a polynomial (in self.state["poly"])
        if "freerelfilename" in self.state:
            assert "poly" in self.state
            # The target file must correspond to the polynomial "poly"
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.name, self.state)
        
        # Get best polynomial found by polyselect
        poly = self.send_request(Request.GET_POLYNOMIAL)
        if not poly:
            raise Exception("FactorBaseOrFreerelTask(): no polynomial "
                            "received from PolyselTask")
        
        # Check if we have already computed the target file for this polynomial
        if "poly" in self.state:
            prevpoly = Polynomial(self.state["poly"].splitlines())
            if poly != prevpoly:
                if "freerelfilename" in self.state:
                    self.logger.info("Received different polynomial, discarding old one")
                    del(self.state["freerelfilename"])
                    del(self.state["renumberfilename"])
                self.state["poly"] = str(poly)
        else:
            self.state["poly"] = str(poly)
        
        if not "freerelfilename" in self.state:
            self.logger.info("Starting")
            # Write polynomial to a file
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            
            # Make file name for factor base/free relations file
            freerelfilename = self.workdir.make_filename("freerel")
            renumberfilename = self.workdir.make_filename("renumber")

            # Run command to generate factor base/free relations file
            p = cadoprograms.FreeRel(poly=polyfilename,
                                     renumber=renumberfilename,
                                     stdout=str(freerelfilename),
                                     **self.progparams[0])
            (identifier, rc, stdout, stderr, output_files) = \
                    self.submit_command(p, "")
            if rc:
                raise Exception("Program failed")
            found = self.parse_file(stderr)
            self.state.update(found)
            self.logger.info("Found %d free relations" % self.state["nfree"])
            
            self.state["freerelfilename"] = freerelfilename.get_relative()
            self.state["renumberfilename"] = renumberfilename.get_relative()
            self.logger.info("Finished")

        self.check_files_exist([self.get_freerel_filename(),
                                self.get_renumber_filename()], "output", 
                               shouldexist=True)

    def parse_file(self, stderr):
        found = {}
        for line in stderr.decode("ascii").splitlines():
            for (key, (regex, datatype)) in self.wanted_regex.items():
                match = re.match(regex, line)
                if match:
                    if key in found:
                        raise Exception("Received two values for %s" % key)
                    found[key] = datatype(match.group(1))
        
        for key in self.wanted_regex:
            if not key in found:
                raise Exception("Received no value for %s" % key)
        return found
    
    def get_freerel_filename(self):
        return self.get_state_filename("freerelfilename")
    
    def get_renumber_filename(self):
        return self.get_state_filename("renumberfilename")
    
    def get_nrels(self):
        return self.state["nfree"]
    
    def get_nprimes(self):
        return self.state["nprimes"]

    def get_minindex(self):
        return self.state["minindex"]

class SievingTask(ClientServerTask, FilesCreator, patterns.Observer):
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
            ("qmin", "qrange", "rels_wanted", "alim")
    # We seek to this many bytes before the EOF to look for the "Total xxx reports" message
    file_end_offset = 1000
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        qmin = int(self.params.get("qmin", 0))
        if "qnext" in self.state:
            self.state["qnext"] = max(self.state["qnext"], qmin)
        else:
            self.state["qnext"] = int(self.params.get("qmin", self.params["alim"]))
        
        self.state.setdefault("rels_found", 0)
        self.state.setdefault("rels_wanted", 0)
        self.params.setdefault("maxwu", "10")
        self.state["rels_wanted"] = max(self.state.get("rels_wanted", 0), 
                                        int(self.params.get("rels_wanted", 0)))
        if self.state["rels_wanted"] == 0:
            # TODO: Choose sensible default value
            pass
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        while self.state["rels_found"] < self.state["rels_wanted"]:
            q0 = self.state["qnext"]
            q1 = q0 + int(self.params["qrange"])
            outputfilename = self.workdir.make_filename("%d-%d" % (q0, q1))
            self.check_files_exist([outputfilename], "output", shouldexist=False)
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            factorbase = self.send_request(Request.GET_FACTORBASE_FILENAME)
            p = cadoprograms.Las(q0=q0, q1=q1,
                                 poly=polyfilename, factorbase=factorbase,
                                 out=outputfilename,
                                 **self.progparams[0])
            self.submit_command(p, "%d-%d" % (q0, q1))
            self.state["qnext"] = q1
        self.logger.info("Reached target of %d relations, now have %d",
                         self.state["rels_wanted"], self.state["rels_found"])
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
                match = re.match(r"# Total (\d+) reports ", line)
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
            # Fixme: don't access self.output_files directly
            return self.output_files[filename]
    
    def request_more_relations(self, target):
        if target > self.state["rels_wanted"]:
            self.state["rels_wanted"] = target
        if self.state["rels_wanted"] > self.state["rels_found"]:
            self.send_notification(Notification.WANT_TO_RUN, None)
            self.logger.info("New goal for number of relations is %d, "
                             "currently have %d. Need to sieve more",
                             self.state["rels_wanted"], self.state["rels_found"])
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
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
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
            if parts != self.nr_slices:
                # TODO: ask interactively (or by -recover) whether to delete 
                # old output files and generate new ones, if input file is 
                # still available
                # If input file is not available but the previously split
                # parts are, we could join them again... not sure if want
                raise Exception("%s was previously split into %d parts, "
                                "now %d parts requested",
                                infile, parts, self.nr_slices)
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
            self.logger.info("Splitting %d new files", len(newfiles))
            self.workdir.make_directories(map(str, range(0, self.nr_slices)))
            # TODO: can we recover from missing input files? Ask Sieving to
            # generate them again? Just ignore the missing ones?
            self.check_files_exist(newfiles, "input", shouldexist = True)
            # Split the new files
            if self.nr_slices == 1:
                # If we should split into only 1 part, we don't actually
                # split at all. We simply write the input file name
                # to the table of output files, so the next stages will 
                # read the original siever output file, thus avoiding 
                # having another copy of the data on disk. Since we don't
                # process the file at all, we need to ask the Siever task
                # for the relation count in this file
                # TODO: pass a list or generator expression in the request
                # here?
                total = 0
                for f in newfiles:
                    count = self.send_request(Request.GET_SIEVER_RELCOUNT, f)
                    total += count
                self.slice_relcounts["0"] += total
                update = dict(zip(newfiles, [self.nr_slices] * len(newfiles)))
                self.already_split_input.update(update)
                update = dict(zip(newfiles, [0] * len(newfiles)))
                self.add_output_files(update)
            else:
                outfilenames = {}
                for f in newfiles:
                    outfilenames.update(self.make_output_filenames(f))
                # TODO: how to recover from existing output files?
                # Simply re-split? Check whether they all exist and assume 
                # they are correct if they do?
                self.check_files_exist(outfilenames.keys(), "output", 
                                        shouldexist=False)
                outputdir = self.workdir.make_dirname()
                if len(newfiles) <= 10:
                    p = cadoprograms.Duplicates1(*newfiles,
                                                 out=outputdir,
                                                 **self.progparams[0])
                else:
                    filelistname = self.workdir.make_filename("filelist")
                    with open(str(filelistname), "w") as filelistfile:
                        filelistfile.write("\n".join(newfiles))
                    p = cadoprograms.Duplicates1(filelist=filelistname,
                                                 **self.progparams[0])
                (identifier, rc, stdout, stderr, output_files) = \
                        self.submit_command(p, "")
                if rc:
                    raise Exception("Program failed")
                    # Check that the output files exist now
                    # TODO: How to recover from error? Presumably a dup1
                    # process failed, but that should raise a return code
                    # exception
                self.check_files_exist(outfilenames.keys(), "output", 
                                       shouldexist=True)
                current_counts = self.parse_slice_counts(stderr)
                for idx in range(self.nr_slices):
                    self.slice_relcounts[str(idx)] += current_counts[idx]
                update = dict(zip(newfiles, [self.nr_slices] * len(newfiles)))
                self.already_split_input.update(update)
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
                r[str(self.workdir.make_filename(basename, use_subdir=True, subdir=str(i)))] = i
            return r
    
    def parse_slice_counts(self, stderr):
        """ Takes lines of text and looks for slice counts as printed by dup1
        """
        counts = [None] * self.nr_slices
        for line in stderr.decode("ascii").splitlines():
            match = re.match(r'# slice (\d+) received (\d+) relations', line)
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
    
    def request_more_relations(self, target):
        self.send_notification(Notification.WANT_MORE_RELATIONS, target)
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
        return super().paramnames + ("nslices_log",)
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.nr_slices = 2**int(self.params["nslices_log"])
        self.already_done_input = self.make_db_dict(self.make_tablename("infiles"))
        self.slice_relcounts = self.make_db_dict(self.make_tablename("counts"))
        self.slice_relcounts.setdefault(
            None, {str(i): 0 for i in range(0, self.nr_slices)})
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        self.logger.info("Starting")
        input_nrel = 0
        for i in range(0, self.nr_slices):
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            renumber_filename = self.send_request(Request.GET_RENUMBER_FILENAME)
            files = self.send_request(Request.GET_DUP1_FILENAMES, i.__eq__)
            rel_count = self.send_request(Request.GET_DUP1_RELCOUNT, i)
            input_nrel += rel_count
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
            # self.workdir.make_directories(str(i)) OBSOLETE
            if len(files) <= 10:
                p = cadoprograms.Duplicates2(*files,
                                             poly=polyfilename,
                                             rel_count=rel_count * 12 // 10,
                                             renumber=renumber_filename,
                                             **self.progparams[0])
            else:
                filelistname = self.workdir.make_filename("filelist")
                with open(str(filelistname), "w") as filelistfile:
                    filelistfile.write("\n".join(files))
                p = cadoprograms.Duplicates2(poly=polyfilename,
                                             rel_count=rel_count * 12 // 10,
                                             renumber=renumber_filename,
                                             filelist=filelistname,
                                             **self.progparams[0])
            (identifier, rc, stdout, stderr, output_files) = \
                    self.submit_command(p, "")
            if rc:
                raise Exception("Program failed")
            nr_rels = self.parse_remaining(stderr.decode("ascii").splitlines())
            # Mark input file names and output file names
            for f in files:
                self.already_done_input[f] = True
            outfilenames = {f:i for f in files}
            self.add_output_files(outfilenames)
            self.logger.info("%d unique relations remain on slice %d", nr_rels, i)
            self.slice_relcounts[str(i)] = nr_rels
        self.update_ratio(input_nrel, self.get_nrels())
        self.state["last_input_nrel"] = input_nrel
        self.logger.info("%d unique relations remain in total", self.get_nrels())
        self.logger.debug("Exit Duplicates2Task.run(" + self.name + ")")
    
    def parse_remaining(self, text):
        # "     112889 remaining relations"
        for line in text:
            match = re.match(r'\s*(\d+) remaining relations', line)
            if match:
                remaining = int(match.group(1))
                return remaining
        raise Exception("Received no value for remaining relation count")

    def get_nrels(self):
        nrels = 0
        for i in range(0, self.nr_slices):
            nrels += self.slice_relcounts[str(i)]
        return nrels
    
    def update_ratio(self, input_nrel, output_nrel):
        last_input_nrel = self.state.get("last_input_nrel", 0)
        last_output_nrel = self.state.get("last_output_nrel", 0)
        new_in = input_nrel - last_input_nrel
        new_out = output_nrel - last_output_nrel
        if new_in < 0:
            self.logger.error("Negative number %d of new relations?", new_in)
            return
        if new_in == 0:
            return
        if new_out > new_in:
            self.logger.error("More new output relations (%d) than input (%d)?",
                              new_out, new_in)
            return
        ratio = new_out / new_in
        self.logger.info("Of %d newly added relations %d were unique (ratio %f)",
                         new_in, new_out, ratio)
        self.state.update({"last_input_nrel": input_nrel,
            "last_output_nrel": output_nrel, "unique_ratio": ratio})
    
    def request_more_relations(self, target):
        nrels = self.get_nrels()
        if target <= nrels:
            return
        additional_out = target - nrels
        ratio = self.state.get("unique_ratio", 1.)
        additional_in = int(additional_out / ratio)
        newtarget = self.state["last_input_nrel"] + additional_in
        self.logger.info("Got request for %d (%d additional) output relations, estimate %d (%d additional) needed in input",
                         target, additional_out, newtarget, additional_in)
        self.send_notification(Notification.WANT_MORE_RELATIONS, newtarget)
        self.send_notification(Notification.WANT_TO_RUN, None)


class PurgeTask(Task):
    """ Removes singletons and computes excess """
    @property
    def name(self):
        return "purge"
    @property
    def title(self):
        return "Filtering - Singleton removal"
    @property
    def programs(self):
        return (cadoprograms.Purge,)
    @property
    def paramnames(self):
        return super().paramnames
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.state.setdefault("input_nrels", 0)
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        nfree = self.send_request(Request.GET_FREEREL_RELCOUNT)
        nunique = self.send_request(Request.GET_UNIQUE_RELCOUNT)
        minindex = self.send_request(Request.GET_RENUMBER_MININDEX)
        nprimes = self.send_request(Request.GET_RENUMBER_PRIMECOUNT)
        if not nunique:
            raise Exception("No unique relation count received")
        input_nrels = nfree + nunique
        
        if "purgedfile" in self.state and input_nrels == self.state["input_nrels"]:
            self.logger.info("Already have a purged file, and no new input "
                             "relations available. Nothing to do")
            return
        
        self.state.pop("purgedfile", None)
        self.state.pop("input_nrels", None)
        
        self.logger.info("Reading %d unique and %d free relations, total %d"
                         % (nunique, nfree, input_nrels))
        purgedfile = self.workdir.make_filename("purged.gz")
        freerel_filename = self.send_request(Request.GET_FREEREL_FILENAME)
        unique_filenames = self.send_request(Request.GET_UNIQUE_FILENAMES)
        files = unique_filenames + [str(freerel_filename)]
        
        if len(files) <= 10:
            p = cadoprograms.Purge(*files,
                                   nrels=input_nrels, out=purgedfile,
                                   minindex=minindex, nprimes=nprimes,
                                   **self.progparams[0])
        else:
            filelistname = self.workdir.make_filename("filelist")
            with open(str(filelistname), "w") as filelistfile:
                filelistfile.write("\n".join(files))
            p = cadoprograms.Purge(nrels=input_nrels,
                                   out=purgedfile, minindex=minindex,
                                   nprimes=nprimes, filelist=filelistname,
                                   **self.progparams[0])
        (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
        if self.parse_stderr(stderr, input_nrels):
            stats = self.parse_stdout(stdout)
            self.logger.info("After purge, %d relations with %d primes remain "
                             "with weight %s and excess %s", *stats)
            # Update both atomically
            self.state.update({"purgedfile": purgedfile.get_relative(),
                               "input_nrels": input_nrels})
            self.logger.info("Have enough relations")
            self.send_notification(Notification.HAVE_ENOUGH_RELATIONS, None)
        else:
            self.logger.info("Not enough relations")
            self.request_more_relations(nunique)
        self.logger.debug("Exit PurgeTask.run(" + self.name + ")")
    
    def request_more_relations(self, nunique):
        r"""
        We want an excess of $e = 0.1 n_p$, 
        with $n_p$ primes among $n_r$ relations in the output file.
        
        We have $\Delta_r$ additional relations in the output file per
        additional input relation, and $\Delta_p$ additional primes in
        the output file per additional input relation.
        
        We need $a$ additional relations so that
        \begin{eqnarray*}
            n_r + a  \Delta_r - (n_p + a  \Delta_p) & = & 0.1  (n_p + a  \Delta_p) \\
            n_r + a  \Delta_r & = & 1.1  (n_p + a  \Delta_p) \\
            n_r - 1.1 n_p & = & 1.1 a \Delta_p - a \Delta_r \\ 
            \frac{n_r - 1.1 n_p}{1.1 \Delta_p - \Delta_r} & = & a \\
        \end{eqnarray*}
        """
        
        additional = nunique * 0.1
        if "delta_r" in self.state:
            excess = self.state["last_input_nrels"] - \
                self.state["last_input_nprimes"]
            n_r = self.state["last_output_nrels"]
            n_p = self.state["last_output_nprimes"]
            delta_r = self.state["delta_r"]
            delta_p = self.state["delta_p"]
            if abs(1.1 * delta_p - delta_r) > 0.01:
                a = (n_r - 1.1 * n_p) / (1.1 * delta_p - delta_r)
                self.logger.info("a = %f, excess = %d", a, excess)
                # Use the negative excess among the input relations as an
                # upper limit on the additional relations needed, to keep
                # a small donominator from causing a huge value
                if excess < 0 and a > -excess:
                    a = -excess
                if a > 10000. and a < nunique * 0.5:
                    additional = int(a)
        # Always request at least 10k more
        additional = max(additional, 10000)
        
        self.logger.info("Requesting %d additional relations", additional)
        self.send_notification(Notification.WANT_MORE_RELATIONS, nunique + additional)
        self.send_notification(Notification.WANT_TO_RUN, None)
    
    def get_purged_filename(self):
        return self.get_state_filename("purgedfile")
    
    def parse_stderr(self, stderr, input_nrels):
        # If stderr ends with 
        # b'excess < 0.10 * #primes. See -required_excess argument.'
        # then we need more relations from filtering and return False
        input_nprimes = None
        have_enough = True
        # not_enough1 = re.compile(r"excess < (\d+.\d+) \* #primes")
        not_enough1 = re.compile(r"\(excess / nprimes\) = \d+.?\d* < \d+.?\d*. See -required_excess argument.")
        not_enough2 = re.compile(r"number of relations <= number of ideals")
        nrels_nprimes = re.compile(r"\s*nrels=(\d+), nprimes=(\d+); excess=(-?\d+)")
        for line in stderr.decode("ascii").splitlines():
            match = not_enough1.match(line)
            if match:
                have_enough = False
                break
            if not_enough2.match(line):
                have_enough = False
                break
            match = nrels_nprimes.match(line)
            if not match is None:
                (nrels, nprimes, excess) = map(int, match.groups())
                assert nrels - nprimes == excess
                # The first occurrence of the message counts input relations
                if input_nprimes is None:
                    assert input_nrels == nrels
                    input_nprimes = nprimes
        
        # At this point we shoud have:
        # input_nrels, input_nprimes: rels and primes among input
        # nrels, nprimes, excess: rels and primes when purging stopped
        if not input_nprimes is None:
            self.update_excess_per_input(input_nrels, input_nprimes, nrels, nprimes)
        return have_enough
    
    def update_excess_per_input(self, input_nrels, input_nprimes, nrels,
                                nprimes):
        if input_nrels == 0:
            return # Nothing sensible that we can do
        last_input_nrels = self.state.get("last_input_nrels", 0)
        last_input_nprimes = self.state.get("last_input_nprimes", 0)
        last_nrels = self.state.get("last_output_nrels", 0)
        last_nprimes = self.state.get("last_output_nprimes", 0)
        if last_input_nrels >= input_nrels:
            self.logger.warn("Previously stored input nrels (%d) is no "
                             "smaller than value from new run (%d)",
                             last_input_nrels, input_nrels)
            return
        if nrels <= last_nrels:
            self.logger.warn("Previously stored nrels (%d) is no "
                             "smaller than value from new run (%d)",
                             last_nrels, nrels)
            return
        if nprimes <= last_nprimes:
            self.logger.warn("Previously stored nprimes (%d) is no "
                             "smaller than value from new run (%d)",
                             last_nprimes, nprimes)
            return
        self.logger.info("Previous run had %d input relations and ended with "
                         "%d relations and %d primes, new run had %d input "
                         "relations and ended with %d relations and %d primes",
                         last_input_nrels, last_nrels, last_nprimes,
                         input_nrels, nrels, nprimes)
        delta_r = (nrels - last_nrels) / (input_nrels - last_input_nrels)
        delta_p = (nprimes - last_nprimes) / (input_nrels - last_input_nrels)
        self.logger.info("Gained %f output relations and %f primes per input relation",
                         delta_r, delta_p)
        update = {"last_output_nrels": nrels, "last_output_nprimes": nprimes,
                  "last_input_nrels": input_nrels, "last_input_nprimes": input_nprimes,
                  "delta_r": delta_r, "delta_p": delta_p}
        self.state.update(update) 
    
    def parse_stdout(self, stdout):
        # Program stdout is expected in the form:
        #   Final values:
        #   nrels=23105 nprimes=22945 excess=160
        #   weight=382433 weight*nrels=8.84e+09
        # but we allow some extra whitespace
        r = {}
        keys = ("nrels", "nprimes", "weight", "excess")
        for line in stdout.decode("ascii").splitlines():
            for key in keys:
                # Match the key at the start of a line, or after a whitespace
                # Note: (?:) is non-capturing group
                match = re.search(r"(?:^|\s)%s\s*=\s*(\d+)" % key, line)
                if match:
                    if key in r:
                        raise Exception("Found multiple values for %s" % key)
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
        return "merge"
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
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        skip = int(self.progparams[0].get("skip", 32))
        self.progparams[0].setdefault("skip", str(skip))
        self.progparams[0].setdefault("keep", str(skip + 128))

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
            historyfile = self.workdir.make_filename("history")
            progname = cadoprograms.Merge.name
            stdoutfilename = self.workdir.make_filename("%s.stdout" % progname)
            stderrfilename = self.workdir.make_filename("%s.stderr" % progname)
            p = cadoprograms.Merge(mat=purged_filename,
                                   out=historyfile,
                                   stdout=str(stdoutfilename), append_stdout=True,
                                   stderr=str(stderrfilename), append_stderr=True,
                                   **self.progparams[0])
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            if rc:
                raise Exception("Program failed")
            
            indexfile = self.workdir.make_filename("index")
            mergedfile = self.workdir.make_filename("small.bin")
            progname = cadoprograms.Replay.name
            stdoutfilename = self.workdir.make_filename("%s.stdout" % progname)
            stderrfilename = self.workdir.make_filename("%s.stderr" % progname)
            p = cadoprograms.Replay(binary=True,
                                    purged=purged_filename,
                                    history=historyfile, index=indexfile,
                                    out=mergedfile, stdout=str(stdoutfilename),
                                    append_stdout=True,
                                    stderr=str(stderrfilename),
                                    append_stderr=True,
                                    **self.progparams[1])
            (identifier, rc, stdout, stderr, output_files) = \
                    self.submit_command(p, "")
            if rc:
                raise Exception("Program failed")
            
            if not os.path.isfile(str(indexfile)):
                raise Exception("Output file %s does not exist" % indexfile)
            if not os.path.isfile(str(mergedfile)):
                raise Exception("Output file %s does not exist" % mergedfile)
            self.state["indexfile"] = indexfile.get_relative()
            self.state["mergedfile"] = mergedfile.get_relative()
            densefilename = self.workdir.make_filename("small.dense.bin")
            if os.path.isfile(str(densefilename)):
                self.state["densefile"] = densefilename.get_relative()
            
        self.logger.debug("Exit MergeTask.run(" + self.name + ")")
    
    def get_index_filename(self):
        return self.get_state_filename("indexfile")
    
    def get_merged_filename(self):
        return self.get_state_filename("mergedfile")

    
    def get_dense_filename(self):
        return self.get_state_filename("densefile")
    

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
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if not "dependency" in self.state:
            self.logger.info("Starting")
            workdir = self.workdir.make_dirname()
            self.workdir.make_directories()
            mergedfile = self.send_request(Request.GET_MERGED_FILENAME)
            progname = cadoprograms.BWC.name
            stdoutfilename = self.workdir.make_filename("%s.stdout" % progname)
            stderrfilename = self.workdir.make_filename("%s.stderr" % progname)
            matrix = os.path.realpath(str(mergedfile))
            wdir = os.path.realpath(str(workdir))
            p = cadoprograms.BWC(complete=True,
                                 matrix=matrix,  wdir=wdir, nullspace="left", 
                                 stdout=str(stdoutfilename), append_stdout=True, 
                                 stderr=str(stderrfilename), append_stderr=True,
                                 **self.progparams[0])
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            if rc:
                raise Exception("Program failed")
            dependencyfilename = self.workdir.make_filename("W", use_subdir = True)
            if not os.path.isfile(str(dependencyfilename)):
                raise Exception("Kernel file %s does not exist" % dependencyfilename)
            self.state["dependency"] =  dependencyfilename.get_relative()
        self.logger.debug("Exit LinAlgTask.run(" + self.name + ")")

    def get_dependency_filename(self):
        return self.get_state_filename("dependency")
    
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
        return super().paramnames + ()

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if not "kernel" in self.state:
            self.logger.info("Starting")
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            kernelfilename = self.workdir.make_filename("kernel")
            
            purgedfilename = self.send_request(Request.GET_PURGED_FILENAME)
            indexfilename = self.send_request(Request.GET_INDEX_FILENAME)
            densefilename = self.send_request(Request.GET_DENSE_FILENAME)
            dependencyfilename = self.send_request(Request.GET_DEPENDENCY_FILENAME)
            
            progname = cadoprograms.Characters.name
            stdoutfilename = self.workdir.make_filename("%s.stdout" % progname)
            stderrfilename = self.workdir.make_filename("%s.stderr" % progname)
            p = cadoprograms.Characters(poly=polyfilename,
                    purged=purgedfilename, index=indexfilename,
                    wfile=dependencyfilename, out=kernelfilename,
                    heavyblock=densefilename, stdout=str(stdoutfilename),
                    append_stdout=True, stderr=str(stderrfilename),
                    append_stderr=True, **self.progparams[0])
            (identifier, rc, stdout, stderr, output_files) = \
                    self.submit_command(p, "")
            if rc:
                raise Exception("Program failed")
            if not os.path.isfile(str(kernelfilename)):
                raise Exception("Output file %s does not exist" % kernelfilename)
            self.state["kernel"] = kernelfilename.get_relative()
        self.logger.debug("Exit CharactersTask.run(" + self.name + ")")
        return
    
    def get_kernel_filename(self):
        return self.get_state_filename("kernel")


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
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.factors = self.make_db_dict(self.make_tablename("factors"))
        self.add_factor(int(self.params["N"]))
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        if not self.is_done():
            self.logger.info("Starting")
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            purgedfilename = self.send_request(Request.GET_PURGED_FILENAME)
            indexfilename = self.send_request(Request.GET_INDEX_FILENAME)
            kernelfilename = self.send_request(Request.GET_KERNEL_FILENAME)
            prefix = self.send_request(Request.GET_LINALG_PREFIX)
            p = cadoprograms.Sqrt(ab = True,
                    poly=polyfilename, purged=purgedfilename,
                    index=indexfilename, kernel=kernelfilename,
                    prefix=prefix, **self.progparams[0])
            (identifier, rc, stdout, stderr, output_files) = self.submit_command(p, "")
            if rc:
                raise Exception("Program failed")
            
            while not self.is_done():
                self.state.setdefault("next_dep", 0)
                dep = self.state["next_dep"]
                p = cadoprograms.Sqrt(ab=False, rat=True,
                        alg=True, gcd=True, dep=dep, poly=polyfilename,
                        purged=purgedfilename, index=indexfilename,
                        kernel=kernelfilename, prefix=prefix,
                        **self.progparams[0])
                (identifier, rc, stdout, stderr, output_files) = \
                    self.submit_command(p, "dep%d" % self.state["next_dep"])
                if rc:
                    raise Exception("Program failed")
                if not stdout.decode("ascii").strip() == "Failed":
                    factorlist = list(map(int, stdout.decode("ascii").split()))
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
    
    def __init__(self, address, port, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.used_ids = {}
        self.pids = self.make_db_dict(self.make_tablename("client_pids"))
        self.hosts = self.make_db_dict(self.make_tablename("client_hosts"))
        assert set(self.pids) == set(self.hosts)
        # Invariants: the keys of self.pids and of self.hosts are the same set.
        # The keys of self.used_ids are a subset of the keys of self.pids.
        # A clientid is in self.used_ids if we know that clientid to be
        # currently running.
        
        if 'scriptpath' in self.params:
            self.progparams[0]['execpath'] = self.params['scriptpath']
        self.server = "http://%s:%d" % (address, port)
        
        # If hostnames are of the form @file, read host names from file,
        # one host name per line
        match = re.match(r"@(.*)", self.params["hostnames"])
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
            self.launch_one_client(host.strip())
        running_clients = [(cid, self.hosts[cid], pid) for (cid, pid) in self.pids.items()]
        s = ", ".join(["%s (Host %s, PID %d)" % t for t in running_clients])
        self.logger.info("Running clients: %s" % s)
        # Check for old clients which we did not mean to start this run
        for cid in set(self.pids) - set(self.used_ids):
            if self.is_alive(cid):
                self.logger.warn("Client id %s (Host %s, PID %d), launched "
                                 "in a previous run and not meant to be "
                                 "launched this time, is still running",
                                 cid, self.hosts[cid], self.pids[cid])
            else:
                self.logger.warn("Client id %s (Host %s, PID %d), launched "
                                 "in a previous run and not meant to be "
                                 "launched this time, seems to have died. "
                                 "I'll forget about this client.",
                                 cid, self.hosts[cid], self.pids[cid])
                del(self.hosts[cid])
                del(self.pids[cid])
    
    def make_unique_id(self, host):
        # Make a unique client id for host
        clientid = host
        i = 1
        while clientid in self.used_ids:
            assert clientid in self.pids
            assert clientid in self.hosts
            i += 1
            clientid = "%s%d" % (host, i)
        return clientid
    
    # Cases:
    # Client was never started. Start it, add to state
    # Client was started, but does not exist any more. Remove from state, then start and add again
    # Client was started, and does still exists. Nothing to do.
    
    def launch_one_client(self, host, clientid = None):
        if clientid is None:
            clientid = self.make_unique_id(host)
        # Check if client is already running
        if clientid in self.pids:
            assert self.hosts[clientid] == host
            if self.is_alive(clientid):
                self.logger.info("Client %s on host %s with PID %d already running",
                                 clientid, host, self.pids[clientid])
                self.used_ids[clientid] = True
                return
            else:
                self.logger.info("Client %s on host %s with PID %d seems to have died",
                                 clientid, host, self.pids[clientid])
                del(self.pids[clientid])
                del(self.hosts[clientid])
        
        self.logger.info("Starting client id %s on host %s", clientid, host)
        wuclient = cadoprograms.WuClient(server=self.server,
                                         clientid=clientid, daemon=True,
                                         **self.progparams[0])
        if host == "localhost":
            process = cadocommand.Command(wuclient)
        else:
            process = cadocommand.RemoteCommand(wuclient, host, self.get_parameters(), 
                                                self.get_param_prefix())
        (rc, stdout, stderr) = process.wait()
        if rc != 0:
            self.logger.warning("Starting client on host %s failed.", host)
            if stdout:
                self.logger.warning("Stdout: %s", stdout.decode("ASCII").strip())
            if stderr:
                self.logger.warning("Stdout: %s", stderr.decode("ASCII").strip())
            return
        match = None
        if not stdout is None:
            match = re.match(r"PID: (\d+)", stdout.decode("ascii"))
        if not match:
            self.logger.warning("Client did not print PID")
            if not stdout is None:
                self.logger.warning("Stdout: %s", stdout.decode("ASCII").strip())
            if not stderr is None:
                self.logger.warning("Stdout: %s", stderr.decode("ASCII").strip())
            return
        self.used_ids[clientid] = True
        self.pids[clientid] = int(match.group(1))
        self.hosts[clientid] = host

    def kill_all_clients(self):
        # Need the list() to make a copy as dict will change in loop body
        for clientid in list(self.pids):
            (rc, stdout, stderr) = self.kill_client(clientid)
            if rc == 0:
                self.logger.info("Stopped client %s (Host %s, PID %d)",
                                 clientid, self.hosts[clientid], self.pids[clientid])
                del(self.pids[clientid])
                del(self.hosts[clientid])
            else:
                self.logger.warning("Stopping client %s (Host %s, PID %d) failed",
                                    clientid, self.hosts[clientid], self.pids[clientid])
                if stdout:
                    self.logger.warning("Stdout: %s", stdout.decode("ASCII").strip())
                if stderr:
                    self.logger.warning("Stdout: %s", stderr.decode("ASCII").strip())
                # Assume that the client is already dead and remove it from
                # the list of running clients
                del(self.pids[clientid])
                del(self.hosts[clientid])
    
    def kill_client(self, clientid, signal = None):
        pid = self.pids[clientid]
        host = self.hosts[clientid]
        kill = cadoprograms.Kill(pid, signal=signal)
        process = cadocommand.RemoteCommand(kill, host, self.parameters, self.get_param_prefix())
        return process.wait()

class Message(object):
    def __init__(self, sender, key, value = None):
        self.sender = sender
        self.key = key
        self.value = value
    def get_sender(self):
        return self.sender
    def get_key(self):
        return self.key
    def get_value(self):
        return self.value
    @classmethod
    def reverse_lookup(cls, reference):
        for key in dir(cls):
            if getattr(cls, key) == reference:
                return key


class Notification(Message):
    FINISHED_POLYNOMIAL_SELECTION = object()
    WANT_MORE_RELATIONS = object()
    HAVE_ENOUGH_RELATIONS = object()
    REGISTER_FILENAME = object()
    UNREGISTER_FILENAME = object()
    SUBMIT_WU = object()
    VERIFY_WU = object()
    WANT_TO_RUN = object()
    SUBSCRIBE_WU_NOTIFICATIONS = object()
    CHECK_TIMEDOUT_WUS = object()

class Request(Message):
    # Lacking a proper enum before Python 3.3, we generate dummy objects
    # which have separate identity and can be used as dict keys
    GET_POLYNOMIAL = object()
    GET_POLYNOMIAL_FILENAME = object()
    GET_FACTORBASE_FILENAME = object()
    GET_FREEREL_FILENAME = object()
    GET_RENUMBER_FILENAME = object()
    GET_FREEREL_RELCOUNT = object()
    GET_RENUMBER_PRIMECOUNT = object()
    GET_RENUMBER_MININDEX = object()
    GET_SIEVER_FILENAMES = object()
    GET_SIEVER_RELCOUNT = object()
    GET_DUP1_FILENAMES = object()
    GET_DUP1_RELCOUNT = object()
    GET_UNIQUE_RELCOUNT = object()
    GET_UNIQUE_FILENAMES = object()
    GET_PURGED_FILENAME = object()
    GET_MERGED_FILENAME = object()
    GET_INDEX_FILENAME = object()
    GET_DENSE_FILENAME = object()
    GET_DEPENDENCY_FILENAME = object()
    GET_LINALG_PREFIX = object()
    GET_KERNEL_FILENAME = object()
    GET_WU_RESULT = object()
    GET_NR_AVAILABLE_WU = object()

class CompleteFactorization(wudb.DbAccess, cadoparams.UseParameters, patterns.Mediator):
    """ The complete factorization, aggregate of the individual tasks """
    @property
    def name(self):
        return "tasks"
    
    CAN_CANCEL_WUS = 0
    
    def __init__(self, db, parameters, path_prefix):
        super().__init__(db = db, parameters = parameters, path_prefix = path_prefix)
        self.logger = logging.getLogger("Complete Factorization")
        self.params = self.myparams(("name", "workdir"))
        self.db_listener = self.make_db_listener()
        self.registered_filenames = self.make_db_dict(self.params["name"] + '_server_registered_filenames')
        self.chores = []
        
        # Init WU DB
        self.wuar = self.make_wu_access()
        self.wuar.create_tables()
        
        # Set up WU server
        serverparams = parameters.myparams(["address", "port"], path_prefix + ["server"])
        serveraddress = serverparams.get("address", socket.gethostname())
        serverport = int(serverparams.get("port", 8001))
        uploaddir = self.params["workdir"].rstrip(os.sep) + os.sep + self.params["name"] + ".upload/"
        threaded = False
        self.server = wuserver.ServerLauncher(serveraddress, serverport, threaded, self.get_db_filename(), 
            self.registered_filenames, uploaddir, bg=True, only_registered=True)
        
        # Init client lists
        self.clients = []
        for (path, key) in self.get_parameters().find(['slaves'], 'hostnames'):
            self.clients.append(StartClientsTask(serveraddress, serverport, 
                                                 mediator = self,
                                                 db = db, 
                                                 parameters = parameters, 
                                                 path_prefix = path))
        
        parampath = self.get_param_path()
        sievepath = parampath + ['sieve']
        filterpath = parampath + ['filter']
        linalgpath = parampath + ['linalg']
        
        self.polysel = PolyselTask(mediator = self,
                                   db = db, 
                                   parameters = parameters, 
                                   path_prefix = parampath)
        self.fb = FactorBaseTask(mediator = self,
                                 db = db, 
                                 parameters = parameters, 
                                 path_prefix = sievepath)
        self.freerel = FreeRelTask(mediator = self,
                                   db = db, 
                                   parameters = parameters, 
                                   path_prefix = sievepath)
        self.sieving = SievingTask(mediator = self,
                                   db = db, 
                                   parameters = parameters, 
                                   path_prefix = sievepath)
        self.dup1 = Duplicates1Task(mediator = self,
                                    db = db, 
                                    parameters = parameters, 
                                    path_prefix = filterpath)
        self.dup2 = Duplicates2Task(mediator = self,
                                    db = db, 
                                    parameters = parameters, 
                                    path_prefix = filterpath)
        self.purge = PurgeTask(mediator = self,
                               db = db, 
                               parameters = parameters, 
                               path_prefix = filterpath)
        self.merge = MergeTask(mediator = self,
                               db = db, 
                               parameters = parameters, 
                               path_prefix = filterpath)
        self.linalg = LinAlgTask(mediator = self,
                                 db = db, 
                                 parameters = parameters, 
                                 path_prefix = linalgpath)
        self.characters = CharactersTask(mediator = self,
                                         db = db, 
                                         parameters = parameters, 
                                         path_prefix = linalgpath)
        self.sqrt = SqrtTask(mediator = self,
                             db = db, 
                             parameters = parameters, 
                             path_prefix = parampath)
        
        # Defines an order on tasks in which tasks that want to run should be run
        self.tasks = (self.polysel, self.fb, self.freerel, self.sieving,
                      self.dup1, self.dup2, self.purge, self.merge,
                      self.linalg, self.characters, self.sqrt)
        
        # Assume that all tasks want to run. Ff they are finished already, 
        # they will just return immediately
        self.tasks_that_want_to_run = list(self.tasks)
        
        self.request_map = {
            Request.GET_POLYNOMIAL: self.polysel.get_poly,
            Request.GET_POLYNOMIAL_FILENAME: self.polysel.get_poly_filename,
            Request.GET_FACTORBASE_FILENAME: self.fb.get_filename,
            Request.GET_FREEREL_FILENAME: self.freerel.get_freerel_filename,
            Request.GET_RENUMBER_FILENAME: self.freerel.get_renumber_filename,
            Request.GET_FREEREL_RELCOUNT: self.freerel.get_nrels,
            Request.GET_RENUMBER_PRIMECOUNT: self.freerel.get_nprimes,
            Request.GET_RENUMBER_MININDEX: self.freerel.get_minindex,
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
            Request.GET_KERNEL_FILENAME: self.characters.get_kernel_filename,
            Request.GET_WU_RESULT: self.db_listener.send_result,
            Request.GET_NR_AVAILABLE_WU: self.wuar.count_available
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
                # self.logger.info("Next task that wants to run: %s", task.title)
                self.tasks_that_want_to_run.remove(task)
                task.run()
                return True
        return False
    
    def do_chores(self):
        while self.chores:
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
    
    def check_timedout_wus(self, sender):
        # Check at most once per minute
        attr = "last_check_timedout_wus"
        if getattr(self, attr, 0.) + 60. > time.time():
            # TODO
            # self.wuar.timeout_wus(self.params.get("wu_timeout", 3600))
            pass
        setattr(self, attr, time.time())
    
    def relay_notification(self, message):
        assert isinstance(message, Notification)
        sender = message.get_sender()
        key = message.get_key()
        value = message.get_value()
        self.logger.message("Received notification from %s, key = %s, value = %s",
                            sender, Notification.reverse_lookup(key), value)
        """ The relay for letting Tasks talk to us and each other """
        if key is Notification.WANT_MORE_RELATIONS:
            if sender is self.purge:
                self.dup2.request_more_relations(value)
            elif sender is self.dup2:
                self.dup1.request_more_relations(value)
            elif sender is self.dup1:
                self.sieving.request_more_relations(value)
            else:
                raise Exception("Got WANT_MORE_RELATIONS from unknown sender")
        elif key is Notification.HAVE_ENOUGH_RELATIONS:
            if sender is self.purge:
                # TODO: cancel only sieving WUs?
                self.chores.append(self.CAN_CANCEL_WUS)
            else:
                raise Exception("Got HAVE_ENOUGH_RELATIONS from unknown sender")
        elif key is Notification.REGISTER_FILENAME:
            if isinstance(sender, ClientServerTask):
                self.register_filename(value)
            else:
                raise Exception("Got REGISTER_FILENAME, but not from a ClientServerTask")
        elif key is Notification.SUBMIT_WU:
            if isinstance(sender, ClientServerTask):
                self.add_wu(value)
            else:
                raise Exception("Got SUBMIT_WU, but not from a ClientServerTask")
        elif key is Notification.VERIFY_WU:
            if isinstance(sender, ClientServerTask):
                self.verify_wu(value)
            else:
                raise Exception("Got VERIFY_WU, but not from a ClientServerTask")
        elif key is Notification.WANT_TO_RUN:
            if sender in self.tasks_that_want_to_run:
                raise Exception("Got request from %s to run, but it was in run queue already",
                                sender)
            else:
                self.tasks_that_want_to_run.append(sender)
        elif key is Notification.SUBSCRIBE_WU_NOTIFICATIONS:
            if isinstance(sender, ClientServerTask):
                return self.db_listener.subscribeObserver(sender)
            else:
                raise Exception("Got SUBSCRIBE_WU_NOTIFICATIONS, but not from a ClientServerTask")
        elif key is Notification.CHECK_TIMEDOUT_WUS:
            if isinstance(sender, ClientServerTask):
                return self.check_timedout_wus(sender)
            else:
                raise Exception("Got CHECK_TIMEDOUT_WUS, but not from a ClientServerTask")
        else:
            raise KeyError("Notification from %s has unknown key %s" % (sender, key))
    
    def answer_request(self, request):
        assert isinstance(request, Request)
        sender = request.get_sender()
        key = request.get_key()
        value = request.get_value()
        self.logger.message("Received request from %s, key = %s, values = %s",
                            sender, Request.reverse_lookup(key), value)
        if not key in self.request_map:
            raise KeyError("Unknown Request key %s from sender %s" %
                           (key, sender))
        if value is None:
            return self.request_map[key]()
        else:
            return self.request_map[key](value)
    
    def handle_message(self, message):
        if isinstance(message, Notification):
            self.relay_notification(Notification)
        elif isinstance(message, Request):
            return self.answer_request(message)
        else:
            raise TypeError("Message is neither Notification nor Request")
