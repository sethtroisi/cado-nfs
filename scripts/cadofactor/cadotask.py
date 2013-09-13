#!/usr/bin/env python3
import re
import os.path
from fractions import gcd
import abc
import random
import time
import datetime
from collections import OrderedDict
from itertools import zip_longest
from math import log, sqrt
import logging
import socket
import patterns
import wudb
import cadoprograms
import cadoparams
import cadocommand
import wuserver
from workunit import Workunit

# Some parameters are provided by the param file but can change during
# the factorization, like rels_wanted. On one hand, we want automatic
# updates to parameters to be stored in the DB, otoh, we want to allow
# externally setting new parameters. Need to distinguish between new
# external parameters that overwrite DB, and old external parameters
# that don't overwrite. Or maybe two ways to specify external params:
# --defaults which does not overwrite, and --forceparam which does

# Pattern for floating-point number
re_fp = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
cap_fp = "(%s)" % re_fp

class PolynomialParseException(Exception):
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
        self.MurphyE = 0.
        poly = {}
        for line in lines:
            # print ("Parsing line: >%s<" % line)
            # If this is a comment line telling the Murphy E value,
            # extract the value and store it
            match = re.match(r"\s*#\s*MurphyE\s*\(.*\)=(.*)$", line)
            if match:
                self.MurphyE = float(match.group(1))
                continue
            # Drop comment, strip whitespace
            line2 = line.split('#', 1)[0].strip()
            # If nothing is left, process next line
            if not line2:
                continue
            # All remaining lines must be of the form "x: y"
            array = line2.split(":")
            if not len(array) == 2:
                raise PolynomialParseException("Invalid line %s" % line)
            key = array[0].strip()
            value = array[1].strip()
            if not key in dict(self.keys):
                raise PolynomialParseException("Invalid key %s in line %s" %
                                (key, line))
            poly[key] = value
        for (key, isrequired) in self.keys:
            if isrequired and not key in poly:
                raise PolynomialParseException("Key %s missing" % key)
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
    
    def setE(self, MurphyE):
        self.MurphyE = float(MurphyE)
    
    def create_file(self, filename, params):
        # Write polynomial to a file, and add lines with parameters such as
        # "alim" if supplied in params
        with open(str(filename), "w") as poly_file:
            poly_file.write(str(self))
            for key in self.paramnames:
                if key in params:
                    poly_file.write(key + ": %s\n" % params[key])


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
    def get_wdir_relative(self):
        return self.filepath
    def isfile(self):
        return os.path.isfile(str(self))
    def isdir(self):
        return os.path.isdir(str(self))
    def mkdir(self):
        os.mkdir(str(self))
    def realpath(self):
        return os.path.realpath(str(self))
    def open(self, *args, **kwargs):
        return open(str(self), *args, **kwargs)


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
        If use_subdir is True and subdir is a string, make a filename of the
        form
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
        dirname = self.make_dirname()
        if not dirname.isdir():
            dirname.mkdir()
        if subdirs:
            for subdir in subdirs:
                dirname = self.make_dirname(subdir)
                if not dirname.isdir():
                    dirname.mkdir()
        return


class Statistics(object):
    """ Class that holds statistics on program execution, and can merge two
    such statistics.
    """
    def __init__(self, conversions):
        self.conversions = conversions
        self.stats = {}
    
    @staticmethod
    def typecast(values, types):
        """ Cast the values in values to the types specified in types """
        return [t(v) for (v, t) in zip(values, types)]
    
    @staticmethod
    def _to_str(stat):
        """ Convert one statistic to a string """
        return " ".join(map(str, stat))
    
    @staticmethod
    def _from_str(string, types):
        """ Convert a string (probably from a state dict) to a statistic """
        return Statistics.typecast(string.split(), types)
    
    def from_dict(self, stats):
        """ Initialise values in self from the strings in the "stats"
        dictionary
        """
        for conversion in self.conversions:
            (msgfmt, key, types, defaults, combine, regex) = conversion
            if key in stats:
                assert not key in self.stats
                self.stats[key] = self._from_str(stats.get(key, defaults),
                                                 types)
                assert not self.stats[key] is None
    
    def parse_line(self, line):
        """ Parse one line of program output and look for statistics.
        
        If they are found, they are added to self.stats.
        """
        for conversion in self.conversions:
            (msgfmt, key, types, defaults, combine, regex) = conversion
            match = regex.match(line)
            if match:
                assert not key in self.stats
                # print (pattern.pattern, match.groups())
                self.stats[key] = self.typecast(match.groups(), types)
                assert not self.stats[key] is None
    
    def merge_one_stat(self, key, new_val, combine):
        if key in self.stats:
            self.stats[key] = combine(self.stats[key], new_val)
        else:
            self.stats[key] = new_val
        assert not self.stats[key] is None
        # print(self.stats)
    
    def merge_stats(self, new_stats):
        """ Merge the stats currently in self with the Statistics in
        "new_stats"
        """
        
        assert self.conversions == new_stats.conversions
        for conversion in self.conversions:
            (msgfmt, key, types, defaults, combine, regex) = conversion
            if key in new_stats.stats:
                self.merge_one_stat(key, new_stats.stats[key], combine)
    
    def as_dict(self):
        return {key:self._to_str(self.stats[key]) for key in self.stats}
    
    def as_strings(self):
        result = []
        for conversion in self.conversions:
            (msgfmt, key, types, defaults, combine, regex) = conversion
            if key in self.stats:
                if len(self.stats[key]) == len(types):
                    result.append(msgfmt % tuple(self.stats[key]))
        return result
    
    # Helper functions for processing statistics.
    # We can't make them @staticmethod or references are not callable
    def add_list(*lists):
        """ Add zero or more lists elementwise.
        
        Short lists are handled as if padded with zeroes.
        
        >>> Statistics.add_list([])
        []
        >>> Statistics.add_list([1])
        [1]
        >>> Statistics.add_list([1,2], [3,7])
        [4, 9]
        >>> Statistics.add_list([1,2], [3,7], [5], [3,1,4,1,5])
        [12, 10, 4, 1, 5]
        """
        return [sum(items) for items in zip_longest(*lists, fillvalue=0)]
    
    def weigh(samples, weights):
        return [sample*weight for (sample, weight) in zip(samples, weights)]
    
    def combine_mean(means, samples):
        """ From two lists, one containing values and the other containing
        the respective sample sizes (i.e., weights of the values), compute
        the combined mean (i.e. the weighted mean of the values).
        The two lists must have equal length.
        """
        assert len(means) == len(samples)
        total_samples = sum(samples)
        weighted_sum = sum(Statistics.weigh(means, samples))
        return [weighted_sum / total_samples, total_samples]
    
    def zip_combine_mean(*lists):
        """ From a list of 2-tuples, each tuple containing a value and a
        weight, compute the weighted average of the values.
        """
        for l in lists:
            assert len(l) == 2
        (means, samples) = zip(*lists)
        return Statistics.combine_mean(means, samples)
    
    def combine_stats(*stats):
        """ Computes the combined mean and std.dev. for the stats
        
        stats is a list of 3-tuples, each containing number of sample points,
        mean, and std.dev.
        Returns a 3-tuple with the combined number of sample points, mean,
        and std. dev.
        """
        
        # Samples is a list containing the first item (number of samples) of
        # each item of stats, means is list of means, stddevs is list of
        # std. dev.s
        for s in stats:
            assert len(s) == 3
        
        (samples, means, stddevs) = zip(*stats)
        
        (total_mean, total_samples) = Statistics.combine_mean(means, samples)
        # t is the E[X^2] part of V(X)=E(X^2) - (E[X])^2
        t = [mean**2 + stddev**2 for (mean, stddev) in zip(means, stddevs)]
        # Compute combined variance
        total_var = Statistics.combine_mean(t, samples)[0] - total_mean**2
        return [total_samples, total_mean, sqrt(total_var)]
    
    def test_combine_stats():
        """ Test function for combine_stats()
        
        >>> Statistics.test_combine_stats()
        True
        """
        
        from random import randrange
        
        def mean(x):
            return float(sum(x))/float(len(x))
        def var(x):
            E = mean(x)
            return mean([(a-E)**2 for a in x])
        def stddev(x):
            return sqrt(var(x))
        
        # Generate between 1 and 5 random integers in [1,100]
        lengths = [randrange(100) + 1 for i in range(randrange(5) + 1)]
        lengths = [1, 10]
        # Generate lists of random integers in [1,100]
        lists = [[randrange(100) for i in range(l)] for l in lengths]
        stats = [(length, mean(l), stddev(l))
            for (length, l) in zip(lengths, lists)]
        
        combined = []
        for l in lists:
            combined += l
        
        combined1 = Statistics.combine_stats(*stats)
        combined2 = [len(combined), mean(combined), stddev(combined)]
        if abs(combined1[2] - combined2[2]) > 0.2 * combined2[2]:
            print("lists = %r" % lists)
            print("combineds = %r" % combined)
            print("stats = %r" % stats)
            print("combined1 = %r" % combined1)
            print("combined2 = %r" % combined2)
            print(combined1[2], combined2[2])
            print(abs(combined1[2] / combined2[2] - 1))
        return combined1[0] == combined2[0] and \
                abs(combined1[1] / combined2[1] - 1) < 1e-10 and \
                abs(combined1[2] - combined2[2]) <= 1e-10 * combined2[2]
    
    def smallest_10(*lists):
        concat = []
        for l in lists:
            concat += l
        concat.sort()
        return concat[0:10]


class HasName(object, metaclass=abc.ABCMeta):
    @abc.abstractproperty
    def name(self):
        # The name of the task in a simple form that can be used as
        # a Python dictionary key, a directory name, part of a file name,
        # part of an SQL table name, etc. That pretty much limits it to
        # alphabetic first letter, and alphanumeric rest.
        pass

class HasTitle(object, metaclass=abc.ABCMeta):
    @abc.abstractproperty
    def title(self):
        # A pretty name for the task, will be used in screen output
        pass

class DoesLogging(HasTitle, metaclass=abc.ABCMeta):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger = logging.getLogger(self.title)

class MakesTablenames(HasName):
    @property
    def tablename_prefix(self):
        """ Prefix string for table names
        
        By default, the table name prefix is the name attribute, but this can
        be overridden
        """
        return self.name
    
    @staticmethod
    def check_tablename(name):
        no_ = name.replace("_", "")
        if not no_[0].isalpha() or not no_[1:].isalnum():
            raise Exception("%s is not valid for an SQL table name" % name)
    
    def make_tablename(self, extra = None):
        """ Return a name for a DB table """
        # Maybe replace SQL-disallowed characters here, like digits and '.' ?
        # Could be tricky to avoid collisions
        name = self.tablename_prefix
        if extra:
            name = name + '_' + extra
        self.check_tablename(name)
        return name

class HasState(MakesTablenames, wudb.HasDbConnection):
    """ Declares that the class has a DB-backed dictionary in which the class
    can store state information.
    
    The dicatonary is available as an instance attribute "state".
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        name = self.make_tablename()
        self.state = self.make_db_dict(name, connection=self.db_connection)


class FilesCreator(MakesTablenames, wudb.HasDbConnection, metaclass=abc.ABCMeta):
    """ A base class for classes that produce a list of output files, with
    some auxiliary information stored with each file (e.g., nr. of relations).
    This info is stored in the form of a DB-backed dictionary, with the file
    name as the key and the auxiliary data as the value.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        tablename = self.make_tablename("outputfiles")
        self.output_files = self.make_db_dict(tablename,
                                              connection=self.db_connection)
    
    def add_output_files(self, filenames, *, commit):
        """ Adds a dict of files to the list of existing output files """
        for filename in filenames:
            if filename in self.output_files:
                raise KeyError("%s already in output files table" % filename)
        self.output_files.update(filenames, commit=commit)
    
    def get_output_filenames(self, condition=None):
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
    
    def forget_output_filenames(self, filenames, *, commit):
        self.output_files.clear(filenames, commit=commit)


class HasStatistics(object, metaclass=abc.ABCMeta):
    @abc.abstractproperty
    def stat_conversions(self):
        pass

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.statistics = Statistics(self.stat_conversions)
    def get_statistics_as_strings(self):
        return self.statistics.as_strings()


class Task(patterns.Colleague, HasState, cadoparams.UseParameters,
           DoesLogging, metaclass=abc.ABCMeta):
    """ A base class that represents one task that needs to be processed.
    
    Sub-classes must define class variables:
    """
    
    # Properties that subclasses need to define
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
    @property
    def param_nodename(self):
        return self.name
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        ''' Sets up a database connection and a DB-backed dictionary for
        parameters. Reads parameters from DB, and merges with hierarchical
        parameters in the parameters argument. Parameters passed in by
        parameters argument do not override values in the DB-backed
        parameter dictionary.
        '''
        
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.logger.debug("Enter Task.__init__(%s)",
                          self.name)
        self.logger.debug("state = %s", self.state)
        # Set default parameters for this task, if any are given
        self.params = self.parameters.myparams(self.paramnames)
        self.logger.debug("self.parameters = %s", self.parameters)
        self.logger.debug("params = %s", self.params)
        # Set default parameters for our programs
        self.progparams = []
        for prog in self.programs:
            progparams = self.parameters.myparams(prog.get_accepted_keys(),
                                                  prog.name)
            self.progparams.append(progparams)
        # FIXME: whether to init workdir or not should not be controlled via
        # presence of a "workdir" parameter, but by class definition
        if "workdir" in self.params:
            self.workdir = WorkDir(self.params["workdir"], self.params["name"],
                               self.name)
        self.init_stats()
        self.logger.debug("Exit Task.__init__(%s)", self.name)
        return
    
    def translate_input_filename(self, filename):
        return filename

    def test_outputfile_exists(self, filename):
        return filename.isfile()
    
    @staticmethod
    def check_files_exist(filenames, filedesc, shouldexist):
        """ Check that the output files in "filenames" exist or don't exist,
        according to shouldexist.
        
        Raise IOError if any check fails, return None
        """
        for f in filenames:
            if isinstance(f, FilePath):
                exists = f.isfile()
            else:
                exists = os.path.isfile(f)
            if shouldexist and not exists:
                raise IOError("%s file %s does not exist" % (filedesc, f))
            elif not shouldexist and exists:
                raise IOError("%s file %s already exists" % (filedesc, f))
        return
    
    # These two function go together, one produces a workunit name from the
    # name of the factorization, the task name, and a task-provided identifier,
    # and the other function splits them again
    wu_paste_char = '_'
    def make_wuname(self, identifier, attempt=None):
        assert not self.wu_paste_char in self.params["name"]
        assert not self.wu_paste_char in self.name
        assert not self.wu_paste_char in identifier
        if attempt is None:
            return self.wu_paste_char.join([self.params["name"], self.name,
                                            identifier])
        else:
            assert isinstance(attempt, int)
            return self.wu_paste_char.join([self.params["name"], self.name,
                                            identifier, str(attempt)])
    
    def split_wuname(self, wuname):
        arr = wuname.split(self.wu_paste_char)
        if len(arr) == 3:
            arr.append(None)
        else:
            arr[3] = int(arr[3])
        return arr
    
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
    
    def submit_command(self, command, identifier, commit=True):
        ''' Run a command.
        Return the result tuple. If the caller is an Observer, also send
        result to updateObserver().
        '''
        wuname = self.make_wuname(identifier)
        process = cadocommand.Command(command)
        (rc, stdout, stderr) = process.wait()
        message = Task.ResultInfo(wuname, rc, stdout, stderr,
                                  command.get_output_files())
        if isinstance(self, patterns.Observer):
            # pylint: disable=E1101
            self.updateObserver(message)
        return message
    
    def filter_notification(self, message):
        wuid = message.get_wu_id()
        rc = message.get_exitcode(0)
        stdout = message.get_stdout(0)
        stderr = message.get_stderr(0)
        output_files = message.get_output_files()
        self.logger.message("%s: Received notification for wuid=%s, rc=%d, "
                            "output_files=[%s]",
                            self.name, wuid, rc, ", ".join(output_files))
        if rc != 0:
            self.logger.error("Return code is: %d", rc)
        if stdout:
            self.logger.debug("stdout is: %s", stdout)
        if stderr:
            self.logger.debug("stderr is: %s", stderr)
        if output_files:
            self.logger.message("Output files are: %s", ", ".join(output_files))
        (name, task, identifier, attempt) = self.split_wuname(wuid)
        if name != self.params["name"] or task != self.name:
            # This notification is not for me
            self.logger.message("Notification %s is not for me", wuid)
            return
        self.logger.message("Notification %s is for me", wuid)
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
    
    def verification(self, message, ok, *, commit):
        pass

    def get_state_filename(self, key):
        if not key in self.state:
            return None
        return self.workdir.path_in_workdir(self.state[key])

    def make_std_paths(self, progname, do_increment=True):
        count = self.state.get("stdiocount", 0)
        if do_increment:
            count += 1
        did_increment = do_increment
        while True:
            try:
                stdoutname = "%s.stdout.%d" % (progname, count)
                stderrname = "%s.stderr.%d" % (progname, count)
                self.check_files_exist((stdoutname, stderrname), "stdio",
                                       shouldexist=False)
            except IOError:
                count += 1
                did_increment = True
                self.logger.warning("Stdout or stderr files with index %d "
                                    "already exist", count)
            else:
                break
        stdoutpath = self.workdir.make_filename(stdoutname)
        stderrpath = self.workdir.make_filename(stderrname)
        if did_increment:
            self.state["stdiocount"] = count
        return (stdoutpath, stderrpath)
    
    def init_stats(self):
        if not isinstance(self, HasStatistics):
            return
        self.statistics.from_dict(self.state)
    
    def print_stats(self):
        if not isinstance(self, HasStatistics):
            return
        stat_msgs =  self.get_statistics_as_strings()
        if stat_msgs:
            self.logger.info("Aggregate statistics:")
            for msg in stat_msgs:
                self.logger.info(msg)
    
    def parse_stats(self, filename, *, commit):
        if not isinstance(self, HasStatistics):
            return
        new_stats = Statistics(self.stat_conversions)
        with open(filename, "r") as inputfile:
            for line in inputfile:
                new_stats.parse_line(line)
        self.logger.debug("Newly arrived stats: %s", new_stats.as_dict())
        self.statistics.merge_stats(new_stats)
        update = self.statistics.as_dict()
        self.logger.debug("Combined stats: %s", update)
        self.state.update(update, commit=commit)


class ClientServerTask(Task, wudb.UsesWorkunitDb, patterns.Observer):
    @abc.abstractproperty
    def paramnames(self):
        return super().paramnames + ("maxwu", "wutimeout")
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.state.setdefault("wu_submitted", 0)
        self.state.setdefault("wu_received", 0)
        self.params.setdefault("maxwu", 10)
        assert self.get_number_outstanding_wus() >= 0
        self.send_notification(Notification.SUBSCRIBE_WU_NOTIFICATIONS, None)
    
    def submit_command(self, command, identifier, commit=True):
        ''' Submit a workunit to the database. '''
        
        while self.get_number_available_wus() >= self.params["maxwu"]:
            self.wait()
        wuid = self.make_wuname(identifier)
        wutext = command.make_wu(wuid)
        for filename in command.get_exec_files() + command.get_input_files():
            basename = os.path.basename(filename)
            self.send_notification(Notification.REGISTER_FILENAME,
                                   {basename:filename})
        
        self.logger.info("Adding workunit %s to database", wuid)
        # print ("WU:\n%s" % wutext)
        key = "wu_submitted"
        self.state.update({key: self.state[key] + 1}, commit=False)
        self.wuar.create(wutext, commit=commit)
    
    def verification(self, message, ok, *, commit):
        wuid = message.get_wu_id()
        ok_str = "ok" if ok else "not ok"
        self.logger.info("Marking workunit %s as %s", wuid, ok_str)
        assert self.get_number_outstanding_wus() >= 1
        key = "wu_received"
        self.state.update({key: self.state[key] + 1}, commit=False)
        self.wuar.verification(message.get_wu_id(), ok, commit=commit)
    
    def cancel_available_wus(self):
        self.logger.info("Cancelling remaining workunits")
        self.wuar.cancel_all_available()
    
    def get_number_outstanding_wus(self):
        return self.state["wu_submitted"] - self.state["wu_received"]

    def get_number_available_wus(self):
        return self.wuar.count_available()

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
            self.resubmit_timed_out_wus()
            time.sleep(1)
    
    def resubmit_one_wu(self, wu, commit=True):
        """ Takes a Workunit instance and adds it to workunits table under
        a modified name.
        """
        wuid = wu.get_id()
        (name, task, identifier, attempt) = self.split_wuname(wuid)
        attempt = 2 if attempt is None else attempt + 1
        new_wuid = self.make_wuname(identifier, attempt)
        wu.set_id(new_wuid)
        self.logger.info("Resubmitting workunit %s as %s", wuid, new_wuid)
        self.wuar.create(str(wu), commit=commit)

    def resubmit_timed_out_wus(self):
        # We don't store the lastcheck in state as we do *not* want to check
        # instantly when we start up - clients should get a chance to upload
        # results first
        now = time.time()
        if not hasattr(self, "last_timeout_check"):
            self.logger.debug("Setting last timeout check to %f", now)
            self.last_timeout_check = now
            return
        
        check_every = 60 # Check every xx seconds
        if self.last_timeout_check + check_every >= now:
            # self.logger.info("It's not time to check yet, now = %f", now)
            return
        self.last_timeout_check = now
        
        timeout = int(self.params.get("wutimeout", 10800))
        delta = datetime.timedelta(seconds=timeout)
        cutoff = str(datetime.datetime.utcnow() - delta)
        self.logger.debug("Doing timeout check, cutoff=%s, and setting last check to %f",
                          cutoff, now)
        results = self.wuar.query(eq={"status":1}, lt={"timeassigned": cutoff})
        if not results:
            self.logger.debug("Found no timed-out workunits")
        for entry in results:
            self.resubmit_one_wu(Workunit(entry["wu"]), commit=False)
            self.wuar.cancel(entry["wuid"], commit=True)


class PolyselTask(ClientServerTask, HasStatistics, patterns.Observer):
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
    # Stat: potential collisions=124.92 (2.25e+00/s)
    # Stat: raw lognorm (nr/min/av/max/std): 132/18.87/21.83/24.31/0.48
    # Stat: optimized lognorm (nr/min/av/max/std): 125/20.10/22.73/24.42/0.69
    # Stat: av. g0/adm2 ratio: 8.594e+04
    # Stat: tried 83 ad-value(s), found 132 polynomial(s), 125 below maxnorm
    # Stat: best logmu: 20.10 21.05 21.41 21.48 21.51 21.57 21.71 21.74 21.76 21.76
    # Stat: total phase took 55.47s
    # Stat: rootsieve took 54.54s
    @property
    def stat_conversions(self):
        return (
            ("potential collisions: %f",
             "stats_collisions",
             (float,),
             "0",
             Statistics.add_list,
             re.compile(r"# Stat: potential collisions=%s" % cap_fp)
            ),
            ("raw lognorm (nr/min/av/max/std): %d/%f/%f/%f/%f",
             "stats_rawlognorm",
             (int, float, float, float, float),
             "0 0 0 0 0",
             PolyselTask.update_lognorms,
             re.compile(r"# Stat: raw lognorm \(nr/min/av/max/std\): (\d+)/%s/%s/%s/%s" % ((cap_fp,) * 4))
            ),
            ("optimized lognorm (nr/min/av/max/std): %d/%f/%f/%f/%f",
             "stats_optlognorm",
             (int, float, float, float, float),
             "0 0 0 0 0",
             PolyselTask.update_lognorms,
             re.compile(r"# Stat: optimized lognorm \(nr/min/av/max/std\): (\d+)/%s/%s/%s/%s" % ((cap_fp,) * 4))
            ),
            ("tried ad-value(s): %d, found polynomial(s): %d, " \
                "below maxnorm: %d",
             "stats_tries",
             (int, )*3,
             "0 0 0",
             Statistics.add_list,
             re.compile(r"# Stat: tried (\d+) ad-value\(s\), found (\d+) polynomial\(s\), (\d+) below maxnorm")
            ),
            # Note for "best logmu" pattern: a regex like (%s )* does not work;
            # the number of the capture group is determined by the parentheses
            # in the regex string, so trying to repeat a group like this will
            # always capture to the *same* group, overwriting previous matches,
            # so that in the end, only the last match is in the capture group.
            ("10 best logmu: %g %g %g %g %g %g %g %g %g %g",
             "stats_logmu",
             (float, )*10,
             "",
             Statistics.smallest_10,
             re.compile(r"# Stat: best logmu: %s %s %s %s %s %s %s %s %s %s"
                        % ((cap_fp, )*10))
            ),
            ("total time: %f",
             "stats_total_time",
             (float,),
             "0",
             Statistics.add_list,
             re.compile(r"# Stat: total phase took %ss" % cap_fp)
            ),
            ("rootsieve time: %f",
             "rootsieve_time",
             (float,),
             "0",
             Statistics.add_list,
             re.compile(r"# Stat: rootsieve took %ss" % cap_fp)
            )
        )
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.state["adnext"] = \
            max(self.state.get("adnext", 0), self.params.get("admin", 0))
        self.bestpoly = None
        if "bestpoly" in self.state:
            self.bestpoly = Polynomial(self.state["bestpoly"].splitlines())
            self.bestpoly.setE(self.state["bestE"])
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.name, self.state)
        
        if not self.bestpoly is None:
            self.logger.info("Best polynomial previously found in %s has "
                             "Murphy_E = %g",
                             self.state["bestfile"], self.bestpoly.MurphyE)
        else:
            self.logger.info("No polynomial was previously found")
        
        if self.is_done():
            self.logger.info("Polynomial selection already finished - "
                             "nothing to do")
            # If the poly file got lost somehow, write it again
            filename = self.get_state_filename("polyfilename")
            if filename is None or not filename.isfile():
                self.logger.warn("Polynomial file disappeared, writing again")
                self.write_poly_file()
            return True
        
        # Submit all the WUs we need to reach admax
        while self.need_more_wus():
            self.submit_one_wu()
        
        # Wait for all the WUs to finish
        while self.get_number_outstanding_wus() > 0:
            self.wait()
        
        if self.bestpoly is None:
            self.logger.error ("No polynomial found. Consider increasing the "
                               "search range bound admax, or maxnorm")
            return False
        self.logger.info("Finished, best polynomial from file %s has Murphy_E "
                         "= %g", self.state["bestfile"] , self.bestpoly.MurphyE)
        self.write_poly_file()
        return True
    
    def is_done(self):
        return not self.bestpoly is None and not self.need_more_wus() and \
            self.get_number_outstanding_wus() == 0
    
    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return
        (filename, ) = message.get_output_files()
        try:
            poly = self.parse_poly(filename)
        except PolynomialParseException as e:
            self.logger.error("Invalid polyselect file %s: %s",
                              filename, e)
            self.verification(message, False, commit=True)
            return
        if not poly is None:
            self.bestpoly = poly
            update = {"bestE": poly.MurphyE, "bestpoly": str(poly),
                      "bestfile": filename}
            self.state.update(update, commit=False)
        self.parse_stats(filename, commit=False)
        self.verification(message, True, commit=True)
    
    def parse_poly(self, filename):
        poly = None
        with open(filename, "r") as polyfile:
            for line in polyfile:
                # A "# WARNING:" line can occur when a bad gcc version was
                # used for compilation
                if re.match("# WARNING", line):
                    self.logger.warn("File %s contains: %s",
                                     filename, line.strip())
                # If there is a "No polynomial found" message in the file,
                # we just skip it. If there is no polynomial in this file,
                # we'll reach EOF next and get the poly==None case below.
                # If the file happens to be the concatenation of several
                # polyselect output files, we should keep looking.
                if re.match("No polynomial found", line):
                    continue
                # If we get the "Best polynomial" marker, we stop reading
                # the file, and let Polynomial.__init__() parse the rest
                if re.match("# Best polynomial found:", line):
                    break
            try:
                poly = Polynomial(polyfile)
            except PolynomialParseException as e:
                self.logger.error("Invalid polyselect file %s: %s",
                                  filename, e)
                return None
        
        if not poly or not poly.is_valid():
            self.logger.info('No polynomial found in %s', filename)
            return None
        if not poly.MurphyE:
            self.logger.warn("Polynomial in file %s has no Murphy E value",
                             filename)
        if self.bestpoly is None or poly.MurphyE > self.bestpoly.MurphyE:
            self.logger.info("New best polynomial from file %s:"
                             " Murphy E = %g" % (filename, poly.MurphyE))
            self.logger.debug("New best polynomial is:\n%s", poly)
            return poly
        else:
            self.logger.info("Best polynomial from file %s with E=%g is "
                             "no better than current best with E=%g",
                             filename, poly.MurphyE, self.bestpoly.MurphyE)
        return None
    
    def update_lognorms(old_lognorm, new_lognorm):
        lognorm = [0, 0, 0, 0, 0]
        # print("update_lognorms: old_lognorm: %s" % old_lognorm)
        # print("update_lognorms: new_lognorm: %s" % new_lognorm)
        # New minimum. Don't use default value of 0 for minimum
        lognorm[1] = min(old_lognorm[1] or new_lognorm[1], new_lognorm[1])
        # New maximum
        lognorm[3] = max(old_lognorm[3], new_lognorm[3])
        # Rest is done by combine_stats(). [0::2] selects indices 0,2,4
        lognorm[0::2] = Statistics.combine_stats(old_lognorm[0::2],
                                                 new_lognorm[0::2])
        return lognorm
    
    def write_poly_file(self):
        filename = self.workdir.make_filename("poly")
        self.bestpoly.create_file(filename, self.params)
        self.state["polyfilename"] = filename.get_wdir_relative()
    
    def get_poly(self):
        if not "bestpoly" in self.state:
            return None
        return Polynomial(self.state["bestpoly"].splitlines())
    
    def get_poly_filename(self):
        return self.get_state_filename("polyfilename")

    def need_more_wus(self):
        return self.state["adnext"] < self.params["admax"]
    
    def submit_one_wu(self):
        adstart = self.state["adnext"]
        adend = adstart + self.params["adrange"]
        adend = min(adend, self.params["admax"])
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
            self.submit_command(p, "%d-%d" % (adstart, adend), commit=False)
        self.state.update({"adnext": adend}, commit=True)

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
            ("alim", "gzip")


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
            raise Exception("FactorBaseTask(): no polynomial "
                            "received from PolyselTask")
        
        # Check if we have already computed the target file for this polynomial
        if "poly" in self.state:
            prevpoly = Polynomial(self.state["poly"].splitlines())
            if poly != prevpoly:
                if "outputfile" in self.state:
                    self.logger.info("Received different polynomial, "
                                     "discarding old one")
                    del(self.state["outputfile"])
                self.state["poly"] = str(poly)
        else:
            self.state["poly"] = str(poly)
        
        if not "outputfile" in self.state:
            self.logger.info("Starting")
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            
            # Make file name for factor base/free relations file
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params.get("gzip", True) else ""
            outputfilename = self.workdir.make_filename("roots" + use_gz)

            # Run command to generate factor base/free relations file
            p = cadoprograms.MakeFB(poly=polyfilename,
                                    stdout = str(outputfilename),
                                    **self.progparams[0])
            message = self.submit_command(p, "")
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            
            self.state["outputfile"] = outputfilename.get_wdir_relative()
            self.logger.info("Finished")

        self.check_files_exist([self.get_filename()], "output",
                               shouldexist=True)
        return True
    
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
            ("alim", "rlim", "gzip")
    wanted_regex = {
        'nfree': (r'# Free relations: (\d+)', int),
        'nprimes': (r'Renumbering struct: nprimes=(\d+)', int)
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
            raise Exception("FreerelTask(): no polynomial "
                            "received from PolyselTask")
        
        # Check if we have already computed the target file for this polynomial
        if "poly" in self.state:
            prevpoly = Polynomial(self.state["poly"].splitlines())
            if poly != prevpoly:
                if "freerelfilename" in self.state:
                    self.logger.info("Received different polynomial, "
                                     "discarding old one")
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
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params.get("gzip", True) else ""
            freerelfilename = self.workdir.make_filename("freerel" + use_gz)
            renumberfilename = self.workdir.make_filename("renumber" + use_gz)

            # Run command to generate factor base/free relations file
            p = cadoprograms.FreeRel(poly=polyfilename,
                                     renumber=renumberfilename,
                                     out=str(freerelfilename),
                                     **self.progparams[0])
            message = self.submit_command(p, "")
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            stderr = message.get_stderr(0)
            found = self.parse_file(stderr.decode("ascii").splitlines())
            self.state.update(found)
            self.logger.info("Found %d free relations" % self.state["nfree"])
            
            self.state["freerelfilename"] = freerelfilename.get_wdir_relative()
            self.state["renumberfilename"] = renumberfilename.get_wdir_relative()
            self.logger.info("Finished")

        self.check_files_exist([self.get_freerel_filename(),
                                self.get_renumber_filename()], "output",
                               shouldexist=True)
        return True

    def parse_file(self, text):
        found = {}
        for line in text:
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


class SievingTask(ClientServerTask, FilesCreator, HasStatistics,
                  patterns.Observer):
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
    @property
    def stat_conversions(self):
        # Average J=1017 for 168 special-q's, max bucket fill 0.737035
        # Total wct time 7.0s [precise timings available only for mono-thread]
        # Total 26198 reports [0.000267s/r, 155.9r/sq]
        return (
            (
                "Average J: %f for %d special-q",
                "stats_avg_J",
                (float, int),
                "0 0",
                Statistics.zip_combine_mean,
                re.compile(r"# Average J=%s for (\d+) special-q's" % cap_fp)
            ),
            (
                "Max bucket fill: %f",
                "stats_max_bucket_fill",
                (float, ),
                "0",
                max,
                re.compile(r"#.*max bucket fill %s" % cap_fp)
            ),
            (
                "Total wall clock time: %f",
                "stats_total_wall_clock_time",
                (float, ),
                "0",
                Statistics.add_list,
                re.compile(r"# Total wct time %ss" % cap_fp)
            )
        )
    # We seek to this many bytes before the EOF to look for the
    # "Total xxx reports" message
    file_end_offset = 1000
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        qmin = self.params.get("qmin", 0)
        if "qnext" in self.state:
            self.state["qnext"] = max(self.state["qnext"], qmin)
        else:
            self.state["qnext"] = self.params.get("qmin", self.params["alim"])
        
        self.state.setdefault("rels_found", 0)
        self.state["rels_wanted"] = self.params.get("rels_wanted", 0)
        self.params.setdefault("maxwu", "10")
        if self.state["rels_wanted"] == 0:
            # TODO: Choose sensible default value
            pass
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        while self.state["rels_found"] < self.state["rels_wanted"]:
            q0 = self.state["qnext"]
            q1 = q0 + self.params["qrange"]
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params.get("gzip", True) else ""
            outputfilename = \
                self.workdir.make_filename("%d-%d%s" % (q0, q1, use_gz))
            self.check_files_exist([outputfilename], "output",
                                   shouldexist=False)
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            factorbase = self.send_request(Request.GET_FACTORBASE_FILENAME)
            p = cadoprograms.Las(q0=q0, q1=q1,
                                 poly=polyfilename, factorbase=factorbase,
                                 out=outputfilename, stats_stderr = True,
                                 **self.progparams[0])
            self.submit_command(p, "%d-%d" % (q0, q1), commit=False)
            self.state.update({"qnext": q1}, commit=True)
        self.logger.info("Reached target of %d relations, now have %d",
                         self.state["rels_wanted"], self.state["rels_found"])
        self.logger.debug("Exit SievingTask.run(" + self.name + ")")
        return True
    
    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return
        output_files = message.get_output_files()
        assert len(output_files) == 1
        stderrfilename = message.get_stderr(0)
        if  stderrfilename is None:
            self.logger.error("No stderr output received for workunit %s,"
                              "ignoring (relations file: %s)",
                              identifier, output_files[0])
            self.verification(message, False, commit=True)
            return
        rels = self.parse_rel_count(stderrfilename)
        if rels is None:
            self.verification(message, False, commit=True)
            self.logger.error("Number of relations message not found in "
                              "file %s", stderrfilename)
            return
        self.state["rels_found"] += rels
        self.add_output_files({output_files[0]: rels}, commit=False)
        self.parse_stats(stderrfilename, commit=False)
        self.verification(message, True, commit=True)
        self.logger.info("Found %d relations in %s, total is now %d",
                         rels, output_files[0], self.state["rels_found"])
    
    def parse_rel_count(self, filename):
        size = os.path.getsize(filename)
        with open(filename, "r") as f:
            f.seek(max(size - self.file_end_offset, 0))
            for line in f:
                match = re.match(r"# Total (\d+) reports ", line)
                if match:
                    rels = int(match.group(1))
                    return rels
        return None
    
    def get_statistics_as_strings(self):
        strings = ["Total number of relations: %d" % self.get_nrels()]
        strings += super().get_statistics_as_strings()
        return strings
    
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
                             self.state["rels_wanted"],
                             self.state["rels_found"])
        else:
            self.logger.info("New goal for number of relations is %d, but "
                             "already have %d. No need to sieve more",
                             self.state["rels_wanted"],
                             self.state["rels_found"])

class Duplicates1Task(Task, FilesCreator, HasStatistics):
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
    @property
    def stat_conversions(self):
        # "End of read: 229176 relations in 0.9s -- 21.0 MB/s -- 253905.7 rels/s"
        # Without leading "# " !
        return (
            (
                "CPU time for dup1: %fs",
                "stats_dup1_time",
                (float, ),
                "0",
                Statistics.add_list,
                re.compile(r"End of read: \d+ relations in %ss" % cap_fp)
            ),
        )
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.nr_slices = 2**self.params["nslices_log"]
        tablename = self.make_tablename("infiles")
        self.already_split_input = self.make_db_dict(tablename,
                                                     connection=self.db_connection)
        self.slice_relcounts = self.make_db_dict(self.make_tablename("counts"),
                                                 connection=self.db_connection)
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
                # for the relation count in the files
                # TODO: pass a list or generator expression in the request
                # here?
                total = self.slice_relcounts["0"]
                for f in newfiles:
                    total += self.send_request(Request.GET_SIEVER_RELCOUNT, f)
                self.slice_relcounts.update({"0": total}, commit=False)
                update1 = dict.fromkeys(newfiles, self.nr_slices)
                self.already_split_input.update(update1, commit=False)
                update2 = dict.fromkeys(newfiles, 0)
                self.add_output_files(update2, commit=True)
            else:
                outputdir = self.workdir.make_dirname()
                run_counter = self.state.get("run_counter", 0)
                prefix = "dup1.%s" % run_counter
                (stdoutpath, stderrpath) = \
                        self.make_std_paths(cadoprograms.Duplicates1.name)
                if len(newfiles) <= 10:
                    p = cadoprograms.Duplicates1(*newfiles,
                                                 prefix=prefix,
                                                 out=outputdir,
                                                 stdout=str(stdoutpath),
                                                 stderr=str(stderrpath),
                                                 **self.progparams[0])
                else:
                    filelistname = self.workdir.make_filename("filelist")
                    with filelistname.open("w") as filelistfile:
                        filelistfile.write("\n".join(newfiles) + "\n")
                    p = cadoprograms.Duplicates1(filelist=filelistname,
                                                 prefix=prefix,
                                                 out=outputdir,
                                                 stdout=str(stdoutpath),
                                                 stderr=str(stderrpath),
                                                 **self.progparams[0])
                message = self.submit_command(p, "")
                if message.get_exitcode(0) != 0:
                    raise Exception("Program failed")
                    # Check that the output files exist now
                    # TODO: How to recover from error? Presumably a dup1
                    # process failed, but that should raise a return code
                    # exception
                assert message.get_stderr(0) is None
                with stderrpath.open("r") as stderrfile:
                    stderr = stderrfile.read()
                outfilenames = self.parse_output_files(stderr)
                self.logger.debug("Output file names: %s", outfilenames)
                self.check_files_exist(outfilenames.keys(), "output",
                                       shouldexist=True)
                self.state.update({"run_counter": run_counter + 1},
                                  commit=False)
                current_counts = self.parse_slice_counts(stderr)
                self.parse_stats(str(stderrpath), commit=False)
                # Add relation count from the newly processed files to the
                # relations-per-slice dict
                update1 = {str(idx): self.slice_relcounts[str(idx)] + 
                                     current_counts[idx]
                           for idx in range(self.nr_slices)}
                self.slice_relcounts.update(update1, commit=False)
                # Add the newly processed input files to the list of already
                # processed input files
                update2 = dict.fromkeys(newfiles, self.nr_slices)
                self.already_split_input.update(update2, commit=False)
                # Add the newly produced output files and commit everything
                self.add_output_files(outfilenames, commit=True)
        totals = ["%d: %d" % (i, self.slice_relcounts[str(i)])
                  for i in range(0, self.nr_slices)]
        self.logger.info("Relations per slice: %s", ", ".join(totals))
        self.logger.debug("Exit Duplicates1Task.run(" + self.name + ")")
        return True
    
    @staticmethod
    def parse_output_files(stderr):
        files = {}
        for line in stderr.splitlines():
            match = re.match(r'# Opening output file for slice (\d+) : (.+)$',
                             line)
            if match:
                (slicenr, filename) = match.groups()
                files[filename] = int(slicenr)
        return files
    
    def parse_slice_counts(self, stderr):
        """ Takes lines of text and looks for slice counts as printed by dup1
        """
        counts = [None] * self.nr_slices
        for line in stderr.splitlines():
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

class Duplicates2Task(Task, FilesCreator, HasStatistics):
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
        return super().paramnames + ("nslices_log", "alim", "rlim")
    @property
    def stat_conversions(self):
        # "End of read: 229176 relations in 0.9s -- 21.0 MB/s -- 253905.7 rels/s"
        # Without leading "# " !
        return (
            (
                "CPU time for dup2: %fs",
                "stats_dup2_time",
                (float, ),
                "0",
                Statistics.add_list,
                re.compile(r"End of read: \d+ relations in %ss" % cap_fp)
            ),
        )
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.nr_slices = 2**self.params["nslices_log"]
        tablename = self.make_tablename("infiles")
        self.already_done_input = self.make_db_dict(tablename, connection=self.db_connection)
        self.slice_relcounts = self.make_db_dict(self.make_tablename("counts"), connection=self.db_connection)
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
            self.forget_output_filenames(self.get_output_filenames(i.__eq__),
                                         commit=True)
            del(self.slice_relcounts[str(i)])
            name = "%s.slice%d" % (cadoprograms.Duplicates2.name, i)
            (stdoutpath, stderrpath) = \
                self.make_std_paths(name, do_increment=(i == 0))
            if len(files) <= 10:
                p = cadoprograms.Duplicates2(*files,
                                             poly=polyfilename,
                                             rel_count=rel_count,
                                             renumber=renumber_filename,
                                             stdout=str(stdoutpath),
                                             stderr=str(stderrpath),
                                             **self.progparams[0])
            else:
                filelistname = self.workdir.make_filename("filelist")
                with filelistname.open("w") as filelistfile:
                    filelistfile.write("\n".join(files) + "\n")
                p = cadoprograms.Duplicates2(poly=polyfilename,
                                             rel_count=rel_count,
                                             renumber=renumber_filename,
                                             filelist=filelistname,
                                             stdout=str(stdoutpath),
                                             stderr=str(stderrpath),
                                             **self.progparams[0])
            message = self.submit_command(p, "")
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            assert message.get_stderr(0) is None
            with stderrpath.open("r") as stderrfile:
                nr_rels = self.parse_remaining(stderrfile)
            # Mark input file names and output file names
            for f in files:
                self.already_done_input[f] = True
            outfilenames = {f:i for f in files}
            self.add_output_files(outfilenames, commit=True)
            # Disabled for now, there are multiple lines of the same format
            # which we can't parse atm.
            # self.parse_stats(str(stderrpath))
            self.logger.info("%d unique relations remain on slice %d",
                             nr_rels, i)
            self.slice_relcounts[str(i)] = nr_rels
        self.update_ratio(input_nrel, self.get_nrels())
        self.state["last_input_nrel"] = input_nrel
        self.logger.info("%d unique relations remain in total",
                         self.get_nrels())
        self.logger.debug("Exit Duplicates2Task.run(" + self.name + ")")
        return True
    
    @staticmethod
    def parse_remaining(text):
        # "     112889 remaining relations"
        for line in text:
            match = re.match(r'At the end:\s*(\d+) remaining relations', line)
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
        self.logger.info("Of %d newly added relations %d were unique "
                         "(ratio %f)", new_in, new_out, ratio)
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
        self.logger.info("Got request for %d (%d additional) output relations, "
                         "estimate %d (%d additional) needed in input",
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
        return super().paramnames + ("alim", "rlim")
    
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
        minlim = min(self.params["alim"], self.params["rlim"])
        minindex = int(2. * minlim / (log(minlim) - 1))
        nprimes = self.send_request(Request.GET_RENUMBER_PRIMECOUNT)
        if not nunique:
            raise Exception("No unique relation count received")
        input_nrels = nfree + nunique
        
        if "purgedfile" in self.state and \
                input_nrels == self.state["input_nrels"]:
            self.logger.info("Already have a purged file, and no new input "
                             "relations available. Nothing to do")
            return True
        
        self.state.pop("purgedfile", None)
        self.state.pop("input_nrels", None)
        
        self.logger.info("Reading %d unique and %d free relations, total %d"
                         % (nunique, nfree, input_nrels))
        purgedfile = self.workdir.make_filename("purged.gz")
        freerel_filename = self.send_request(Request.GET_FREEREL_FILENAME)
        unique_filenames = self.send_request(Request.GET_UNIQUE_FILENAMES)
        files = unique_filenames + [str(freerel_filename)]
        (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.Purge.name)
        
        if len(files) <= 10:
            p = cadoprograms.Purge(*files,
                                   nrels=input_nrels, out=purgedfile,
                                   minindex=minindex, nprimes=nprimes,
                                   stdout=str(stdoutpath),
                                   stderr=str(stderrpath),
                                   **self.progparams[0])
        else:
            filelistname = self.workdir.make_filename("filelist")
            with filelistname.open("w") as filelistfile:
                filelistfile.write("\n".join(files))
            p = cadoprograms.Purge(nrels=input_nrels,
                                   out=purgedfile, minindex=minindex,
                                   nprimes=nprimes, filelist=filelistname,
                                   stdout=str(stdoutpath),
                                   stderr=str(stderrpath),
                                   **self.progparams[0])
        message = self.submit_command(p, "")
        assert message.get_stdout(0) is None
        assert message.get_stderr(0) is None
        with stderrpath.open("r") as stderrfile:
            stderr = stderrfile.read()
        with stdoutpath.open("r") as stdoutfile:
            stdout = stdoutfile.read()
        if self.parse_stderr(stderr, input_nrels):
            stats = self.parse_stdout(stdout)
            self.logger.info("After purge, %d relations with %d primes remain "
                             "with weight %s and excess %s", *stats)
            self.state.update({"purgedfile": purgedfile.get_wdir_relative(),
                               "input_nrels": input_nrels})
            self.logger.info("Have enough relations")
            self.send_notification(Notification.HAVE_ENOUGH_RELATIONS, None)
        else:
            self.logger.info("Not enough relations")
            self.request_more_relations(nunique)
        self.logger.debug("Exit PurgeTask.run(" + self.name + ")")
        return True
    
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
        self.send_notification(Notification.WANT_MORE_RELATIONS,
                               nunique + additional)
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
        not_enough1 = re.compile(r"\(excess / nprimes\) = \d+.?\d* < \d+.?\d*. "
                                 r"See -required_excess argument.")
        not_enough2 = re.compile(r"number of relations <= number of ideals")
        nrels_nprimes = re.compile(r"\s*nrels=(\d+), nprimes=(\d+); "
                                   r"excess=(-?\d+)")
        for line in stderr.splitlines():
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
            self.update_excess_per_input(input_nrels, input_nprimes, nrels,
                                         nprimes)
        return have_enough
    
    def update_excess_per_input(self, input_nrels, input_nprimes, nrels,
                                nprimes):
        if input_nrels == 0:
            return # Nothing sensible that we can do
        last_input_nrels = self.state.get("last_input_nrels", 0)
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
        self.logger.info("Gained %f output relations and %f primes per input "
                         "relation", delta_r, delta_p)
        update = {"last_output_nrels": nrels, "last_output_nprimes": nprimes,
                  "last_input_nrels": input_nrels, "delta_r": delta_r,
                  "delta_p": delta_p}
        self.state.update(update)
    
    def parse_stdout(self, stdout):
        # Program stdout is expected in the form:
        #   Final values:
        #   nrels=23105 nprimes=22945 excess=160
        #   weight=382433 weight*nrels=8.84e+09
        # but we allow some extra whitespace
        r = {}
        keys = ("nrels", "nprimes", "weight", "excess")
        for line in stdout.splitlines():
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
            ("skip", "forbw", "coverNmax", "keep", "maxlevel", "ratio", "gzip")
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        skip = self.progparams[0].get("skip", 32)
        self.progparams[0].setdefault("skip", skip)
        self.progparams[0].setdefault("keep", skip + 128)

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
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params.get("gzip", True) else ""
            historyfile = self.workdir.make_filename("history" + use_gz)
            (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.Merge.name)
            p = cadoprograms.Merge(mat=purged_filename,
                                   out=historyfile,
                                   stdout=str(stdoutpath),
                                   stderr=str(stderrpath),
                                   **self.progparams[0])
            message = self.submit_command(p, "")
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            
            indexfile = self.workdir.make_filename("index" + use_gz)
            mergedfile = self.workdir.make_filename("small.bin")
            (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.Replay.name)
            p = cadoprograms.Replay(binary=True,
                                    purged=purged_filename,
                                    history=historyfile, index=indexfile,
                                    out=mergedfile, stdout=str(stdoutpath),
                                    stderr=str(stderrpath),
                                    **self.progparams[1])
            message = self.submit_command(p, "")
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            
            if not indexfile.isfile():
                raise Exception("Output file %s does not exist" % indexfile)
            if not mergedfile.isfile():
                raise Exception("Output file %s does not exist" % mergedfile)
            self.state["indexfile"] = indexfile.get_wdir_relative()
            self.state["mergedfile"] = mergedfile.get_wdir_relative()
            densefilename = self.workdir.make_filename("small.dense.bin")
            if densefilename.isfile():
                self.state["densefile"] = densefilename.get_wdir_relative()
            
        self.logger.debug("Exit MergeTask.run(" + self.name + ")")
        return True
    
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
            (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.BWC.name)
            matrix = mergedfile.realpath()
            wdir = workdir.realpath()
            p = cadoprograms.BWC(complete=True,
                                 matrix=matrix,  wdir=wdir, nullspace="left",
                                 stdout=str(stdoutpath),
                                 stderr=str(stderrpath),
                                 **self.progparams[0])
            message = self.submit_command(p, "")
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            dependencyfilename = self.workdir.make_filename("W", use_subdir = True)
            if not dependencyfilename.isfile():
                raise Exception("Kernel file %s does not exist" % dependencyfilename)
            self.state["dependency"] =  dependencyfilename.get_wdir_relative()
        self.logger.debug("Exit LinAlgTask.run(" + self.name + ")")
        return True

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
            
            (stdoutpath, stderrpath) = \
                    self.make_std_paths(cadoprograms.Characters.name)
            p = cadoprograms.Characters(poly=polyfilename,
                    purged=purgedfilename, index=indexfilename,
                    wfile=dependencyfilename, out=kernelfilename,
                    heavyblock=densefilename, stdout=str(stdoutpath),
                    stderr=str(stderrpath),
                    **self.progparams[0])
            message = self.submit_command(p, "")
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            if not kernelfilename.isfile():
                raise Exception("Output file %s does not exist" % kernelfilename)
            self.state["kernel"] = kernelfilename.get_wdir_relative()
        self.logger.debug("Exit CharactersTask.run(" + self.name + ")")
        return True
    
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
        self.factors = self.make_db_dict(self.make_tablename("factors"), connection=self.db_connection)
        self.add_factor(self.params["N"])
    
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
            message = self.submit_command(p, "")
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            
            while not self.is_done():
                self.state.setdefault("next_dep", 0)
                dep = self.state["next_dep"]
                self.logger.info("Trying dependency %d", dep)
                p = cadoprograms.Sqrt(ab=False, rat=True,
                        alg=True, gcd=True, dep=dep, poly=polyfilename,
                        purged=purgedfilename, index=indexfilename,
                        kernel=kernelfilename, prefix=prefix,
                        **self.progparams[0])
                message = self.submit_command(p, "dep%d" % self.state["next_dep"])
                if message.get_exitcode(0) != 0:
                    raise Exception("Program failed")
                stdout = message.get_stdout(0)
                lines = stdout.decode("ascii").splitlines()
                # Skip last factor which cannot produce a new split on top
                # of what the smaller factors did
                for line in lines[:-1]:
                    if line == "Failed":
                        break
                    self.add_factor(int(line))
                self.state["next_dep"] += 1
            self.logger.info("finished")
        self.logger.info("Factors: %s" % " ".join(self.get_factors()))
        self.logger.debug("Exit SqrtTask.run(" + self.name + ")")
        return True
    
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
        if number <= 3:
            return number >= 2
        if number % 2 == 0:
            return False
        for i in range(0, passes):
            # random.randrange(n) produces random integer in [0, n-1].
            # We want [2, n-2]
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
    @property
    def param_nodename(self):
        return None
    
    def __init__(self, address, port, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator = mediator, db = db, parameters = parameters,
                         path_prefix = path_prefix)
        self.used_ids = {}
        self.pids = self.make_db_dict(self.make_tablename("client_pids"), connection=self.db_connection)
        self.hosts = self.make_db_dict(self.make_tablename("client_hosts"), connection=self.db_connection)
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
                self.hosts_to_launch = [line.strip() for line in f]
        else:
            self.hosts_to_launch = [host.strip() for host in
                    self.params["hostnames"].split(",")]

        if "nrclients" in self.params:
            self.hosts_to_launch = self.make_multiplicity(self.hosts_to_launch,
                    self.params["nrclients"])

    @staticmethod
    def make_multiplicity(names, multi):
        """ Produce a list in which each unique entry of the list "names"
        occurs "multi" times. The order of elements in names is preserved.
        
        >>> names = ['a', 'b', 'a', 'c', 'c', 'a', 'a']
        >>> StartClientsTask.make_multiplicity(names, 1)
        ['a', 'b', 'c']
        >>> StartClientsTask.make_multiplicity(names, 2)
        ['a', 'a', 'b', 'b', 'c', 'c']
        """
        result = []
        # Use OrderedDict to get unique names, preserving order
        for name in OrderedDict.fromkeys(names, None):
            result.extend([name] * multi)
        return result

    def is_alive(self, clientid):
        # Simplistic: just test if process with that pid exists and accepts
        # signals from us. TODO: better testing here, probably with ps|grep
        # or some such
        rc = self.kill_client(clientid, signal=0)
        return (rc == 0)
    
    def launch_clients(self):
        for host in self.hosts_to_launch:
            self.launch_one_client(host.strip())
        running_clients = [(cid, self.hosts[cid], pid) for (cid, pid) in
            self.pids.items()]
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
    # Client was started, but does not exist any more. Remove from state,
    #   then start and add again
    # Client was started, and does still exists. Nothing to do.
    
    def launch_one_client(self, host, clientid = None):
        if clientid is None:
            clientid = self.make_unique_id(host)
        # Check if client is already running
        if clientid in self.pids:
            assert self.hosts[clientid] == host
            if self.is_alive(clientid):
                self.logger.info("Client %s on host %s with PID %d already "
                                 "running",
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
            process = cadocommand.RemoteCommand(wuclient, host, self.parameters)
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
        if host == "localhost":
            process = cadocommand.Command(kill)
        else:
            process = cadocommand.RemoteCommand(kill, host, self.parameters)
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
    WANT_TO_RUN = object()
    SUBSCRIBE_WU_NOTIFICATIONS = object()

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

class CompleteFactorization(wudb.DbAccess, cadoparams.UseParameters,
                            DoesLogging, patterns.Mediator):
    """ The complete factorization, aggregate of the individual tasks """
    @property
    def name(self):
        return "tasks"
    @property
    def param_nodename(self):
        return self.name
    @property
    def title(self):
        return "Complete Factorization"
    
    def __init__(self, db, parameters, path_prefix):
        super().__init__(db = db, parameters = parameters, path_prefix = path_prefix)
        self.params = self.parameters.myparams(("name", "workdir"))
        self.db_listener = self.make_db_listener()
        self.registered_filenames = self.make_db_dict('server_registered_filenames')
        
        # Init WU BD
        self.wuar = self.make_wu_access()
        self.wuar.create_tables()
        
        # Set up WU server
        serverparams = parameters.myparams(["address", "port"], path_prefix + ["server"])
        serveraddress = serverparams.get("address", socket.gethostname())
        serverport = serverparams.get("port", 8001)
        uploaddir = self.params["workdir"].rstrip(os.sep) + os.sep + self.params["name"] + ".upload/"
        threaded = False
        self.server = wuserver.ServerLauncher(serveraddress, serverport, threaded, self.get_db_filename(),
            self.registered_filenames, uploaddir, bg=True, only_registered=True)
        
        # Init client lists
        self.clients = []
        for (path, key) in self.parameters.get_parameters().find(['slaves'], 'hostnames'):
            self.clients.append(StartClientsTask(serveraddress, serverport,
                                                 mediator = self,
                                                 db = db,
                                                 parameters = self.parameters,
                                                 path_prefix = path))
        
        parampath = self.parameters.get_param_path()
        sievepath = parampath + ['sieve']
        filterpath = parampath + ['filter']
        linalgpath = parampath + ['linalg']
        
        self.polysel = PolyselTask(mediator = self,
                                   db = db,
                                   parameters = self.parameters,
                                   path_prefix = parampath)
        self.fb = FactorBaseTask(mediator = self,
                                 db = db,
                                 parameters = self.parameters,
                                 path_prefix = sievepath)
        self.freerel = FreeRelTask(mediator = self,
                                   db = db,
                                   parameters = self.parameters,
                                   path_prefix = sievepath)
        self.sieving = SievingTask(mediator = self,
                                   db = db,
                                   parameters = self.parameters,
                                   path_prefix = sievepath)
        self.dup1 = Duplicates1Task(mediator = self,
                                    db = db,
                                    parameters = self.parameters,
                                    path_prefix = filterpath)
        self.dup2 = Duplicates2Task(mediator = self,
                                    db = db,
                                    parameters = self.parameters,
                                    path_prefix = filterpath)
        self.purge = PurgeTask(mediator = self,
                               db = db,
                               parameters = self.parameters,
                               path_prefix = filterpath)
        self.merge = MergeTask(mediator = self,
                               db = db,
                               parameters = self.parameters,
                               path_prefix = filterpath)
        self.linalg = LinAlgTask(mediator = self,
                                 db = db,
                                 parameters = self.parameters,
                                 path_prefix = linalgpath)
        self.characters = CharactersTask(mediator = self,
                                         db = db,
                                         parameters = self.parameters,
                                         path_prefix = linalgpath)
        self.sqrt = SqrtTask(mediator = self,
                             db = db,
                             parameters = self.parameters,
                             path_prefix = parampath)
        
        # Defines an order on tasks in which tasks that want to run should be
        # run
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
            Request.GET_DEPENDENCY_FILENAME: \
                self.linalg.get_dependency_filename,
            Request.GET_LINALG_PREFIX: self.linalg.get_prefix,
            Request.GET_KERNEL_FILENAME: self.characters.get_kernel_filename,
            Request.GET_WU_RESULT: self.db_listener.send_result,
        }
    
    def run(self):
        self.server.serve()
        
        try:
            for clients in self.clients:
                clients.launch_clients()
            
            while self.run_next_task():
                pass
            
            for task in self.tasks:
                task.print_stats()
            
        except KeyboardInterrupt:
            self.logger.fatal("Received KeyboardInterrupt. Terminating")
        
        for c in self.clients:
            c.kill_all_clients()
        
        self.server.shutdown()
    
    def run_next_task(self):
        for task in self.tasks:
            if task in self.tasks_that_want_to_run:
                # self.logger.info("Next task that wants to run: %s",
                #                  task.title)
                self.tasks_that_want_to_run.remove(task)
                return task.run()
        return False
    
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
        """ The relay for letting Tasks talk to us and each other """
        assert isinstance(message, Notification)
        sender = message.get_sender()
        key = message.get_key()
        value = message.get_value()
        self.logger.message("Received notification from %s, key = %s, value = %s",
                            sender, Notification.reverse_lookup(key), value)
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
                self.sieving.cancel_available_wus()
            else:
                raise Exception("Got HAVE_ENOUGH_RELATIONS from unknown sender")
        elif key is Notification.REGISTER_FILENAME:
            if isinstance(sender, ClientServerTask):
                self.register_filename(value)
            else:
                raise Exception("Got REGISTER_FILENAME, but not from a ClientServerTask")
        elif key is Notification.WANT_TO_RUN:
            if sender in self.tasks_that_want_to_run:
                raise Exception("Got request from %s to run, but it was in run queue already",
                                sender)
            else:
                self.tasks_that_want_to_run.append(sender)
        elif key is Notification.SUBSCRIBE_WU_NOTIFICATIONS:
            return self.db_listener.subscribeObserver(sender)
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
            result = self.request_map[key]()
        else:
            result = self.request_map[key](value)
        self.logger.message("Completed request from %s, key = %s, values = %s, result = %s",
                            sender, Request.reverse_lookup(key), value, result)
        return result
    
    def handle_message(self, message):
        if isinstance(message, Notification):
            self.relay_notification(Notification)
        elif isinstance(message, Request):
            return self.answer_request(message)
        else:
            raise TypeError("Message is neither Notification nor Request")
