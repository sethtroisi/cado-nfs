import os
import platform
import abc
import inspect
import hashlib
import logging
import cadocommand
import cadologger
from cadoparams import BoolParam


class InspectType(type):
    """ Meta-class that adds a class attribute "init_signature" with the
    signature of the __init__() method, as produced by
    inspect.getfullargspec()
    """
    def __init__(cls, name, bases, dct):
        super().__init__(name, bases, dct)
        # inspect.getfullargspec() produces an object with attributes:
        # args, varargs, kwonlyargs, annotations (among others)
        # where args is a list that contains the names of positional parameters
        # (including those with default values), varargs contains the name of
        # the list of the variable-length positional parameters (i.e., the name
        # of the * catch-all), and kwonlyargs contains a list of keyword-only
        # parameters.
        # E.g.: def f(a:"A", b:"B"=1, *args:"*ARGS", c:"C", d:"D"=3,
        #             **kwargs:"**KWARGS"):
        # produces an object with attributes args=['a', 'b'], varargs='args',
        # varkw='kwargs', defaults=(1,), kwonlyargs=['c', 'd'],
        # kwonlydefaults={'d': 3}, annotations={'a': 'A', 'c': 'C', 'b': 'B',
        # 'd': 'D', 'args': '*ARGS', 'kwargs': '**KWARGS'}
        cls.init_signature = inspect.getfullargspec(cls.__init__)


class Option(object, metaclass=abc.ABCMeta):
    ''' Base class for command line options that may or may not take parameters
    '''
    # TODO: decide whether we need '-' or '/' as command line option prefix
    # character under Windows. Currently: use "-'
    if False and platform.system() == "Windows":
        prefix = '/'
    else:
        prefix = '-'

    def __init__(self, arg=None, prefix=None, is_input_file=False,
                 is_output_file = False, checktype = None):
        """ Define a mapping from a parameter name and value to command line
        parameters.

        arg is the command line parameter to which this should map, e.g.,
            'v' or 'thr'. If not specified, the default name supplied to the
            map() method will be used.
        prefix is the command line parameter prefix to use; the default is the
            class variable which currently is '-'. Some programs, e.g. bwc.pl,
            want some parameters with a different prefix, e.g., ':complete'
        is_input_file must be set to True if this parameter gives the filename
            of an input file to the command. This will be used to generate FILE
            lines in workunits.
        is_output_file must be set to True if this parameter gives the filename
            of an output file of the command. This will be used to generate
            RESULT lines in workunits.
        if checktype is specified, map() will assert that the value is an
            instance of checktype.
        """
        self.arg = arg
        if prefix:
            self.prefix = prefix # hides the class variable
        self.is_input_file = is_input_file
        self.is_output_file = is_output_file
        self.checktype = checktype
        self.defaultname = None

    def set_defaultname(self, defaultname):
        self.defaultname = defaultname

    def get_arg(self):
        assert not self.defaultname is None
        return self.defaultname if self.arg is None else self.arg

    def map(self, value):
        """ Public class that converts the Option instance into an array of
        strings with command line parameters. It also checks the type, if
        checktype was specified in the constructor.
        """
        if not self.checktype is None:
            assert isinstance(value, self.checktype)
        return self._map(value)

    @abc.abstractmethod
    def _map(self, value):
        """ Private class that does the actual translation to an array of string
        with command line parameters. Subclasses need to implement it.

        A simple case of the Template Method pattern.
        """
        pass

class PositionalParameter(Option):
    ''' Positional command line parameter '''
    def _map(self, value):
        return [str(value)]

class Parameter(Option):
    ''' Command line option that takes a parameter '''
    def _map(self, value):
        return [self.prefix + self.get_arg(), str(value)]

class ParameterEq(Option):
    ''' Command line option that takes a parameter, used on command line in
    the form of "key=value" '''
    def _map(self, value):
        return ["%s=%s" % (self.get_arg(), value)]

class Toggle(Option):
    ''' Command line option that does not take a parameter.
    value is interpreted as a truth value, the option is either added or not
    '''
    def _map(self, value):
        if isinstance(value, BoolParam):
            value = bool(value)
        if value is True:
            return [self.prefix + self.get_arg()]
        elif value is False:
            return []
        else:
            raise ValueError("Toggle.map() for %s requires a boolean "
                             "type argument" % self.get_arg())

class Sha1Cache(object):
    """ A class that computes SHA1 sums for files and caches them, so that a
    later request for the SHA1 sum for the same file is not computed again.
    File identity is checked via the file's realpath and the file's inode, size,
    and modification time.
    """
    def __init__(self):
        self._sha1 = {}

    @staticmethod
    def _read_file_in_blocks(file_object):
        blocksize = 65536
        while True:
            data = file_object.read(blocksize)
            if not data:
                break
            yield data

    def get_sha1(self, filename):
        realpath = os.path.realpath(filename)
        stat = os.stat(realpath)
        file_id = (stat.st_ino, stat.st_size, stat.st_mtime)
        # Check whether the file on disk changed
        if realpath in self._sha1 and not self._sha1[realpath][1] == file_id:
            logger = logging.getLogger("Sha1Cache")
            logger.warn("File %s changed! Discarding old SHA1 sum", realpath)
            del(self._sha1[realpath])
        if not realpath in self._sha1:
            logger = logging.getLogger("Sha1Cache")
            logger.debug("Computing SHA1 for file %s", realpath)
            sha1 = hashlib.sha1()
            with open(realpath, "rb") as inputfile:
                for data in self._read_file_in_blocks(inputfile):
                    sha1.update(data)
            self._sha1[realpath] = (sha1.hexdigest(), file_id)
            logger.debug("SHA1 for file %s is %s", realpath,
                         self._sha1[realpath])
        return self._sha1[realpath][0]

sha1cache = Sha1Cache()

class Program(object, metaclass=InspectType):
    ''' Base class that represents programs of the CADO suite

    The Program base class is oblivious to how programs get input data,
    or how they provide output data. It does, however, handle redirecting
    stdin/stdout/stderr from/to files and accepts file names to/from which
    to redirect. This is done so that workunits can be generated from
    Program instances; in the workunit text, shell redirection syntax needs
    to be used to connect stdio to files, so a Program instance needs to be
    aware of which file names should be used for stdio redirection.

    Sub-classes correspond to individual programs (such as las) and must
    define these class variables:
    binary: a string with the name of the binary executable file
    name: a name used internally in the Python scripts for this program, must
      start with a letter and contain only letters and digits; if the binary
      file name is of this form, then that can be used
    params: a mapping that tells how to translate configuration parameters to
      command line options. The keys of the mapping are the keys as in the
      configuration file, e.g., "verbose" in a configuartion file line
      tasks.polyselect.verbose = 1
      The values of the mapping are instances of subclasses of Option,
      initialised with the command line parameter they should map to, e.g.,
      "verbose" should map to an instance Toggle("v"), and "threads" should
      map to an instance Parameter("t")

    >>> p = Ls()
    >>> p.make_command_line()
    '/bin/ls'
    >>> p = Ls(stdout='foo')
    >>> p.make_command_line()
    '/bin/ls > foo'
    >>> p = Ls(stderr='foo')
    >>> p.make_command_line()
    '/bin/ls 2> foo'
    >>> p = Ls(stdout='foo', stderr='bar')
    >>> p.make_command_line()
    '/bin/ls > foo 2> bar'
    >>> p = Ls(stdout='foo', stderr='foo')
    >>> p.make_command_line()
    '/bin/ls > foo 2>&1'
    >>> p = Ls(stdout='foo', append_stdout=True)
    >>> p.make_command_line()
    '/bin/ls >> foo'
    >>> p = Ls(stderr='foo', append_stderr=True)
    >>> p.make_command_line()
    '/bin/ls 2>> foo'
    >>> p = Ls(stdout='foo', append_stdout=True, stderr='bar', append_stderr=True)
    >>> p.make_command_line()
    '/bin/ls >> foo 2>> bar'
    >>> p = Ls(stdout='foo', append_stdout=True, stderr='foo', append_stderr=True)
    >>> p.make_command_line()
    '/bin/ls >> foo 2>&1'
    >>> p = Ls('foo', 'bar')
    >>> p.make_command_line()
    '/bin/ls foo bar'
    >>> p = Ls('foo', 'bar', long = True)
    >>> p.make_command_line()
    '/bin/ls -l foo bar'
    '''

    path = '.'
    subdir = ""
    paramnames = ("execpath", "execsubdir", "execbin", "runprefix")
    # These should be abstract properties, but we want to reference them as
    # class attributes, which properties can't. Ergo dummy variables
    binary = None

    # This class variable definition should not be here. It gets overwritten
    # when the InspectType meta-class creates the class object. The only purpose
    # is to make pylint shut up about the class not having an init_signature
    # attribute
    init_signature = None

    def __init__(self, options, stdin = None,
                 stdout = None, append_stdout = False, stderr = None,
                 append_stderr = False, background = False, execpath = None,
                 execsubdir = None, execbin = None, runprefix=None):
        ''' Takes a dict of of command line options. Defaults are filled in
        from the cadoaprams.Parameters instance parameters.

        The stdin, stdout, and stderr parameters accept strings. If a string is
        given, it is interpreted as the file name to use for redirecting that
        stdio stream. In workunit generation, the file name is used for shell
        redirection.
        '''

        # The function annotations define how to convert each value in
        # self.parameters to command line arguments, but the annotations are
        # stored in an unordered dict. We would like to put the command line
        # parameters in a particular order, so that, for example, positional
        # parameters are placed correctly. We also need to handle a variable-
        # length positional parameter list.
        # That is, we need three lists:
        # - Individual positional parameters
        # - Variable-length positional parameter (e.g., list of filenames)
        # - Keyword parameters

        # First off, we set the default names for Option subclass instances.
        # Since the Option instances are always the same (the annotation values
        # are evaluated once at class creation time), we'll set the defaultnames
        # on the same Option instances again each time a Program subclass gets
        # instantiated, but that does not hurt.
        for (name, option) in self.__init__.__annotations__.items():
            option.set_defaultname(name)

        # "options" contains all local symbols defined in the sub-class'
        # __init__() method. This includes self, __class__, and a few others.
        # Filter them so that only those that correspond to annotated
        # parameters remain.
        filtered = self._filter_annotated_keys(options)

        # The default value for all __init__() parameters that correspond to
        # command line parameters is always None. Remove those entries.
        # We could in principle look up the default values from the .defaults
        # attribute given by inspect.getfullargspec(), but it's easier to
        # use the convention that the default is always None.
        filtered = {key:options[key] for key in filtered \
                    if not options[key] is None}

        self.parameters = filtered

        self.stdin = stdin
        self.stdout = stdout
        self.append_stdout = append_stdout
        self.stderr = stderr
        self.append_stderr = append_stderr
        self.runprefix = runprefix

        # If we are to run in background, we add " &" to the command line.
        # We require that stdout and stderr are redirected to files.
        self.background = background
        if self.background and (stdout is None or stderr is None):
            raise Exception("Programs to run in background must redirect "
                            "stdout and stderr")

        # Look for location of the binary executable at __init__, to avoid
        # calling os.path.isfile() multiple times
        path = str(execpath or self.path)
        subdir = str(execsubdir or self.subdir)
        binary = str(execbin or self.binary)
        execfile = os.path.normpath(os.sep.join([path, binary]))
        execsubfile = os.path.normpath(os.sep.join([path, subdir, binary]))
        if execsubfile != execfile and os.path.isfile(execsubfile):
            self.execfile = execsubfile
        elif os.path.isfile(execfile):
            self.execfile = execfile
        else:
            self.execfile = binary
            # Raising an exception here makes it impossible to run doctests
            # for each Program subclass
            # raise Exception("Binary executable file %s not found" % execfile)

    @classmethod
    def _get_option_annotations(cls):
        """ Extract the elements of this class' __init__() annotations where
        the annotation is an Option instance
        """
        return {key:val for (key, val) in cls.__init__.__annotations__.items()
                if isinstance(val, Option)}

    @classmethod
    def _filter_annotated_keys(cls, keys):
        """ From the list of keys given in "keys", return those that are
        parameters of the __init__() method of this class and annotated with an
        Option instance.
        """
        return [key for key in keys if key in cls._get_option_annotations()]

    @classmethod
    def get_accepted_keys(cls):
        """ Return all parameter file keys which can be used by this program,
        including those that don't directly map to command line parameters,
        like those specifying search paths for the binary executable file.

        This list is the union of annotated positional and keyword parameters to
        the constructor, plus the Program class defined parameters. Note that an
        *args parameter (stored in cls.init_signature.varargs) should not be
        included here; there is no way to specify a variable-length set of
        positional parameters in the parameter file.
        """
        # The fact that we want to exclude varargs is why we need the inspect
        # info here.
        parameters = cls._filter_annotated_keys(cls.init_signature.args +
                                                cls.init_signature.kwonlyargs)
        return parameters + list(Program.paramnames)

    def get_stdout(self):
        return self.stdout

    def get_stderr(self):
        return self.stderr

    def get_stdio(self):
        """ Returns a 3-tuple with information on the stdin, stdout, and stderr
        streams for this program.

        The first entry is for stdin, and is either a file name or None.
        The second entry is for stdout and is a 2-tuple, whose first entry is
            either a file name or None, and whose second entry is True if we
            should append to the file and False if we should not.
        The third entry is like the second, but for stderr.
        """
        return (self.stdin,
                (self.stdout, self.append_stdout),
                (self.stderr, self.append_stderr)
            )

    def _get_files(self, is_output):
        """ Helper method to get a list of input/output files for this Program.
        This method returns the list of filenames without considering files for
        stdio redirection. If "is_output" evaluates to True, the list of output
        files is generated, otherwise the list of input files.
        """
        files = []
        parameters = self.init_signature.args + self.init_signature.kwonlyargs
        if not self.init_signature.varargs is None:
            parameters.append(self.init_signature.varargs)
        parameters = self._filter_annotated_keys(parameters)
        for param in parameters:
            ann = self.__init__.__annotations__[param]
            if param in self.parameters and \
                    (ann.is_output_file if is_output else ann.is_input_file):
                if param == self.init_signature.varargs:
                    # vararg is a list; we need to convert each entry to string
                    # and append this list of strings to files
                    files.extend(map(str, self.parameters[param]))
                else:
                    files.append(str(self.parameters[param]))
        return files

    def get_input_files(self):
        input_files = self._get_files(is_output = False)
        if isinstance(self.stdin, str):
            input_files.append(self.stdin)
        return input_files

    def get_regular_output_files(self):
        """ Returns a list of output files, excluding files for stdout/stderr
        redirection.
        """
        return self._get_files(is_output = True)

    def get_output_files(self):
        """ Returns a list of output files, including files for stdout/stderr
        redirection if such redirection is used
        """
        output_files = self.get_regular_output_files()
        for filename in (self.stdout, self.stderr):
            if isinstance(filename, str):
                output_files.append(filename)
        return output_files

    def get_exec_file(self):
        return self.execfile

    def get_exec_files(self):
        return [self.get_exec_file()]

    @staticmethod
    def translate_path(filename, path = None):
        if not path:
            return filename
        else:
            return str(path).rstrip(os.sep) + os.sep + \
                    os.path.basename(str(filename))

    def make_command_array(self, binpath = None, inputpath = None,
                           outputpath = None):
        # Begin command line with program to execute
        command = []
        if not self.runprefix is None:
            command.append(self.runprefix)
        command.append(self.translate_path(self.get_exec_file(), binpath))

        # Add keyword command line parameters, then positional parameters
        parameters = self.init_signature.kwonlyargs + self.init_signature.args
        parameters = self._filter_annotated_keys(parameters)
        for key in parameters:
            if key in self.parameters:
                ann = self.__init__.__annotations__[key]
                value = self.parameters[key]
                # If this is an input or an output file name, we may have to
                # translate it, e.g., for workunits
                assert not (ann.is_input_file and ann.is_output_file)
                if ann.is_input_file:
                    value = self.translate_path(value, inputpath)
                elif ann.is_output_file:
                    value = self.translate_path(value, outputpath)
                command += ann.map(value)

        # Add positional command line parameters
        key = self.init_signature.varargs
        if not key is None:
            paramlist = self.parameters[key]
            ann = self.__init__.__annotations__.get(key, None)
            for value in paramlist:
                command += ann.map(value)
        return command

    def make_command_line(self, binpath=None, inputpath=None, outputpath=None):
        """ Make a shell command line for this program.

        If files are given for stdio redirection, the corresponding redirection
        tokens are added to the command line.
        """
        cmdarr = self.make_command_array(binpath, inputpath, outputpath)
        cmdline = " ".join(map(cadocommand.shellquote, cmdarr))
        if isinstance(self.stdin, str):
            translated = self.translate_path(self.stdin, inputpath)
            cmdline += ' < ' + cadocommand.shellquote(translated)
        if isinstance(self.stdout, str):
            redir = ' >> ' if self.append_stdout else ' > '
            translated = self.translate_path(self.stdout, outputpath)
            cmdline += redir + cadocommand.shellquote(translated)
        if not self.stderr is None and self.stderr is self.stdout:
            cmdline += ' 2>&1'
        elif isinstance(self.stderr, str):
            redir = ' 2>> ' if self.append_stderr else ' 2> '
            translated = self.translate_path(self.stderr, outputpath)
            cmdline += redir + cadocommand.shellquote(translated)
        if self.background:
            cmdline += " &"
        return cmdline

    def make_wu(self, wuname):
        def append_file(wu, key, filename):
            wu.append('%s %s' % (key, os.path.basename(filename)))
            wu.append('CHECKSUM %s' % sha1cache.get_sha1(filename))
        workunit = ['WORKUNIT %s' % wuname]
        for filename in self.get_input_files():
            append_file(workunit, 'FILE', filename)
        append_file(workunit, 'EXECFILE', self.get_exec_file())
        cmdline = self.make_command_line(binpath = "${EXECDIR}",
            inputpath = "${DLDIR}", outputpath = "${WORKDIR}")
        workunit.append('COMMAND %s' % cmdline)
        for filename in self.get_output_files():
            workunit.append('RESULT %s' % os.path.basename(filename))
        workunit.append("") # Make a trailing newline
        return '\n'.join(workunit)


class Polyselect2l(Program):
    """
    >>> p = Polyselect2l(P=5, N=42, degree=4, verbose=True)
    >>> p.make_command_line()
    'polyselect2l -N 42 -degree 4 -v 5'
    >>> p = Polyselect2l(P=5, N=42, degree=4, verbose=True)
    >>> p.make_command_line()
    'polyselect2l -N 42 -degree 4 -v 5'
    """
    binary = "polyselect2l"
    name = binary
    subdir = "polyselect"

    def __init__(self, *,
                 P : Parameter(), 
                 N : Parameter(),
                 degree : Parameter(),
                 verbose : Toggle("v") = None,
                 quiet : Toggle("q") = None,
                 sizeonly : Toggle("r") = None,
                 threads : Parameter("t") = None,
                 admin : Parameter() = None,
                 admax : Parameter() = None,
                 incr : Parameter() = None,
                 nq : Parameter() = None,
                 save : Parameter(is_output_file = True) = None,
                 resume : Parameter(is_input_file = True) = None,
                 maxtime : Parameter() = None,
                 out : Parameter() = None,
                 printdelay : Parameter("s") = None,
                 keep: Parameter() = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)


class MakeFB(Program):
    """
    >>> p = MakeFB("foo.poly")
    >>> p.make_command_line()
    'makefb -poly foo.poly'
    >>> p = MakeFB(poly="foo.poly", nopowers=True, maxbits=5, stdout="foo.roots")
    >>> p.make_command_line()
    'makefb -poly foo.poly -nopowers -maxbits 5 > foo.roots'
    """
    binary = "makefb"
    name = binary
    subdir = "sieve"

    def __init__(self, *,
                 poly: Parameter(is_input_file = True),
                 alim: Parameter(),
                 maxbits: Parameter() = None,
                 out: Parameter() = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)


class FreeRel(Program):
    """
    >>> p = FreeRel("foo.poly", "foo.renumber")
    >>> p.make_command_line()
    'freerel -poly foo.poly -renumber foo.renumber'
    >>> p = FreeRel(poly="foo.poly", renumber="foo.renumber", badideals="foo.bad", verbose=True, pmin=123, pmax=234)
    >>> p.make_command_line()
    'freerel -poly foo.poly -renumber foo.renumber -badideals foo.bad -v -pmin 123 -pmax 234'
    """
    binary = "freerel"
    name = binary
    subdir = "sieve"
    def __init__(self, *,
                 poly: Parameter(is_input_file = True),
                 renumber: Parameter(is_output_file = True),
                 lpbr: Parameter(),
                 lpba: Parameter(),
                 out: Parameter(is_output_file = True),
                 badideals: Parameter(is_output_file = True) = None,
                 pmin: Parameter() = None,
                 pmax: Parameter() = None,
                 dlp: Toggle("addfullcol") = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Las(Program):
    binary = "las"
    name = binary
    subdir = "sieve"
    def __init__(self,
                 I: Parameter(),
                 poly: Parameter(is_input_file = True),
                 factorbase: Parameter("fb", is_input_file = True),
                 q0: Parameter(),
                 q1: Parameter() = None,
                 rho: Parameter() = None,
                 tdthresh: Parameter() = None,
                 bkthresh: Parameter() = None,
                 rlim: Parameter() = None,
                 alim: Parameter() = None,
                 lpbr: Parameter() = None,
                 lpba: Parameter() = None,
                 mfbr: Parameter() = None,
                 mfba: Parameter() = None,
                 rlambda: Parameter() = None,
                 alambda: Parameter() = None,
                 skewness: Parameter("S") = None,
                 verbose: Toggle("v") = None,
                 out: Parameter(is_output_file = True) = None,
                 threads: Parameter("mt") = None,
                 ratq: Toggle() = None,
                 stats_stderr: Toggle("stats-stderr") = None,
                 # Let's make fbcache neither input nor output file. It should
                 # not be distributed to clients, nor sent back to the server.
                 # It's a local temp file, but re-used between different runs.
                 fbcache: Parameter("fbc") = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)


class Duplicates1(Program):
    binary = "dup1"
    name = binary
    subdir = "filter"
    def __init__(self,
                 *args: PositionalParameter(is_input_file = True),
                 prefix : Parameter(),
                 out: Parameter() = None,
                 outfmt: Parameter() = None,
                 bzip: Toggle("bz") = None,
                 only: Parameter() = None,
                 nslices_log: Parameter("n") = None,
                 lognrels: Parameter() = None,
                 filelist: Parameter(is_input_file = True) = None,
                 basepath: Parameter() = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)


class Duplicates2(Program):
    binary = "dup2"
    name = binary
    subdir = "filter"
    def __init__(self,
                 *args: PositionalParameter(is_input_file = True),
                 rel_count: Parameter("nrels"),
                 poly: Parameter(is_input_file = True),
                 renumber: Parameter(is_input_file = True),
                 filelist: Parameter(is_input_file = True) = None,
                 badidealinfo: Parameter(is_input_file = True) = None,
                 dlp: Toggle("dl") = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Purge(Program):
    binary = "purge"
    name = binary
    subdir = "filter"
    def __init__(self,
                 *args: PositionalParameter(is_input_file = True),
                 out: Parameter(is_output_file = True),
                 filelist: Parameter(is_input_file = True) = None,
                 basepath: Parameter() = None,
                 subdirlist: Parameter() = None,
                 nrels: Parameter() = None,
                 outdel: Parameter(is_output_file = True) = None,
                 sos: Parameter(is_output_file = True) = None,
                 keep: Parameter() = None,
                 minindex: Parameter() = None,
                 nprimes: Parameter() = None,
                 raw: Toggle() = None,
                 threads: Parameter("npthr") = None,
                 inprel: Parameter(is_input_file = True) = None,
                 outrel: Parameter(is_output_file = True) = None,
                 npass: Parameter() = None,
                 required_excess: Parameter() = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Merge(Program):
    binary = "merge"
    name = binary
    subdir = "filter"
    def __init__(self,
                 mat: Parameter(is_input_file = True),
                 out: Parameter(is_output_file = True),
                 maxlevel: Parameter() = None,
                 keep: Parameter() = None,
                 skip: Parameter() = None,
                 forbw: Parameter() = None,
                 ratio: Parameter() = None,
                 coverNmax: Parameter() = None,
                 nbmergemax: Parameter() = None,
                 resume: Parameter() = None,
                 mkztype: Parameter() = None,
                 wmstmax: Parameter() = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class MergeDLP(Program):
    binary = "merge-dl"
    name = binary
    subdir = "filter"
    def __init__(self,
                 mat: Parameter(is_input_file = True),
                 out: Parameter(is_output_file = True),
                 maxlevel: Parameter() = None,
                 keep: Parameter() = None,
                 skip: Parameter() = None,
                 forbw: Parameter() = None,
                 ratio: Parameter() = None,
                 coverNmax: Parameter() = None,
                 nbmergemax: Parameter() = None,
                 resume: Parameter() = None,
                 mkztype: Parameter() = None,
                 wmstmax: Parameter() = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

# Todo: define is_input_file/is_output_file for remaining programs
class Replay(Program):
    binary = "replay"
    name = binary
    subdir = "filter"
    def __init__(self,
                 skip: Parameter() = None,
                 purged: Parameter() = None,
                 history: Parameter("his") = None,
                 index: Parameter() = None,
                 out: Parameter() = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class ReplayDLP(Program):
    binary = "replay-dl"
    name = binary
    subdir = "filter"
    def __init__(self,
                 purged: Parameter() = None,
                 ideals: Parameter() = None,
                 history: Parameter("his") = None,
                 index: Parameter() = None,
                 out: Parameter() = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class MagmaLinalg(Program):
    binary = "magma-linalg-wrapper.sh"
    name = binary
    subdir = "scripts"
    def __init__(self,
                 gorder: Parameter("ell"),
                 nmaps: Parameter(),
                 sparsemat: Parameter(),
                 sm: Parameter(),
                 ker: Parameter() = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class BWC(Program):
    binary = "bwc.pl"
    name = "bwc"
    subdir = "linalg/bwc"
    def __init__(self,
                 complete: Toggle(prefix=":") = None,
                 wipeout: Toggle(prefix=":") = None,
                 dryrun: Toggle("d") = None,
                 verbose: Toggle("v") = None,
                 mpi: ParameterEq() = None,
                 threads: ParameterEq("thr") = None,
                 mn: ParameterEq() = None,
                 nullspace: ParameterEq() = None,
                 interval: ParameterEq() = None,
                 ys: ParameterEq() = None,
                 matrix: ParameterEq() = None,
                 wdir: ParameterEq() = None,
                 mpiexec: ParameterEq() = None,
                 hosts: ParameterEq() = None,
                 hostfile: ParameterEq() = None,
                 interleaving: ParameterEq() = None,
                 shuffled_product: ParameterEq() = None,
                 bwc_bindir: ParameterEq() = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class SM(Program):
    binary = "sm"
    name = binary
    subdir = "filter"
    def __init__(self, *,
                 poly: Parameter(),
                 purged: Parameter(),
                 index: Parameter(),
                 out: Parameter(),
                 gorder: Parameter(),
                 smexp: Parameter(),
                 threads: Parameter("mt") = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)
 
class Characters(Program):
    binary = "characters"
    name = binary
    subdir = "linalg"
    def __init__(self, *,
                 poly: Parameter(),
                 purged: Parameter(),
                 index: Parameter(),
                 heavyblock: Parameter(),
                 out: Parameter(),
                 wfile: Parameter("ker"),
                 lpbr: Parameter(),
                 lpba: Parameter(),
                 nchar: Parameter() = None,
                 threads: Parameter("t") = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Sqrt(Program):
    binary = "sqrt"
    name = binary
    subdir = "sqrt"
    def __init__(self, *,
                 poly: Parameter(),
                 prefix: Parameter(),
                 purged: Parameter(),
                 index: Parameter(),
                 kernel: Parameter("ker"),
                 dep: Parameter() = None,
                 ab: Toggle() = None,
                 rat: Toggle() = None,
                 alg: Toggle() = None,
                 gcd: Toggle() = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class WuClient(Program):
    binary = "wuclient2.py"
    name = "wuclient"
    subdir = "scripts/cadofactor"
    def __init__(self,
                 server: Parameter(prefix='--'),
                 daemon: Toggle(prefix='--') = None,
                 keepoldresult: Toggle(prefix='--') = None,
                 nosha1check: Toggle(prefix='--') = None,
                 dldir: Parameter(prefix='--') = None,
                 workdir: Parameter(prefix='--') = None,
                 bindir: Parameter(prefix='--') = None,
                 clientid: Parameter(prefix='--') = None,
                 basepath: Parameter(prefix='--') = None,
                 getwupath: Parameter(prefix='--') = None,
                 loglevel: Parameter(prefix='--') = None,
                 postresultpath: Parameter(prefix='--') = None,
                 downloadretry: Parameter(prefix='--') = None,
                 logfile: Parameter(prefix='--') = None,
                 debug: Parameter(prefix='--') = None,
                 niceness: Parameter(prefix='--') = None,
                 wu_filename: Parameter(prefix='--') = None,
                 arch: Parameter(prefix='--') = None,
                 certsha1: Parameter(prefix='--') = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class SSH(Program):
    binary = "ssh"
    name = binary
    path = "/usr/bin"
    def __init__(self,
                 host: PositionalParameter(),
                 *args: PositionalParameter(),
                 compression: Toggle("C") = None,
                 verbose: Toggle("v") = None,
                 cipher: Parameter("c") = None,
                 configfile: Parameter("F") = None,
                 identity_file: Parameter("i") = None,
                 login_name: Parameter("l") = None,
                 port: Parameter("p") = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class RSync(Program):
    binary = "rsync"
    name = binary
    path = "/usr/bin"
    def __init__(self,
                 sourcefile: PositionalParameter(),
                 remotefile: PositionalParameter(),
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Ls(Program):
    binary = "ls"
    name = binary
    path = "/bin"
    def __init__(self,
                 *args : PositionalParameter(),
                 long : Toggle('l') = None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Kill(Program):
    binary = "kill"
    name = binary
    path = "/bin"
    def __init__(self,
                 *args: PositionalParameter(),
                 signal: Parameter("s"),
                 **kwargs):
        super().__init__(locals(), **kwargs)


if __name__ == '__main__':
    PROGRAM = Ls('foo', 'bar')
    CMDLINE = PROGRAM.make_command_line()
    print(CMDLINE)
