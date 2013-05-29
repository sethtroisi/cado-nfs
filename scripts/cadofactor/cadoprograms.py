import os
import sys
import platform
import subprocess
import abc
import logging
import cadocommand

class Option(metaclass=abc.ABCMeta):
    ''' Base class for command line options that may or may not take parameters
    '''
    # TODO: decide whether we need '-' or '/' as command line option prefix 
    # character under Windows. Currently: use "-'
    if False and platform.system() == "Windows":
        prefix = '/'
    else:
        prefix = '-'
    
    def __init__(self, config, arg = None, prefix  = None, is_input_file = False, is_output_file = False):
        """ Define a mapping from a parameter dictionary to command line parameters.
        
        config is the key name of the parameter in the parameter dictionary, e.g., 'verbose', or 'threads'
        arg is the command line parameter to which this should map, e.g., 'v' or 'thr'. If not specified,
            it defaults to the same string as config.
        prefix is the command line parameter prefix to use; the default is the class variable which currently is '-'.
            Some programs, e.g. bwc.pl, want some parameters with a different prefix, e.g., ':complete'
        is_input_file must be set to True if this parameter gives the filename of an input file to the command.
            This will be used to generate FILE lines in workunits.
        is_output_file must be set to True if this parameter gives the filename of an output file of the command.
            This will be used to generate RESULT lines in workunits.
        """
        self.config = config
        self.is_input_file = is_input_file
        self.is_output_file = is_output_file
        if arg is None:
            self.arg = config
        else:
            self.arg = arg
        if prefix:
            self.prefix = prefix # hides the class variable
    
    def get_key(self):
        return self.config
    
    @abc.abstractmethod
    def  map(self, value):
        pass


class PositionalParameter(Option):
    ''' Positional command line parameter '''
    def map(self, value):
        return [value]

class Parameter(Option):
    ''' Command line option that takes a parameter '''
    def map(self, value):
        return [self.prefix + self.arg, value]

class ParameterEq(Option):
    ''' Command line option that takes a parameter, used on command line in 
    the form of "key=value" '''
    def map(self, value):
        return [self.arg + "=" + value]

class Toggle(Option):
    ''' Command line option that does not take a parameter. 
    value is interpreted as a truth value, the option is either added or not 
    '''
    def map(self, value):
        if value.lower() in ["yes", "true", "on", "1"]:
            return [self.prefix + self.arg]
        elif value.lower() in ["no", "false", "off", "0"]:
            return []
        else:
            raise ValueError("Toggle.map() requires a boolean type argument")


class Program(object):
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
    params: a mapping that tell how to translate configuration parameters to 
      command line options. The keys of the mapping are the keys as in the 
      configuration file, e.g., "verbose" in a configuartion file line
      tasks.polyselect.verbose = 1
      The values of the mapping are instances of subclasses of Option, 
      initialised with the command line parameter they should map to, e.g.,
      "verbose" should map to an instance Toggle("-v"), and "threads" should 
      map to an instance Parameter("-t")
    '''
    
    params_list = ("execpath", "execsubdir", "execbin")
    
    @staticmethod
    def __shellquote(s):
        ''' Quote a command line argument
        
        Currently does it the hard way: always encloses the argument in single
        quotes, and escapes any single quotes that are part of the argument
        '''
        return "'" + s.replace("'", "'\\''") + "'"
    
    def __init__(self, args, parameters, 
                 stdin = None, stdout = subprocess.PIPE, 
                 stderr = subprocess.PIPE):
        ''' Takes a list of positional parameters and a dictionary of command 
        line parameters 
        
        The stdin, stdout, and stderr parameters accept the same parameters 
        as subprocess.popen(), but also accept strings. If a string is given,
        it is interpreted as the file name to use for redirecting that stdio
        stream. In direct execution, the file is opened for reading or writing,
        resp., and the file handle is passed to popen(). In workunit 
        generation, the file name is used for shell redirection.
        '''
        
        self.args = args
        self.parameters = parameters
        
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr
        
        # Look for location of the binary executable at __init__, to avoid 
        # calling os.path.isfile() multiple times
        path = self.parameters.get("execpath", ".")
        subdir = self.parameters.get("execsubdir", self.subdir)
        binary = self.parameters.get("execbin", self.binary)
        execfile = os.path.normpath(os.sep.join([path, binary]))
        execsubfile = os.path.normpath(os.sep.join([path, subdir, binary]))
        if execsubfile != execfile and os.path.isfile(execsubfile):
            # print ("Found %s in %s" % (binary, execsubfile))
            self.execfile = execsubfile
        else:
            self.execfile = execfile
    
    @classmethod
    def get_param_keys(cls):
        """ Return the config file keys which map to command line arguments for this program """
        return [opt.get_key() for opt in cls.params_list]
    
    @classmethod
    def get_config_keys(cls):
        """ Return all config file keys which can be used by this program, including those that don't
        directly map to command line parameters, like those specifying search paths.
        """
        l = list(Program.params_list) + cls.get_param_keys()
        return l
    
    def get_input_files(self):
        input_files = []
        if isinstance(self.stdin, str):
            input_files.append(self.stdin)
        for p in self.params_list:
            key = p.get_key()
            if key in self.parameters:
                if p.is_input_file:
                    input_files.append(self.parameters[key])
        return input_files
    
    def get_output_files(self):
        output_files = []
        if isinstance(self.stdout, str):
            output_files.append(self.stdout)
        if isinstance(self.stderr, str):
            output_files.append(self.stderr)
        for p in self.params_list:
            key = p.get_key()
            if key in self.parameters:
                if p.is_output_file:
                    output_files.append(self.parameters[key])
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
            return path.rstrip(os.sep) + os.sep + os.path.basename(filename)
    
    def make_command_array(self, binpath = None, inputpath = None, outputpath = None):
        # Begin command line with program to execute
        command = [self.translate_path(self.get_exec_file(), binpath)]
        
        # Add keyword command line parameters
        for p in self.params_list:
            key = p.get_key()
            if key in self.parameters:
                assert not (p.is_input_file and p.is_output_file)
                argument = self.parameters[key]
                # If this is an input or an output file name, we may have to
                # translate it, e.g., for workunits
                if p.is_input_file:
                    argument = self.translate_path(self.parameters[key], inputpath)
                elif p.is_output_file:
                    argument = self.translate_path(self.parameters[key], outputpath)
                command += p.map(argument)
        
        # Add positional command line parameters
        # FIXME: how to identify input/output files here?
        if self.args:
            command += self.args
        return command
    
    def make_command_line(self, binpath = None, inputpath = None, outputpath = None):
        """ Make a shell command line for this program.
        
        If files are given for stdio redirection, the corresponding redirection tokens
        are added to the command line.
        """
        cmdarr = self.make_command_array(binpath, inputpath, outputpath)
        cmdline = " ".join([Program.__shellquote(arg) for arg in cmdarr])
        if isinstance(self.stdin, str):
            cmdline += ' < ' + Program.__shellquote(self.translate_path(self.stdin, inputpath))
        if isinstance(self.stdout, str):
            cmdline += ' > ' + Program.__shellquote(self.translate_path(self.stdout, outputpath))
        if isinstance(self.stderr, str):
            cmdline += ' 2> ' + Program.__shellquote(self.translate_path(self.stdout, outputpath))
        if self.stderr is subprocess.STDOUT:
            cmdline += ' 2>&1'
        return cmdline
    
    def make_wu(self, wuname):
        wu = ['WORKUNIT %s' % wuname]
        for f in self.get_input_files():
            wu.append('FILE %s' % os.path.basename(f))
        wu.append('EXECFILE %s' % os.path.basename(self.get_exec_file()))
        wu.append('COMMAND %s' % self.make_command_line(binpath = "${DLDIR}", inputpath = "${DLDIR}", outputpath = "${WORKDIR}"))
        for f in self.get_output_files():
            wu.append('RESULT %s' % os.path.basename(f))
        wu.append("") # Make a trailing newline
        return '\n'.join(wu)
        
    @staticmethod
    def _open_or_not(fn, mode):
        """ If fn is a string, opens a file handle to a file with fn as 
        the name, using mode as the file mode. Otherwise returns fn.
        """
        if isinstance(fn, str):
            return open(fn, mode)
        else:
            return fn

    def run(self):
        ''' Runs the command and waits for termination '''

        # If we run a command locally, and file names were given for stdin,
        # stdout or stderr, we open the corresponding file and pass it 
        # to the command
        self.infile = self._open_or_not(self.stdin, "r")
        self.outfile = self._open_or_not(self.stdout, "w")
        self.errfile = self._open_or_not(self.stderr, "w")
        
        # print (self.make_command_array()
        # print ("%s.Program.run(): cmdline = %s" % (__file__, self.make_command_line()))
        # print ("Input files: %s" % ", ".join(self.get_input_files()))
        # print ("Output files: %s" % ", ".join(self.get_output_files()))
        self.child = cadocommand.Command(self.make_command_array(), 
                                         stdin=self.infile,
                                         stdout=self.outfile, 
                                         stderr=self.errfile)
    
    def wait(self):
        (rc, stdout, stderr) = self.child.wait()
        
        if isinstance(self.stdin, str):
            self.infile.close()
        if isinstance(self.stdout, str):
            self.outfile.close()
        if isinstance(self.stderr, str):
            self.errfile.close()
        
        return (rc, stdout, stderr)
    
    
class Polyselect2l(Program):
    binary = "polyselect2l"
    name = binary
    subdir = "polyselect"
    params_list = (
        Toggle("verbose", "v"), 
        Toggle("quiet", "q"), 
        Toggle("sizeonly", "r"), 
        Parameter("N"), 
        Parameter("threads", "t"), 
        Parameter("admin"), 
        Parameter("admax"), 
        Parameter("incr"), 
        Parameter("degree"), 
        Parameter("nq"), 
        Parameter("save", is_output_file = True), 
        Parameter("resume", is_input_file = True), 
        Parameter("maxnorm"), 
        Parameter("maxtime"), 
        Parameter("out"), 
        Parameter("printdelay", "s"),
        PositionalParameter("P")
        )


class MakeFB(Program):
    binary = "makefb"
    name = binary
    subdir = "sieve"
    params_list = (
        Toggle("nopowers"), 
        Parameter("poly", is_input_file = True), 
        Parameter("maxbits")
        )


class FreeRel(Program):
    binary = "freerel"
    name = binary
    subdir = "sieve"
    params_list = (
        Toggle("verbose", "v"), 
        Parameter("pmin"), 
        Parameter("pmax"),
        Parameter("poly", is_input_file = True)
        )

class Las(Program):
    binary = "las"
    name = binary
    subdir = "sieve"
    params_list = (
        Parameter("I"), 
        Parameter("poly", is_input_file = True), 
        Parameter("factorbase", "fb", is_input_file = True),
        Parameter("q0"), 
        Parameter("q1"), 
        Parameter("rho"), 
        Parameter("tdthresh"),
        Parameter("bkthresh"),
        Parameter("rlim"),
        Parameter("alim"),
        Parameter("lpbr"),
        Parameter("lpba"),
        Parameter("mfbr"),
        Parameter("mfba"),
        Parameter("rlambda"),
        Parameter("alambda"),
        Parameter("skewness", "S"),
        Toggle("verbose", "v"),
        Parameter("out", is_output_file = True),
        Parameter("threads", "mt"),
        Toggle("ratq")
        )

class Duplicates1(Program):
    binary = "dup1"
    name = binary
    subdir = "filter"
    params_list = (
        Parameter("out"), 
        Parameter("outfmt"), 
        Toggle("bzip", "bz"), 
        Parameter("only"), 
        Parameter("nslices_log", "n"), 
        Parameter("filelist", is_input_file = True),
        Parameter("basepath"),
        )


class Duplicates2(Program):
    binary = "dup2"
    name = binary
    subdir = "filter"
    params_list = (
        Toggle("remove", "rm"),
        Parameter("output_directory", "out"), 
        Parameter("filelist", is_input_file = True), 
        Parameter("rel_count", "K")
        )

class Purge(Program):
    binary = "purge"
    name = binary
    subdir = "filter"
    params_list = (
        Parameter("poly", is_input_file = True),
        Parameter("out", is_output_file = True), 
        Parameter("nrels"), 
        Parameter("outdel", is_output_file = True), 
        Parameter("sos", is_output_file = True), 
        Parameter("keep"), 
        Parameter("minpa"), 
        Parameter("minpr"), 
        Parameter("nprimes"), 
        Toggle("raw"), 
        Parameter("npthr"), 
        Parameter("inprel", is_input_file = True), 
        Parameter("outrel", is_output_file = True), 
        Parameter("npass"), 
        Parameter("required_excess")
        )

class Merge(Program):
    binary = "merge"
    name = binary
    subdir = "filter"
    params_list = (
        Parameter("mat", is_input_file = True), 
        Parameter("out", is_output_file = True), 
        Parameter("maxlevel"), 
        Parameter("keep"), 
        Parameter("skip"), 
        Parameter("forbw"), 
        Parameter("ratio"), 
        Parameter("coverNmax"), 
        Parameter("itermax"), 
        Parameter("resume"), 
        Parameter("mkztype"), 
        Parameter("wmstmax")
        )

# Todo: define is_input_file/is_output_file for remaining programs
class Replay(Program):
    binary= "replay"
    name = binary
    subdir = "filter"
    params_list = (
        Toggle("binary", prefix = "--"),
        Parameter("skip"),
        Parameter("purged"),
        Parameter("history", "his"),
        Parameter("index"),
        Parameter("out")
    )

class BWC(Program):
    binary = "bwc.pl"
    name = "bwc"
    subdir = "linalg/bwc"
    params_list = (
        Toggle("complete", prefix=":"),
        Toggle("wipeout", prefix=":"),
        Toggle("dryrun", "d"), 
        Toggle("verbose", "v"),
        ParameterEq("mpi"),
        ParameterEq("threads", "thr"),
        ParameterEq("mn"),
        ParameterEq("nullspace"),
        ParameterEq("interval"),
        ParameterEq("ys"),
        ParameterEq("matrix"),
        ParameterEq("wdir"),
        ParameterEq("mpiexec"),
        ParameterEq("hosts"),
        ParameterEq("hostfile"),
        ParameterEq("bwc_bindir")
        )

class Characters(Program):
    binary= "characters"
    name = binary
    subdir = "linalg"
    params_list = (
        Parameter("poly"),
        Parameter("purged"),
        Parameter("index"),
        Parameter("heavyblock"),
        Parameter("nchar"),
        Parameter("nthchar", "t"),
        Parameter("out"),
        PositionalParameter("wfile")
    )

class Sqrt(Program):
    binary = "sqrt"
    name = binary
    subdir = "sqrt"
    params_list = (
        Toggle("ab"),
        Toggle("rat"),
        Toggle("alg"),
        Toggle("gcd"),
        Parameter("poly"),
        Parameter("prefix"),
        Parameter("dep"),
        Parameter("purged"),
        Parameter("index"),
        Parameter("kernel", "ker")
        )
