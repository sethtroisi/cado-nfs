import os
import sys
import platform
import subprocess
import abc
import cadocommand
import cadologger

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
    
    params_list = ("execpath", "execbin")
    path = "programs"
    
    @staticmethod
    def __shellquote(s):
        ''' Quote a command line argument
        
        Currently does it the hard way: always encloses the argument in single
        quotes, and escapes any single quotes that are part of the argument
        '''
        return "'" + s.replace("'", "'\\''") + "'"
    
    def __init__(self, args = None, kwargs = None, 
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
        
        path = kwargs.get("execpath", self.path)
        binary = kwargs.get("execbin", self.binary)
        
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr
        # Begin command line with program to execute
        self.exec_files = [path.rstrip(os.sep) + os.sep + binary]
        self.command = [self.exec_files[0]]
        self.input_files = []
        self.output_files = []
        if isinstance(self.stdin, str):
            self.input_files.append(self.stdin)
        if isinstance(self.stdout, str):
            self.output_files.append(self.stdout)
        if isinstance(self.stderr, str):
            self.output_files.append(self.stderr)
        
        # Add keyword command line parameters
        for p in self.params_list:
            key = p.get_key()
            if kwargs and key in kwargs:
                self.command += p.map(kwargs[key])
                if p.is_input_file:
                    self.input_files.append(kwargs[key])
                if p.is_output_file:
                    self.output_files.append(kwargs[key])
        
        # Add positional command line parameters
        # FIXME: how to identify input/output files here?
        if args:
            self.command += args
    
    @classmethod
    def get_params_list(cls):
        """ Return the accepted parameters as list of config file keywords """
        l = list(Program.params_list) + [opt.get_key() for opt in cls.params_list]
        return l
    
    def __str__(self):
        ''' Returns the command line as a string '''
        return " ".join([Program.__shellquote(arg) for arg in self.command])
    
    def as_array(self):
        ''' Returns the command line as a string array '''
        return self.command
    
    def get_input_files(self):
        return self.input_files
    
    def get_output_files(self):
        return self.output_files
    
    def get_exec_files(self):
        return self.exec_files
    
    def make_cmdline(self):
        cmdline = str(self)
        if isinstance(self.stdin, str):
            cmdline += ' < ' + self.stdin
        if isinstance(self.stdout, str):
            cmdline += ' > ' + self.stdout
        if isinstance(self.stderr, str):
            cmdline += ' 2> ' + self.stdout
        if self.stderr is subprocess.STDOUT:
            cmdline += ' 2>&1'
        return cmdline
    
        # TODO: Make workunit text from a program instance
        # This allows running program instances either directly with run(),
        # or adding them to the WU database table. 
        # Making a WU will require knowledge of which input file the program
        # needs, which output files it produces, and which command line needs 
        # to be run. Input files should probably be mandatory parameters to
        # __init__(). 
        # When we run a command locally: 
        #   We give the input files directly on the command line, with the 
        #   correct path (produced by whichever program generated the input 
        #   file, probably the path points at the working directory). 
        #   The binary file is called with the program path and the binary 
        #   file name. 
        #   Output files are placed in the working directory, under the given 
        #   output file name.
        # When we generate a WU: 
        #   We have to give the input file in a FILE line, place the file in 
        #   the upload directory, and use it on the command line with path
        #   ${DLDIR}/. 
        #   Binary files are given with EXECFILE lines, and are used on the 
        #   commandline with path ${DLDIR}/.
        #   Output files are produced in the client workdir, and are uploaded
        #   to the server upload dir. On the command line, they must be 
        #   referenced by ${WORKDIR}/, and on the server, the file name given
        #   to subsequent tasks must be that of the uploaded file in the 
        #   server's upload directory.
        pass
    
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
        
        # print (self.as_array())
        print ("%s.Program.run(): cmdline = %s" % (__file__, self.make_cmdline()))
        print ("Input files: %s" % ", ".join(self.input_files))
        print ("Output files: %s" % ", ".join(self.output_files))
        self.child = cadocommand.Command(self.as_array(), 
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
    params_list = (
        Toggle("nopowers"), 
        Parameter("poly", is_input_file = True), 
        Parameter("maxbits")
        )


class FreeRel(Program):
    binary = "freerel"
    name = binary
    params_list = (
        Toggle("verbose", "v"), 
        Parameter("pmin"), 
        Parameter("pmax"),
        Parameter("poly", is_input_file = True)
        )

class Las(Program):
    binary = "las"
    name = binary
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
    params_list = (
        Toggle("remove", "rm"),
        Parameter("output_directory", "out"), 
        Parameter("filelist", is_input_file = True), 
        Parameter("rel_count", "K")
        )

class Purge(Program):
    binary = "purge"
    name = binary
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

class Replay(Program):
    binary= "replay"
    name = binary
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
