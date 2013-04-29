import os
import sys
import platform
import subprocess
import cadocommand
import cadologger

class Option(object):
    ''' Base class for command line options that may or may not take parameters
    '''
    # TODO: decide whether we need '-' or '/' as command line option prefix 
    # character under Windows. Currently: use "-'
    if False and platform.system() == "Windows":
        prefix_char = '/'
    else:
        prefix_char = '-'
    
    def __init__(self, config, arg = None):
        # Set config to the name of the parameter in the configuration file,
        # e.g., "verbose", and arg to the command line parameter, e.g., "v"
        # If arg is not given, its default is the same as config
        self.config = config
        if arg is None:
            self.arg = config
        else:
            self.arg = arg
    
    def get_key(self):
        return self.config

class PositionalParameter(Option):
    ''' Positional command line parameter '''
    def map(self, value):
        return [value]

class Parameter(Option):
    ''' Command line option that takes a parameter '''
    def map(self, value):
        return [self.prefix_char + self.arg, value]

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
            return [self.prefix_char + self.arg]
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
    
    path = "programs"
    
    @staticmethod
    def __shellquote(s):
        ''' Quote a command line argument
        
        Currently does it the hard way: always encloses the argument in single
        quotes, and escapes any single quotes that are part of the argument
        '''
        return "'" + s.replace("'", "'\\''") + "'"
    
    def __init__(self, args = None, kwargs = None, path = None, binary = None, 
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
        
        if not path is None:
            self.path = path
        if not binary is None:
            self.binary = binary
        
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr
        # Begin command line with program to execute
        self.command = [self.path.rstrip(os.sep) + os.sep + self.binary]
        
        # Add keyword command line parameters
        for p in self.params_list:
            key = p.get_key()
            if kwargs and key in kwargs:
                self.command += p.map(kwargs[key])
        
        # Add positional command line parameters
        if args:
            self.command += args
    
    @classmethod
    def params_dict(cls):
        """ Return the accepted parameters as a mapping from config file 
        keywords to Option instances  """
        return {p.get_key():p for p in cls.params_list}
    
    def __str__(self):
        ''' Returns the command line as a string '''
        return " ".join([Program.__shellquote(arg) for arg in self.command])
    
    def as_array(self):
        ''' Returns the command line as a string array '''
        return self.command

    def as_wu():
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

        self.child = cadocommand.Command(self.as_array(), 
                                         stdin=self.infile,
                                         stdout=self.outfile, 
                                         stderr=self.errfile)
    
    def wait(self):
        (rc, stdout, stderr) = self.child.wait()
        if stdout:
            print("Stdout: " + str(stdout))
        if stderr:
            print("Stderr: " + str(stderr))
        sys.stdout.flush()
        
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
        Parameter("save"), 
        Parameter("resume"), 
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
        Parameter("poly"), 
        Parameter("maxbits")
        )

class FreeRel(Program):
    binary = "freerel"
    name = binary
    params_list = (
        Toggle("verbose", "v"), 
        Parameter("pmin"), 
        Parameter("pmax"),
        Parameter("poly")
        )

class Las(Program):
    binary = "las"
    name = binary
    params_list = (
        Parameter("I"), 
        Parameter("poly"), 
        Parameter("factorbase", "fb"),
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
        Parameter("out"),
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
        Parameter("filelist"),
        Parameter("basepath"),
        )
    
    # cmd("$param{'bindir'}/filter/dup1 ".
    #     "-n $param{'nslices_log'} ".
    #     "-out $param{'prefix'}.nodup ".
    #     "-filelist $param{'prefix'}.newfilelist ".
    #     "-basepath $param{'wdir'} ",
    #     { cmdlog => 1, kill => 1,   
    #      logfile=>"$param{'prefix'}.dup1.log" });


class Duplicates2(Program):
    binary = "dup2"
    name = binary
    params_list = (
        Toggle("remove", "rm"),
        Parameter("output_directory", "out"), 
        Parameter("filelist"), 
        Parameter("rel_count", "K")
        )

class Purge(Program):
    binary = "purge"
    name = binary
    params_list = (
        Parameter("poly"),
        Parameter("out"), 
        Parameter("nrels"), 
        Parameter("outdel"), 
        Parameter("sos"), 
        Parameter("keep"), 
        Parameter("minpa"), 
        Parameter("minpr"), 
        Parameter("nprimes"), 
        Toggle("raw"), 
        Parameter("npthr"), 
        Parameter("inprel"), 
        Parameter("outrel"), 
        Parameter("npass"), 
        Parameter("required_excess")
        )

class Merge(Program):
    binary = "merge"
    name = binary
    params_list = (
        Parameter("mat"), 
        Parameter("out"), 
        Parameter("maxlevel"), 
        Parameter("keep"), 
        Parameter("skip"), 
        Parameter("forbw"), 
        Parameter("ratio"), 
        Parameter("coverNmax"), 
        Parameter("itermax"), 
        Parameter("resume"), 
        Parameter("mkztype"), 
        Parameter("wmstmax"), 
        )

class BWC(Program):
    binary = "bwc.pl"
    name = "bwc"
    params_list = (
        Toggle("dryrun", "d"), 
        Toggle("verbose", "v"), 
        ParameterEq("mpi"),
        ParameterEq("mn"),
        ParameterEq("nullspace"),
        ParameterEq("interval"),
        ParameterEq("ys"),
        ParameterEq("matrix"),
        ParameterEq("wdir"),
        ParameterEq("mpiexec"),
        ParameterEq("hosts"),
        ParameterEq("hostfile"),
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
        Parameter("ker")
        )
