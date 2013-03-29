import os
import sys
import platform
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
        else:
            return []


class Program(object):
    ''' Base class that represents programs of the CADO suite

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
    
    def __init__(self, args, kwargs, path = None, binary = None):
        ''' Takes a list of positional parameters and a dictionary of command 
        line parameters 
        '''
        
        if not path is None:
            self.path = path
        if not binary is None:
            self.binary = binary
        
        # Start command line with program to execute
        sep = os.sep;
        if self.path[-1] == sep:
            sep = ""
        self.command = [self.path + sep + self.binary]

        # Add keyword command line parameters
        params_dict = {p.config:p for p in params_list}
        for key in kwargs:
            self.command += params_dict[key].map(kwargs[key])
        # Add positional command line parameters
        self.command += args
    
    def __str__(self):
        ''' Returns the command line as a string '''
        return " ".join([Program.__shellquote(arg) for arg in self.command])
    
    def as_array(self):
        ''' Returns the command line as a string array '''
        return self.command
    
    def run(self):
        ''' Runs the command and waits for termination '''
        self.child = cadocommand.Command(self.as_array())
    
    def wait(self):
        r = self.child.wait()
        if len(r[1]) > 0:
            print("Stdout: " + str(r[1]))
        if len(r[2]) > 0:
            print("Stderr: " + str(r[2]))
        sys.stdout.flush()
        return r

class Polyselect2l(Program):
    binary = "polyselect2l"
    name = binary
    params_list = (
        Toggle("verbose", "v"), 
        Toggle("quiet", "q"), 
        Toggle("sizeonly", "r"), 
        Parameter("threads", "t"), 
        Parameter("admin"), 
        Parameter("admax"), 
        Parameter("incr"), 
        Parameter("N"), 
        Parameter("degree"), 
        Parameter("nq"), 
        Parameter("save"), 
        Parameter("resume"), 
        Parameter("maxnorm"), 
        Parameter("maxtime"), 
        Parameter("out"), 
        Parameter("printdelay", "s")
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
    # etc.

class Duplicates1(Program):
    binary = "dup1"
    name = binary
    params_list = (
        Parameter("out"), 
        Parameter("outfmt"), 
        Toggle("bzip", "bz"), 
        Parameter("only"), 
        Parameter("nslices", "n"), 
        )

class Duplicates2(Program):
    binary = "dup2"
    name = binary
    params_list = (
        Toggle("remove", "rm"),
        Parameter("out"), 
        Parameter("filelist"), 
        Parameter("K")
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
