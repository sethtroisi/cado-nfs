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
    def __init__(self, arg):
        self.arg = arg

class Parameter(Option):
    ''' Command line option that takes a parameter '''
    def map(self, value):
        return [self.prefix_char + self.arg, value]

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
    ''' Sub-classes must define class variables 
    binary: a string with the name of the binary executable file
    path: path to the directory of the binary file
    name: a name used internally in the script for this program, must  
      start with a letter and contain only letters and digits; if the binary 
      file name is of this form, it can be used
    params: a mapping that tell how to translate configuration parameters to 
      command line options
    '''
    
    path = "programs"
    
    @staticmethod
    def __shellquote(s):
        ''' Quote a command line argument
        
        Currently does it the hard way: always encloses the argument in single
        quotes, and escapes any single quotes that are part of the argument
        '''
        return "'" + s.replace("'", "'\\''") + "'"
    
    def __init__(self, args, kwargs):
        ''' Takes a list of positional parameters and a dictionary of command 
        line parameters 
        '''
        sep = os.sep;
        if self.path[-1] == sep:
            sep = ""
        self.command = [self.path + sep + self.binary]
        for key in kwargs:
            self.command += self.params[key].map(kwargs[key])
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
    params = {"verbose": Toggle("v"), 
              "quiet": Toggle("q"), 
              "sizeonly": Toggle("r"), 
              "threads": Parameter("t"), 
              "admin": Parameter("admin"), 
              "admax": Parameter("admax"), 
              "incr": Parameter("incr"), 
              "N": Parameter("N"), 
              "degree": Parameter("degree"), 
              "nq": Parameter("nq"), 
              "save": Parameter("save"), 
              "resume": Parameter("resume"), 
              "maxnorm": Parameter("maxnorm"), 
              "maxtime": Parameter("maxtime"), 
              "out": Parameter("out"), 
              "printdelay": Parameter("s")}

class MakeFB(Program):
    binary = "makefb"
    name = binary
    params = {}

class FreeRel(Program):
    binary = "freerel"
    name = binary
    params = {}

class Las(Program):
    binary = "las"
    name = binary
    params = {"q0": Parameter("q0"), 
              "q1": Parameter("q1")}
    # etc.

class Duplicates(Program):
    name = "dup"
    binary = "dup"
    name = binary
    params = {}

class Singletons(Program):
    binary = "singleton"
    name = binary
    params = {}

class Merge(Program):
    binary = "merge"
    name = binary
    params = {}

class BWC(Program):
    binary = "bwc.pl"
    name = "bwc"
    params = {}

class Sqrt(Program):
    binary = "sqrt"
    name = binary
    params = {}

