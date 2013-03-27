import os
import cadocommand
import cadologger

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
            self.command = self.command + ["-" + self.params[key], kwargs[key]]
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
        return self.child.wait()

class Polyselect2l(Program):
    binary = "polyselect2l"
    path = "."
    name = binary
    params = {"verbose": "v", 
              "quiet": "q", 
              "sizeonly": "r", 
              "threads": "t", 
              "admin": "admin", 
              "admax": "admax", 
              "incr": "incr", 
              "N": "N", 
              "degree": "degree", 
              "nq": "nq", 
              "save": "save", 
              "resume": "resume", 
              "maxnorm": "maxnorm", 
              "maxtime": "maxtime", 
              "out": "out", 
              "printdelay": "s"}

class MakeFB(Program):
    binary = "makefb"
    path = "."
    name = binary
    params = {}

class FreeRel(Program):
    binary = "freerel"
    path = "."
    name = binary
    params = {}

class Las(Program):
    binary = "las"
    path = "."
    name = binary
    params = {"q0": "q0", "q1": "q1"}
    # etc.

class Duplicates(Program):
    name = "dup"
    binary = "dup"
    path = "."
    name = binary
    params = {}

class Singletons(Program):
    binary = "singleton"
    path = "."
    name = binary
    params = {}

class Merge(Program):
    binary = "merge"
    path = "."
    name = binary
    params = {}

class BWC(Program):
    binary = "bwc.pl"
    path = "."
    name = "bwc"
    params = {}

class Sqrt(Program):
    binary = "sqrt"
    path = "."
    name = binary
    params = {}

