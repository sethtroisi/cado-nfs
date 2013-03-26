import cadocommand
import cadologger

class Program(object):
    ''' Sub-classes must define class variables 'binary' and 'accepts', 
    which store the file name of the executable, and a list of strings 
    of command line parameters the program accepts, respectively.
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
        self.command = [self.binary]
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

class Las(Program):
    binary = "las"
    params = {"q0": "q0", "q1": "q1"}
    # etc.
