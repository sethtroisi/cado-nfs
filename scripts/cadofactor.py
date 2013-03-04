#!/usr/bin/env python3
import logging
import subprocess

class ANSI(object):
    """ Class defining some ANSI control sequences, for example for 
    changing text colour """
    CSI = '\x1b[' # ANSI Control Sequence Introducer. Not the TV show
    SGR = 'm' # Set Graphics Rendition code
    NORMAL = CSI + '0' + SGR
    BLACK = CSI + '30' + SGR 
    GREY = CSI + '30;1'
    RED = CSI + '31' + SGR 
    BRIGHTRED = CSI + '31;1' + SGR
    GREEN = CSI + '32' + SGR
    BRIGHTGREEN = CSI + '32;1' + SGR
    YELLOW = CSI + '33' + SGR
    BRIGHTYELLOW = CSI + '33;1' + SGR
    BLUE = CSI + '34' + SGR
    BRIGHTBLUE = CSI + '34;1' + SGR
    MAGENTA = CSI + '35' + SGR
    BRIGHTMAGENTA = CSI + '35;1' + SGR
    CYAN = CSI + '36' + SGR
    BRIGHTCYAN = CSI + '36;1' + SGR
    # There is another grey with code "37" - white without intensity
    # Not sure if it is any different from "30;1" aka "bright black"
    WHITE = CSI + '37;1' + SGR

class ScreenFormatter(logging.Formatter):
    """ Class for formatting logger records for screen output, optionally
    with colorized logger level name (like cadofct.pl used to). """
    colours = {
        logging.INFO : ANSI.BRIGHTGREEN,
        logging.WARNING : ANSI.BRIGHTYELLOW,
        logging.ERROR : ANSI.BRIGHTRED
    }

    # Format string that switches to a different colour (with ANSI code 
    # specified in the 'colour' key of the log record) for the log level name, 
    # then back to default text rendition (ANSI code in 'nocolour')
    colourformatstr = \
        '%(padding)s%(colour)s%(levelnametitle)s%(nocolour)s:%(message)s'
    # Format string that does not use colour changes
    nocolourformatstr = \
        '%(padding)s%(levelnametitle)s:%(message)s'

    def __init__(self, colour = True):
        if colour:
            super().__init__(fmt=self.__class__.colourformatstr)
        else:
            super().__init__(fmt=self.__class__.nocolourformatstr)

    def format(self, record):
        # Add attributes to record that our format string expects
        if record.levelno in self.__class__.colours:
            record.colour = self.__class__.colours[record.levelno]
        else:
            record.colour = ANSI.NORMAL
        record.levelnametitle = record.levelname.title()
        record.nocolour = ANSI.NORMAL
        if hasattr(record, "indent"):
            record.padding = " " * record.indent
        else:
            record.padding = ""
        return super().format(record)

class FileFormatter(logging.Formatter):
    """ Class for formatting a log record for writing to a log file. No colours 
    here, but we add the process ID and a time stamp """
    formatstr = \
       'PID%(process)s %(asctime)s %(levelnametitle)s:%(message)s' 

    def format(self, record):
        record.levelnametitle = record.levelname.title()
        return super().format(record)

    def __init__(self):
        super().__init__(fmt=self.__class__.formatstr)

class CmdFileFormatter(logging.Formatter):
    """ Class for formatting a log record for writing to a command log file. 
        No colours here, but we add the process ID, the spanwed process id, 
        and a time stamp """
    formatstr = \
       '# PPID%(process)s PID%(childpid)s %(asctime)s\n%(message)s' 

    def __init__(self):
        super().__init__(fmt=self.__class__.formatstr)

class ScreenHandler(logging.StreamHandler):
    def __init__(self, lvl = logging.INFO, colour = True, **kwargs):
        super().__init__(**kwargs)
        self.setLevel(lvl)
        self.setFormatter(ScreenFormatter(colour = colour))

class FileHandler(logging.FileHandler):
    def __init__(self, filename, lvl = logging.DEBUG, **kwargs):
        super().__init__(filename, **kwargs)
        self.setLevel(lvl)
        self.setFormatter(FileFormatter())

class CmdFileHandler(logging.FileHandler):
    def __init__(self, filename , **kwargs):
        super().__init__(filename, **kwargs)
        self.setLevel(0) # FIXME: how do I make it process only record with level EQUAL to some value?
        self.setFormatter(CmdFileFormatter())

class Logger(object):
    """ Class which gets a logger with name equal to the module name (i.e., as 
        stored in __name__) upon instantiation and sets the logging level to 
        DEBUG. Other method calls are passed though to the logger """
    CMDLEVEL = 51
    def __init__(self):
        # We mustn't instantiate logging.Logger, but get a reference to a
        # pre-existing instance via getLogger(). Hence no inheritance from 
        # logging.Logger
        self.logger = logging.getLogger(__name__)
        # Use level of NOTSET, so handlers get to see everything.
        # They do the filtering by themselves
        self.logger.setLevel(logging.DEBUG)
    
    def cmd(self, msg, *args, **kwargs):
        """ Log a message with a level of Logger.CMDLEVEL """
        self.log(self.__class__.CMDLEVEL, msg, *args, **kwargs)
    
    # Delegate all other method calls to the logging.Logger instance we 
    # have referenced in self.logger
    def __getattr__(self, name):
        return getattr(self.logger, name)

class Command(object):
    def __init__(self, args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.logger = Logger()
        # Convert args array to a string for printing if necessary
        if isinstance(self.args, str):
            self.cmdline = self.args
        else:
            self.cmdline = " ".join(self.args)

        self.child = subprocess.Popen(self.args, stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE, **self.kwargs)
        
        self.logger.info("Running command: " + self.cmdline)
        self.logger.cmd(self.cmdline, extra={"pid": self.child.pid})

    def wait(self):
        # Wait for command to finish executing, capturing stdout and stderr 
        # in output tuple
        (self.stdout, self.stderr) = self.child.communicate()

        if self.child.returncode == 0:
            logger.info("Process with PID " + str(self.child.pid) + " finished successfully")
        else:
            logger.error("Process with PID " + str(self.child.pid) + " finished with return code " + str(self.child.returncode))
        self.returncode = self.child.returncode
        return self.returncode

class RemoteCommand(Command):
    ssh="/usr/bin/ssh"
    ssh_options = {
        "ConnectTimeout": 30,
        "ServerAliveInterval": 10,
        "PasswordAuthentication": "no"
    }
    def __init__(self, command, host, port = None, ssh_options = None, **kwargs):
        ssh_command = [self.__class__.ssh]
        options = self.__class__.ssh_options.copy()
        if not ssh_options is None:
            options.update(ssh_options)
        if not port is None:
            ssh_command += ["-p", str(port)];
        for (opt, val) in options.items():
            if not val is None:
                ssh_command += ["-o", opt + "=" + str(val)]
        ssh_command.append(host)
        if isinstance(command, str):
            ssh_command.append(command)
        else:
            ssh_command += command
        super().__init__(ssh_command, **kwargs)

class SendFile(Command):
    rsync="/usr/bin/rsync"
    rsync_options = []
    def __init__(self, localfile, hostname, hostpath, port = None, rsync_options = None, **kwargs):
        if hostname != "localhost":
            target = hostname + ":"
        if not port is None:
            target += str(port)
        target += hostpath
        copy_command = [self.__class__.rsync] + self.__class__.rsync_options + [localfile, target]
        super().__init__(copy_command, **kwargs)

if __name__ == '__main__':
    logger = Logger()
    logger.addHandler(ScreenHandler(lvl = logging.INFO))
    logger.addHandler(FileHandler(filename = "log", lvl = logging.DEBUG))

    logger.info("An Info Center!")
    logger.warn("Beware")
    logger.error("All hope abandon", extra={"indent" : 4})

    c = Command(["ls", "/"])
    rc = c.wait()
    print("Stdout: " + str(c.stdout, encoding="utf-8"))
    print("Stderr: " + str(c.stderr, encoding="utf-8"))
    del(c)

    c = RemoteCommand(["ls", "/"], "localhost")
    rc = c.wait()
    print("Stdout: " + str(c.stdout, encoding="utf-8"))
    print("Stderr: " + str(c.stderr, encoding="utf-8"))
    del(c)
