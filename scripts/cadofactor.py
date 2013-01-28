#!/usr/bin/env python3
import logging
import subprocess

class ANSI(object):
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
        record.colour = self.__class__.colours[record.levelno]
        record.levelnametitle = record.levelname.title()
        record.nocolour = ANSI.NORMAL
        if hasattr(record, "indent"):
            record.padding = " " * record.indent
        else:
            record.padding = ""
        return super().format(record)

class FileFormatter(logging.Formatter):
    formatstr = \
       'PID%(process)s %(asctime)s %(levelnametitle)s:%(message)s' 

    def format(self, record):
        record.levelnametitle = record.levelname.title()
        return super().format(record)

    def __init__(self):
        super().__init__(fmt=self.__class__.formatstr)

class HandlerRoot(object):
    # Root class to strip off the logger argument before we reach object()
    def __init__(self, logger, **kwargs):
        super().__init__(**kwargs)

class ScreenHandler(HandlerRoot):
    def __init__(self, logger, lvl = logging.INFO, colour = True, **kwargs):
        h = logging.StreamHandler()
        h.setLevel(lvl)
        h.setFormatter(ScreenFormatter(colour = colour))
        logger.addHandler(h)
        super().__init__(logger, **kwargs)

class FileHandler(HandlerRoot):
    def __init__(self, logger, filelvl = logging.DEBUG, filename = None, **kwargs):
        if not filename is None:
            h = logging.FileHandler(filename)
            h.setLevel(filelvl)
            h.setFormatter(FileFormatter())
            logger.addHandler(h)
        super().__init__(logger, **kwargs)

class Logger(HandlerRoot):
    # We mustn't instantiate logging.Logger, but get a reference to a
    # pre-existing instance via getLogger(). Hence no inheritance from 
    # logging.Logger
    def __init__(self, **kwargs):
        self.logger = logging.getLogger(__name__)
        # Lowest possible threshold, so handlers get to see everything.
        # They do the filtering by themselves
        self.logger.setLevel(logging.DEBUG)
        # Init the various handlers which may exist as sibling classes, and
        # tell them what our logging.Logger instance is
        super().__init__(logger = self.logger, **kwargs)
    # Delegate all other method calls to the logging.Logger instance we 
    # have referenced in self.logger
    def __getattr__(self, name):
        return getattr(self.logger, name)

# Put the pieces together
class MyLogger(Logger, ScreenHandler, FileHandler):
    """ Logger that outputs to both screen and disk file. """
    pass

class Command(object):
    def __init__(self, command, logfile=None):
        # Run the command
        self.command = command
        if not logfile is None:
            f = open(logfile, "a")
            f.write(self.command + "\n")
        self.child = subprocess.Popen(self.command, stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE)
        if not logfile is None:
            f.write("# Child process has PID " + str(self.child.pid) + "\n")
            f.close()
        logger.info("Running command (PID=" + str(self.child.pid) + "): " + self.command)

    def wait(self):
        # Wait for command to finish executing, capturing stdout and stderr 
        # in output tuple
        (self.stdout, self.stderr) = self.child.communicate()
        logger.info("Exit status " + str(self.child.returncode) + " for PID " + str(self.child.pid))
        return self.child.returncode

if __name__ == '__main__':
    logger = MyLogger(filename = "log", filelvl = logging.DEBUG, lvl=logging.INFO)
#    logger = MyLogger(lvl=logging.INFO)
    logger.info("An Info Center!")
    logger.warn("Beware")
    logger.error("All hope abandon", extra={"indent" : 4})
    c = Command("ls", logfile = "commands")
    c.wait()
    print(c.stdout)
