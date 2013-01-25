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

    colourformatstr = \
        '%(padding)s%(colour)s%(levelnametitle)s%(nocolour)s:%(message)s'
    nocolourformatstr = \
        '%(padding)s%(levelnametitle)s:%(message)s'

    def __init__(self, colour = True):
        if colour:
            super().__init__(fmt=self.__class__.colourformatstr)
        else:
            super().__init__(fmt=self.__class__.nocolourformatstr)

    def format(self, record):
        record.colour = ScreenFormatter.colours[record.levelno]
        record.levelnametitle = record.levelname.title()
        record.nocolour = ANSI.NORMAL
        if hasattr(record, "indent"):
            assert isinstance(record.indent, int)
            record.padding = " " * record.indent
        else:
            record.padding = ""
        return super().format(record)

class FileFormatter(logging.Formatter):
    formatstr = \
       'PID%(process)s %(asctime)s %(levelnametitle)s:%(message)s' 

    def __init__(self):
        super().__init__(fmt=self.__class__.formatstr)

class Logger(object):
    @staticmethod
    def getLogger(lvl = logging.INFO, filename=None, filelvl = logging.INFO, colour=True):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(lvl)
        screenformatter = ScreenFormatter(colour=colour)
        ch.setFormatter(screenformatter)
        logger.addHandler(ch)
        if not filename is None:
            fh = logging.FileHandler(filename)
            fh.setLevel(filelvl)
            fileformatter = FileFormatter()
            fh.setFormatter(fileformatter)
            logger.addHandler(fh)
        return logger

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

logger = Logger.getLogger(filename = "log", filelvl = logging.DEBUG)
logger.info("An Info Center!")
logger.warn("Beware")
logger.error("All hope abandon", extra={"indent" : 4})
c = Command("ls", logfile = "commands")
c.wait()
print(c.stdout)
