#!/usr/bin/env python3
import logging

class ANSI(object):
    CSI = '\x1b[' # ANSI Control Sequence Introducer. Not the TV show
    SGR = 'm' # Set Graphics Rendition code
    NORMAL = CSI + '0' + SGR
    BLACK = CSI + '30' + SGR 
    GREY = CSI + '30;1'
    RED= CSI + '31' + SGR 
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
    colors = {
        logging.INFO : ANSI.BRIGHTGREEN,
        logging.WARNING : ANSI.BRIGHTYELLOW,
        logging.ERROR : ANSI.BRIGHTRED
        }

    formatstr = \
        '%(padding)s%(colour)s%(levelnametitle)s%(nocolour)s:%(message)s'

    def __init__(self):
        super().__init__(fmt=self.__class__.formatstr)

    def format(self, record):
        record.colour = ScreenFormatter.colors[record.levelno]
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
       '%(asctime)s PID%(process)s %(levelnametitle)s:%(message)s' 

    def __init__(self):
        super().__init__(fmt=self.__class__.formatstr)

class Logger(object):
    @staticmethod
    def getLogger(lvl = logging.INFO, filename=None):
        logger = logging.getLogger(__name__)
        logger.setLevel(lvl)
        ch = logging.StreamHandler()
        screenformatter = ScreenFormatter()
        ch.setFormatter(screenformatter)
        logger.addHandler(ch)
        if not filename is None:
            fh = logging.FileHandler(filename)
            fileformatter = FileFormatter()
            fh.setFormatter(fileformatter)
            logger.addHandler(fh)
        return logger

logger = Logger.getLogger(filename = "log")
logger.info("An Info Center!")
logger.warn("Beware")
logger.error("All hope abandon", extra={"indent" : 4})
