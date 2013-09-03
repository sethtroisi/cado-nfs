#!/usr/bin/env python3

# CGI script to handle file uploads

DEBUG = 1

import cgi, os
import sys
from tempfile import mkstemp
from shutil import copyfileobj
import wudb

def diag(level, text, var = None):
    if DEBUG > level:
        if var == None:
            print (text, file=sys.stderr)
        else:
            print ("%s%s" % (text, var), file=sys.stderr)
        sys.stderr.flush()

def analyze(level, name, obj):
    """ Dump tons of internal data about an object """
    if DEBUG > level:
        diag (level, "*** Content dump of ", name)
        diag (level, "type(%s): " % name, type(obj))
        diag (level, "dir(%s): " % name, dir(obj))
        diag (level, "%s.__str__() = " % name, obj)
        diag (level, "%s.__repr__() = %r" % (name, obj))
        for name2 in dir(obj):
            diag (level, "%s.%s = " % (name, name2), getattr(obj, name2))

# Global variable in this module so that other Python modules can import
# it and store the path to the upload directory in the shell environment
# variable specified here
UPLOADDIRKEY = "UPLOADDIR"
DBFILENAMEKEY = "DBFILENAME"

def do_upload(dbfilename, inputfp = sys.stdin, output = sys.stdout):
    diag(1, "Command line arguments:", sys.argv)
    diag(2, "Environment:", os.environ)

    try: # Windows needs stdio set for binary mode.
        # pylint: disable=F0401
        import msvcrt
        # pylint: disable=E1101
        msvcrt.setmode (0, os.O_BINARY) # stdin  = 0
        msvcrt.setmode (1, os.O_BINARY) # stdout = 1
    except ImportError:
        pass

    diag(1, "Reading POST data")
    form = cgi.FieldStorage(fp = inputfp)
    diag(1, "Finished reading POST data")
    analyze(2, "form", form)

    charset = "utf-8"
    header = "Content-Type: text/plain; charset=" + charset + "\r\n\r\n"

    # A nested FieldStorage instance holds the file
    message = None
    if "WUid" not in form:
        message = 'No "WUid" key found in POST data'
    if "clientid" not in form:
        message = 'No "clientid" key found in POST data'
    elif UPLOADDIRKEY not in os.environ:
        message = 'Script error: Environment variable %s not set' \
                  % UPLOADDIRKEY
    elif not os.path.isdir(os.environ[UPLOADDIRKEY]):
        message = 'Script error: %s is not a directory' \
                  % os.environ[UPLOADDIRKEY]
    else:
        wuid = form['WUid']
        clientid = form['clientid']
        if 'errorcode' in form:
            errorcode = form['errorcode'].value
            diag(1, "errorcode = ", errorcode)
        else:
            errorcode = None
        if 'failedcommand' in form:
            failedcommand = form['failedcommand'].value
            diag(1, "failedcommand = ", failedcommand)
        else:
            failedcommand = None

        # Test if wuid and clientid was set:
        if not wuid.value:
            message = 'No workunit was specified'
        elif not clientid.value:
            message = 'No client id was specified'

        diag(1, "wuid = ", wuid.value)
        diag(1, "clientid = ", clientid.value)

    if not message:
        filetuples = []
        if 'results' in form:
            fileitems = form['results']
            if isinstance(fileitems, cgi.FieldStorage):
                fileitems = [fileitems] # Make it iterable
        else:
            fileitems = []
            diag(1, 'No "results" form found')
        analyze (2, "fileitems", fileitems)

        message = ""
        for fileitem in fileitems:
            analyze (2, "f", fileitem)
            diag(1, "Processing file ", fileitem.filename)
            # strip leading path from file name to avoid directory traversal
            # attacks
            basename = os.path.basename(fileitem.filename)
            # Split extension from file name. We need to preserve the 
            # file extension so that, e.g., gzipped files can be identified
            # as such
            (basename, suffix) = os.path.splitext(basename)
            # Make a file name which does not exist yet and create the file
            (filedesc, filename) = mkstemp(prefix=basename + '.',
                suffix=suffix, dir=os.environ[UPLOADDIRKEY])
            diag(1, "output filename = ", filename)
            filestuple = (fileitem.filename, filename)
            if False:
                filestuple = (fileitem.filename, os.path.basename(filename))
            filetuples.append(filestuple)
            
            # fd is a file descriptor, make a file object from it
            diag(1, "Getting file object for temp file")
            file = os.fdopen(filedesc, "wb")
            diag(1, "Writing data to temp file")
            copyfileobj(fileitem.file, file)
            nr_bytes = file.tell()
            diag(1, "Wrote %d bytes" % nr_bytes)
            diag(1, "Closing file")
            file.close()
            
            # Example output:
            # upload.py: The file "testrun.polyselect.0-5000" for workunit
            # testrun_polyselect_0-5000 was uploaded successfully by client
            # localhost and stored as /localdisk/kruppaal/work/testrun.upload/
            # testrun.polyselect.0-5000.kcudj7, received 84720 bytes.
            message = message + 'The file "%s" for workunit %s was uploaded ' \
            'successfully by client %s and stored as %s, received %d bytes.\n' \
            % (basename, wuid.value, clientid.value, filename, nr_bytes)
        diag(1, "Getting WuAccess object")
        wuar = wudb.WuAccess(dbfilename)
        diag(1, "Got WuAccess object. Calling .result()")
        try:
            wuar.result(wuid.value, clientid.value, filetuples, errorcode,
                        failedcommand)
        except wudb.StatusUpdateError:
            message = 'Workunit ' + wuid.value + 'was not currently assigned'
        else:
            message = message + 'Workunit ' + wuid.value + ' completed.\n'
        diag(1, "Finished .result()")

    diag (0, sys.argv[0] + ': ', message.rstrip("\n"))
    if output == sys.stdout:
        output.write(header + message)
    else:
        output.write((header + message).encode(charset))


# If this file is run directly by Python, call do_upload()
if __name__ == '__main__':
    if DEBUG > 0:
        import cgitb
        cgitb.enable()

    if DBFILENAMEKEY not in os.environ:
        print ('Script error: Environment variable %s not set'
               % DBFILENAMEKEY)
        sys.exit(1)

    DBFILENAME = os.environ[DBFILENAMEKEY]
    diag(1, "About to call do_upload()")
    do_upload(DBFILENAME)
