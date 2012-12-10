#!/usr/bin/env python3

# CGI script to handle file uploads

debug=1

import cgi, os
if debug > 0:
  import cgitb; cgitb.enable()
import sys
from tempfile import mkstemp
import wudb

def diag(level, text, var = None):
    if debug > level:
        if var == None:
            print (text, file=sys.stderr)
        else:
            print (text + str(var), file=sys.stderr)

# Global variable in this module so that other Python modules can import
# it and store the path to the upload directory in the shell environment 
# variable specified here
UPLOADDIRKEY="UPLOADDIR"
DBFILENAMEKEY="DBFILENAME"

def do_upload(db, input = sys.stdin, output = sys.stdout):
    diag(1, "Command line arguments:", sys.argv)
    diag(2, "Environment:", os.environ)

    try: # Windows needs stdio set for binary mode.
        import msvcrt
        msvcrt.setmode (0, os.O_BINARY) # stdin  = 0
        msvcrt.setmode (1, os.O_BINARY) # stdout = 1
    except ImportError:
        pass

    diag(1, "Reading POST data\n", "")
        
    form = cgi.FieldStorage(fp = input)

    diag (2, "Attributes for form: ", dir(form))
    diag (2, "form = ", form)

    header = "Content-Type: text/plain\r\n\r\n"

    # A nested FieldStorage instance holds the file
    message = None
    if "WUid" not in form:
        message = 'No "WUid" key found in POST data'
    if "clientid" not in form:
        message = 'No "clientid" key found in POST data'
    elif "results" not in form:
        message = 'No "results" key found in POST data'
    elif UPLOADDIRKEY not in os.environ:
        message = 'Script error: Environment variable ' + UPLOADDIRKEY + ' not set'
    elif not os.path.isdir(os.environ[UPLOADDIRKEY]):
        message = 'Script error: ' + os.environ[UPLOADDIRKEY] + ' is not a directory'
    else:
        WUid = form['WUid']
        clientid = form['clientid']
        fileitem = form['results']
        if 'errorcode' in form:
            errorcode = form['errorcode'].value
        else:
            errorcode = 0
        diag (2, "Attributes of fileitem: ", dir(fileitem))
        diag (2, "fileitem.name = ", fileitem.name);
        diag (2, "fileitem.filename = ", fileitem.filename);
        diag (2, "fileitem.type = ", fileitem.type);
        diag (2, "fileitem.type_options = ", fileitem.type_options);
        diag (2, "fileitem.disposition = ", fileitem.disposition);
        diag (2, "fileitem.disposition_options = ", fileitem.disposition_options);
        diag (2, "fileitem.headers = ", fileitem.headers);
        diag (2, "fileitem.encoding = ", fileitem.encoding);
        # Test if the result was uploaded and WUid was set:
        if not fileitem.filename:
            message = 'No file was uploaded'
        elif not WUid.value:
            message = 'No work unit was specified'
        elif not clientid.value:
            message = 'No client id was specified'

    if not message:
        filetuples = []
        # strip leading path from file name to avoid directory traversal attacks
        basename = os.path.basename(fileitem.filename)
        # Make a file name which does not exist yet and create the file
        (fd, filename) = mkstemp(suffix='', prefix=basename, 
            dir=os.environ[UPLOADDIRKEY])
        filestuple = (fileitem.filename, filename)
        filetuples.append(filestuple)
        
        wu = wudb.WuActiveRecord(db)
        try:
            wu.result(WUid.value, clientid.value, errorcode, filetuples)
        except wudb.StatusUpdateError:
            message = 'Workunit ' + WUid.value + 'was not currently assigned'

        # fd is a file descriptor, make a file object from it
        file = os.fdopen(fd, "wb")
        file.write(fileitem.file.read())
        bytes = file.tell()
        file.close()
        message = 'The file "' + basename + '" for work unit ' + WUid.value + \
            ' was uploaded successfully by client ' + clientid.value + \
            ' and stored as ' + filename + ', received ' + str(bytes) + ' bytes.'

    diag (0, sys.argv[0] + ': ', message)
    if output == sys.stdout:
        output.write(header + message)
    else:
        output.write((header + message).encode("ascii"))


# If this file is run directly by Python, call do_upload()
if __name__ == '__main__':
    if DBFILENAMEKEY not in os.environ:
        message = 'Script error: Environment variable ' + DBFILENAMEKEY + ' not set'
    dbfilename = os.environ[DBFILENAMEKEY]
    db = wudb.WuDb(dbfilename)
    do_upload(db)
