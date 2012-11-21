#!/usr/bin/env python3

# CGI script to handle file uploads

debug=1

import cgi, os
if debug > 0:
  import cgitb; cgitb.enable()
from sys import stderr
from sys import argv
from tempfile import mkstemp

# Global variable in this module so that other Python modules can import
# it and store the path to the upload directory in the shell environment 
# variable specified here
UPLOADDIRKEY="UPLOADDIR"

def do_upload():
    if debug > 1:
        stderr.write("Command line arguments:\n")
        stderr.write(str(argv))
    if debug > 2:
        stderr.write("\nEnvironment:\n")
        for e in os.environ.items():
            stderr.write(str(e) + "\n")
        stderr.write("\n")

    try: # Windows needs stdio set for binary mode.
        import msvcrt
        msvcrt.setmode (0, os.O_BINARY) # stdin  = 0
        msvcrt.setmode (1, os.O_BINARY) # stdout = 1
    except ImportError:
        pass

    form = cgi.FieldStorage()

    print ("Content-Type: text/plain\r\n\r\n",)

    # A nested FieldStorage instance holds the file
    message = None
    if "WUid" not in form:
        message = 'No "WUid" key found in POST data'
    elif "results" not in form:
        message = 'No "results" key found in POST data'
    elif UPLOADDIRKEY not in os.environ:
        message = 'Script error: Environment variable ' + UPLOADDIRKEY + ' not set'
    elif not os.path.isdir(os.environ[UPLOADDIRKEY]):
        message = 'Script error: ' + os.environ[UPLOADDIRKEY] + ' is not a directory'
    else:
        WUid = form['WUid']
        fileitem = form['results']
        # Test if the result was uploaded and WUid was set:
        if not fileitem.filename:
            message = 'No file was uploaded'
        elif not WUid.value:
            message = 'No work unit was specified'

    if not message:
        # strip leading path from file name to avoid directory traversal attacks
        basename = os.path.basename(fileitem.filename)
        # Make a file name which does not exist yet and create the file
        (fd, filename) = mkstemp(suffix='', prefix=basename, dir='upload/')
        # fd is a file descriptor, make a file object from it
        file = os.fdopen(fd, "wb")
        file.write(fileitem.file.read())
        bytes = file.tell()
        file.close()
        message = 'The file "' + basename + '" for work unit ' + WUid.value + \
            ' was uploaded successfully and stored as ' + filename + \
            ', received ' + str(bytes) + ' bytes.'

    stderr.write(argv[0] + ': ' + message + '\n')
    print(message)


# If this file is run directly by Python, call do_upload()
if __name__ == '__main__':
    do_upload()
