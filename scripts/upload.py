# CGI script to handle file uploads

#!/usr/bin/env python3
import cgi, os
import cgitb; cgitb.enable()
from tempfile import mkstemp

try: # Windows needs stdio set for binary mode.
    import msvcrt
    msvcrt.setmode (0, os.O_BINARY) # stdin  = 0
    msvcrt.setmode (1, os.O_BINARY) # stdout = 1
except ImportError:
    pass

form = cgi.FieldStorage()

# A nested FieldStorage instance holds the file
if "file" not in form:
    message = 'No "file" key found'
else:
    fileitem = form['file']

# Test if the file was uploaded
    if fileitem.filename:
        # strip leading path from file name to avoid directory traversal attacks
        basename = os.path.basename(fileitem.filename)
        # Make a file name which does not exist yet and create/open the file
        (fd, filename) = mkstemp(suffix='', prefix=basename, dir='upload/')
        file = os.fdopen(fd, "wb")
        # self.log_message("Received POST of " + str(content_length) + " bytes, storing in file " + filename)
        file.write(fileitem.file.read())
        bytes = file.tell();
        file.close()
        message = 'The file "' + basename + '" was uploaded successfully, received ' + \
            str(bytes) + ' bytes.'
    else:
        message = 'No file was uploaded'

print ("""\
Content-Type: text/html\r\n\r\n
<html><body>
<p>%s</p>
</body></html>
""" % (message,))
