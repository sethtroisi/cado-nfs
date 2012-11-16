#!/usr/bin/env python3

import http.server
import socketserver
import select
import threading
import pickle
import os
from shutil import copyfileobj
from tempfile import mkstemp

upload_keywords = ['myupload']

newfile_lock = threading.Lock()

class ThreadedHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    """Handle requests in a separate thread."""

# What do we do with excess data? Apache interprets it as an additional 
# request, which usually does not look like a proper request (exploitable?) 
# and responds with '501 Method Not Implemented'. We just keep the data 
# according to 'Content-Length' and assume that's the data the client wanted 
# us to have
def copy_data(input, output, length):
    """ Copies at most length bytes from input to output """
    buffser_size = 65536
    timeout = 30.0 # 30 seconds
    copied = 0
    while copied < length:
        if input.closed:
            break
        # Wait until some data has arrived to avoid busy waiting here
        s = select.select([input], [], [], timeout)
        if not s[0]:
            # Timeout reached. Return to caller
            break
        copy_now = min(buffser_size, length - copied)
        data = input.read(copy_now)
        # print ("Read " + str(len(data)) + " bytes")
        # If we got 0 bytes even tho select returned, we assume we're at EOS
        if len(data) == 0:
            break
        output.write(data)
        copied = copied + len(data)
    return copied


class MyHandler(http.server.CGIHTTPRequestHandler):
    def do_GET(self):
        """Generates a work unit if request is cgi-bin/getwu, otherwise calls
           parent class' do_GET()"""
        if is_getwu():
            make_wu()
        else:
            http.server.CGIHTTPRequestHandler.do_POST(self)

    def do_POST(self):
        """Set environment variable telling the upload directory 
           and call CGI handler to run upload CGI script"""
        http.server.CGIHTTPRequestHandler.do_POST(self)

    def is_upload(self):
        """Test whether request is a file upload."""
        splitpath = http.server._url_collapse_path_split(self.path)
        if self.command == 'POST' and self.is_cgi() and splitpath[1] in upload_keywords:
            return True
        return False

    def is_getwu(self):
        """Test whether request is for a WU."""
        splitpath = http.server._url_collapse_path_split(self.path)
        if self.command == 'GET' and self.is_cgi() and splitpath[1] in ['getwu']:
            return True
        return False

    def discard_excess_data(self):
        """ Discard excess data on stream, like CGIHTTPRequestHandler.run_cgi()
             does. Returns number of bytes that were discarded. """
        discarded = 0
        while select.select([self.rfile], [], [], 0)[0]:
            if not self.rfile.read(1):
                break
            discarded = discarded + 1
        return discarded


    # Totally ignores the HTTP protocol standard for POST data. Do not use.
    def upload(self):
        if not 'Content-Length' in self.headers:
            self.send_error(411) # 'Content-Length missing' error
            return
        try:
            content_length = int(self.headers['Content-Length'])
        except ValueError:
            # Bad request error
            self.send_error(400, "Invalid Content-Length: " + 
                self.headers['Content-Length'])
            return
        (fd, filename) = mkstemp(suffix='', prefix='infile', dir='upload/')
        file = os.fdopen(fd, "wb")
        self.log_message("Received POST of " + str(content_length) + " bytes, "
            "storing in file " + filename)
        copied_length = copy_data(self.rfile, file, content_length)
        file.close()
        discarded = self.discard_excess_data()
        if discarded:
            self.log_error(str(discarded) + " surplus bytes after POST data")
        # print(data)
        if copied_length < content_length:
            # Apache sends '500 Internal Server Error' when length of data  
            # in POST is shorter than 'Content-Length' indicates.
            # FIXME: Should we remove the file?
            msg = "Received only " + str(copied_length) + \
                " bytes of file data when Content-Length: indicated " + \
                str(content_length)
            self.send_error(500, msg)
            return
        file.close()
        self.send_response(200)
        self.send_header('Content-type', b'text/plain')
        body = 'Received ' + str(content_length) + ' bytes of file data.\r\n'
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(bytes(body, encoding='utf-8'))

if __name__ == '__main__':
    from sys import argv
    PORT = 8001
    HTTP = "" # address on which to listen
    want_threaded = 1

    if want_threaded:
        ServerClass = ThreadedHTTPServer
    else:
        ServerClass = http.server.HTTPServer

    if argv[1:]:
        PORT = int(argv[1])

    HandlerClass = MyHandler
    httpd = ServerClass((HTTP, PORT), HandlerClass)
    httpd.server_name = "test"

    print ("serving at " + HTTP + ":" + str(PORT))
    httpd.serve_forever()
