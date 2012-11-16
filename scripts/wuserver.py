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
