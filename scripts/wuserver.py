#!/usr/bin/env python3

import http.server
import socketserver
import os
from Workunit import Workunit

# Get the shell environment variable name in which we should store the path 
# to the upload directory
from upload import UPLOADDIRKEY

upload_keywords = ['upload.py']
uploaddir='upload/' # Upload CGI script puts files in this directory

class ThreadedHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    """Handle requests in a separate thread."""

class Workunits():
    WU_counter = 0
    
    def get_new(self):
        """Returns the contents of a WU file for an available WU"""
        if self.WU_counter > 500:
            return None
        WU = "WORKUNIT WU" + str(self.WU_counter) + "\r\n"
        WU = WU + "EXECFILE ecm\r\n"
        WU = WU + "FILE c200\r\n"
        WU = WU + "COMMAND $DLDIR/ecm 1e5 < $DLDIR/c200 > $WORKDIR/output\r\n"
        WU = WU + "COMMAND gzip -9 $WORKDIR/output\r\n"
        WU = WU + "RESULT output.gz\r\n"
        self.WU_counter = self.WU_counter + 1
        return WU;

WU = Workunits()

class MyHandler(http.server.CGIHTTPRequestHandler):
    def send_body(self, body):
        self.wfile.write(bytes(body, "utf-8"))
        self.wfile.flush()

    def do_GET(self):
        """Generates a work unit if request is cgi-bin/getwu, otherwise calls
           parent class' do_GET()"""
        if self.is_getwu():
            self.send_WU()
        elif self.is_cgi():
            self.send_error(404, "GET for CGI scripts allowed only for work unit request")
        else:
            http.server.CGIHTTPRequestHandler.do_GET(self)

    def send_WU(self):
        new_WU = WU.get_new()
        if new_WU == None:
            return self.send_error(404)
        self.send_response(200)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", str(len(new_WU)))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.send_body(new_WU)

    def do_POST(self):
        """Set environment variable telling the upload directory 
           and call CGI handler to run upload CGI script"""
        os.environ[UPLOADDIRKEY] = uploaddir
        if self.is_upload():
            http.server.CGIHTTPRequestHandler.do_POST(self)
        else:
            self.send_error(404, "POST request allowed only for uploads")

    def is_upload(self):
        """Test whether request is a file upload."""
        splitpath = http.server._url_collapse_path_split(self.path)
        if self.command == 'POST' and self.is_cgi() and splitpath[1] in upload_keywords:
            return True
        return False

    def is_getwu(self):
        """Test whether request is for a new WU."""
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
