#!/usr/bin/env python3

import http.server
import socketserver
import os
import sys
from workunit import Workunit
import wudb
import upload

# Get the shell environment variable name in which we should store the path 
# to the upload directory

upload_keywords = ['upload.py']
uploaddir='upload/' # Upload CGI script puts files in this directory
dbfilename='wudb'

class ThreadedHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    """Handle requests in a separate thread."""

debug = 1
def diag(level, text, var = None):
    if debug > level:
        if var is None:
            print (text, file=sys.stderr)
        else:
            print (text + str(var), file=sys.stderr)
        sys.stderr.flush()

class MyHandler(http.server.CGIHTTPRequestHandler):
    def send_body(self, body):
        self.wfile.write(bytes(body, "utf-8"))
        self.wfile.flush()

    def do_GET(self):
        """Generates a work unit if request is cgi-bin/getwu, otherwise calls
           parent class' do_GET()"""
        if self.is_cgi():
            if self.is_getwu():
                self.send_WU()
            elif self.is_getstatus():
                self.send_status()
            else:
                self.send_error(404, "GET for CGI scripts allowed only for work unit request")
        else:
            super().do_GET(self)
        sys.stdout.flush()

    def do_POST(self):
        """Set environment variable telling the upload directory 
           and call CGI handler to run upload CGI script"""
        if self.is_upload():
            os.environ[upload.UPLOADDIRKEY] = uploaddir
            os.environ[upload.DBFILENAMEKEY] = dbfilename
            if False:
                self.send_response(200, "Script output follows")
                upload.do_upload(db, input = self.rfile, output = self.wfile)
            else:
                http.server.CGIHTTPRequestHandler.do_POST(self)
        else:
            self.send_error(404, "POST request allowed only for uploads")
        sys.stdout.flush()

    def is_upload(self):
        """Test whether request is a file upload."""
        splitpath = http.server._url_collapse_path_split(self.path)
        if self.command == 'POST' and self.is_cgi() and splitpath[1] in upload_keywords:
            return True
        return False

    def is_getwu(self):
        """Test whether request is for a new WU."""
        filename=self.cgi_info[1].split("?", 1)[0]
        return self.command == 'GET' and filename in ['getwu']
    def is_getstatus(self):
        """Test whether request is for a a status page."""
        filename=self.cgi_info[1].split("?", 1)[0]
        return self.command == 'GET' and filename in ['status']

    def send_WU(self):
        filename = self.cgi_info[1]
        if not "?" in filename:
            return self.send_error(400, "No query string given")
        (filename, query) = self.cgi_info[1].split("?", 1)
        if query.count("=") != 1 or "?" in query or "&" in query:
            return self.send_error(400, "Bad query string in request")
        (key, clientid) = query.split("=")
        if key != "clientid":
            return self.send_error(400, "No client id specified")
        if not clientid.isalnum():
            return self.send_error(400, "Malformed client id specified")
        
        # wu = wudb.WuActiveRecord(db)
        wu_text = db_pool.assign(clientid)
        if not wu_text:
            return self.send_error(404, "No work available")
        
        self.log_message("Sending work unit " + Workunit(wu_text).get_id() + " to client " + clientid)
        # wu_text = wu.get_wu()
        self.send_response(200)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", len(wu_text))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.send_body(wu_text)

    def send_status():
        filename = self.cgi_info[1]
        (filename, query) = self.cgi_info[1].split("?", 1)

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

    db_pool = wudb.DbThreadPool(dbfilename, 1)

    HandlerClass = MyHandler
    httpd = ServerClass((HTTP, PORT), HandlerClass)
    httpd.server_name = "test"

    print ("serving at " + HTTP + ":" + str(PORT))
    httpd.serve_forever()
    db.terminate()
