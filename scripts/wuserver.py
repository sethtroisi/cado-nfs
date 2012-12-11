#!/usr/bin/env python3

import http.server
import socketserver
import threading
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

class DbThread():
    def __init__(self, dbfilename):
        self.cv = threading.Condition()
        self.lock = threading.Lock()
        self.dbfilename = dbfilename
        self.serving = False
        self.terminate = False
        self.thread = threading.Thread(target = self.serve_forever)
        self.thread.start()
    
    def assign(self, clientid):
        self.lock.acquire()
        self.cv.acquire()
        self.clientid = clientid
        self.cv.notify()
        self.cv.wait()
        # Make a thread-local copy of the WU while we are in locked region
        local = threading.local()
        if self.wu:
            local.wu = self.wu[:]
        else:
            local.wu = None
        self.cv.release()
        self.lock.release()
        return local.wu
    
    def serve_forever(self):
        self.cv.acquire()
        if self.serving:
            self.cv.release()
            # FIXME: find good exception to raise here
            raise Exception
        self.serving = True
        db = wudb.WuDb(self.dbfilename)
        wuar = wudb.WuActiveRecord(db)
        while True:
            self.cv.wait() # Give up lock and wait for notify, then get lock again
            if self.terminate == True:
                break
            if not wuar.assign(self.clientid):
                self.wu = None
            else:
                self.wu = wuar.get_wu()
            self.cv.notify()
        db.close()
        cv.release()    
    
    def finish(self):
        self.cv.acquire()
        self.terminate = True
        self.cv.notify()
        self.cv.wait()
        self.thread.join()


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
        wu_text = db.assign(clientid)
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

    def is_upload(self):
        """Test whether request is a file upload."""
        splitpath = http.server._url_collapse_path_split(self.path)
        if self.command == 'POST' and self.is_cgi() and splitpath[1] in upload_keywords:
            return True
        return False

    def is_getwu(self):
        """Test whether request is for a new WU."""
        if self.is_cgi():
            filename=self.cgi_info[1].split("?", 1)[0]
            if self.command == 'GET' and filename in ['getwu']:
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

    db = DbThread(dbfilename)

    HandlerClass = MyHandler
    httpd = ServerClass((HTTP, PORT), HandlerClass)
    httpd.server_name = "test"

    print ("serving at " + HTTP + ":" + str(PORT))
    httpd.serve_forever()
    db.finish()
