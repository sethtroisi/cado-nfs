#!/usr/bin/env python3

import http.server
import socketserver
import os
import sys
import re
import io
import logging
import urllib.parse
from workunit import Workunit
import datetime
import wudb
import upload

# Get the shell environment variable name in which we should store the path 
# to the upload directory

class ThreadedHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    """Handle requests in a separate thread."""

class HttpServerLogger(object):
    def __init__(self, level):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(level)
        formatter = logging.Formatter(
            fmt='%(address_string)s - - [%(asctime)s] %(message)s')
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
    
    def log(self, lvl, *args, **kwargs):
        self.logger.log(lvl, *args, **kwargs)

class HtmlGen(io.BytesIO):
    def __init__(self, encoding = None):
        super().__init__()
        if encoding is None:
            self.encoding = 'utf-8'
        else:
            self.encoding = encoding

    def header(self):
        self.write(
            b'<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" ' + 
            b'"http://www.w3.org/TR/html4/strict.dtd">\n' + 
            b'<html>\n' + 
            b'<head>\n' + 
            b'<meta http-equiv="content-type" content="text/html; ' + 
              b'charset=' + self.encoding.encode("ascii") + b'">\n' 
            b'<title>List of workunits</title>\n' + 
            b'</head>\n' + 
            b'<body>')

    def finish(self):
        self.write(b'</body>')

    def __bytes__(self):
        return self.getvalue()

    def get_len(self):
        return len(self.getvalue())

    def append(self, str):
        self.write(str.encode(self.encoding))

    def start_table(self, fields):
        self.append('<table border="1">\n<tr>')
        for h in fields:
            self.append('<th>' + h + '</th>')
        self.append('</tr>\n')

    def add_table_row(self, row):
        self.append('<tr>')
        for d in row:
            self.append('<td>' + str(d) + '</td>')
        self.append('</tr>\n')

    def end_table(self):
        self.append('</table>\n')

    def wu_row(self, wu, fields, cwd):
        arr = []
        for k in fields:
            if k == "files" and not wu["files"] is None:
                s = ""
                for f in wu["files"]:
                    path = f["path"]
                    if path.startswith(cwd):
                        path = path[len(cwd):]
                    s = s + '<a href="' + path + '">' + f["filename"] + \
                    '</a><br>'
                arr.append(s)
            else:
                arr.append(wu[k])
        self.add_table_row(arr)


class MyHandler(http.server.CGIHTTPRequestHandler):
    upload_keywords = ['/upload.py']

    def log(self, lvl, format, *args, **kwargs):
        """ Interface to the logger class. 
            We add the client address (as a string) to the log record so the 
            logger can print that """
        e = kwargs.copy()
        e["address_string"] = self.address_string()
        logger.log(lvl, format, *args, extra=e)

    # These three methods overwrite the corresponding methods from 
    # http.server.BaseHTTPRequestHandler
    # They just call self.log() with a numerical logging level added
    def log_message(self, format, *args, **kwargs):
        self.log(logging.INFO, format, *args, **kwargs)

    def log_request(self, code='-', size='-'):
        self.log(logging.INFO, '"%s" %s %s', self.requestline, str(code), 
                 str(size))

    def log_error(self, format, *args):
        self.log(logging.ERROR, format, *args)

    def send_body(self, body):
        self.wfile.write(body)
        self.wfile.flush()

    def do_GET(self):
        """Generates a work unit if request is cgi-bin/getwu, otherwise calls
           parent class' do_GET()"""
        self.cwd = os.getcwd()
        if self.is_cgi():
            if self.is_getwu():
                self.send_WU()
            elif self.is_getstatus():
                self.send_status()
            else:
                self.send_error(404, "GET for CGI scripts allowed only " + 
                                "for work unit or status page request")
        else:
            super().do_GET()
        sys.stdout.flush()

    def do_POST(self):
        """Set environment variable telling the upload directory 
           and call CGI handler to run upload CGI script"""
        self.cwd = os.getcwd()
        if self.is_upload():
            if False:
                self.send_response(200, "Script output follows")
                upload.do_upload(self.dbfilename, input = self.rfile, output = self.wfile)
            else:
                http.server.CGIHTTPRequestHandler.do_POST(self)
        else:
            self.send_error(404, "POST request allowed only for uploads")
        sys.stdout.flush()

    def is_upload(self):
        """Test whether request is a file upload."""
        splitpath = urllib.parse.urlsplit(self.path)
        if self.command == 'POST' and self.is_cgi() and \
                splitpath.path in self.upload_keywords:
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

    def guess_type(self, path):
        type = super().guess_type(path)
        # Use text/plain for files in upload, unless the type was properly 
        # identified
        # FIXME: make path identification more robust
        if type == "application/octet-stream" and path.startswith(self.cwd + '/upload/'):
            return "text/plain"
        return type

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
        
        # wu = wudb.WuAccess(self.dbfilename)
        wu_text = db_pool.assign(clientid)
        if not wu_text:
            return self.send_error(404, "No work available")
        
        self.log_message("Sending work unit " + Workunit(wu_text).get_id() + 
                         " to client " + clientid)
        # wu_text = wu.get_wu()
        self.send_response(200)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", len(wu_text))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        # FIXME: is ASCII enough for workunits? Is there any shell syntax
        # that needs more, or should we allow non-ASCII workunit names?
        self.send_body(bytes(wu_text, "ascii"))

    def send_status(self):
        self.send_query()

    def send_query(self):
        logging.debug("self.cgi_info = "  + str(self.cgi_info))
        filename = self.cgi_info[1]
        if "#" in filename:
            # Get rid of fragment part
            filename = filename.split("#", 1)[0]
        conditions = {}
        if "?" in filename:
            # Parse query part into SELECT conditions
            (filename, query) = filename.split("?", 1)
            print("Query = " + query)
            conditions = {}
            # Now look at individual key=value pairs
            for q in query.split("&"):
                q = urllib.parse.unquote_plus(q)
                logging.debug("Processing token " + str(q))
                for (name, op) in wudb.MyCursor.name_to_operator.items():
                    if op in q:
                        (key, value) = q.split(op, 1)
                        if not name in conditions:
                            conditions[name] = {}
                        # If value is of the form "now(-123)", convert it to a 
                        # time stamp of 123 minutes ago
                        r = re.match(r"now\((-?\d+)\)", value)
                        if r:
                            minutes_ago = int(r.group(1))
                            td = datetime.timedelta(minutes = minutes_ago)
                            value = str(datetime.datetime.now() + td)
                        conditions[name][key] = value
                        break
        wus = db_pool.query(**conditions)

        body = HtmlGen()

        body.append('<a href="/index.html">Back to index</a>')
        body.append("<p>Query for conditions = " + str(conditions) + "</p>")

        if not wus is None and len(wus) > 0:
            body.append(str(len(wus)) + " records match.")
            keys = wus[0].keys()
            body.start_table(keys)
            for wu in wus:
                body.wu_row(wu, keys, self.cwd)
            body.end_table()
        else:
            body.append("No records match.")
        body.finish()
        
        self.send_response(200)
        self.send_header("Content-Type", "text/html")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Content-Length", body.get_len())
        self.end_headers()
        self.send_body(body.__bytes__())

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-address", help="Listen address", default="localhost")
    parser.add_argument("-port", help="Listen port", default="8001")
    parser.add_argument("-uploaddir", help="Upload directory", default="upload/")
    parser.add_argument("-dbfile", help="Database file name", required=True)
    parser.add_argument("-threaded", help="Use threaded server", action="store_true", default=False)
    args = parser.parse_args()

    if args.threaded:
        print("Using threaded server")
        ServerClass = ThreadedHTTPServer
    else:
        print("Not using threaded server")
        ServerClass = http.server.HTTPServer

    logger = HttpServerLogger(logging.INFO)

    PORT = int(args.port)
    HTTP = args.address
    dbfilename = args.dbfile
    # Set shell environment variables which the upload.py script spawned as 
    # subprocess needs
    os.environ[upload.DBFILENAMEKEY] = dbfilename
    os.environ[upload.UPLOADDIRKEY] = args.uploaddir

    if not os.path.isdir(args.uploaddir):
        os.mkdir(args.uploaddir)

    db_pool = wudb.DbThreadPool(dbfilename, 1)

    HandlerClass = MyHandler
    HandlerClass.dbfilename = dbfilename

    httpd = ServerClass((HTTP, PORT), HandlerClass)
    httpd.server_name = "test"

    print ("serving at " + HTTP + ":" + str(PORT))
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        db_pool.terminate()
    else:
        raise
