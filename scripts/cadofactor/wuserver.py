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
    # Overrides http.server.CGIHTTPRequestHandler.cgi_directories
    clientid_pattern = re.compile("^[\w.-]*$")
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.no_work_available = False
    
    def log(self, lvl, format, *args, **kwargs):
        """ Interface to the logger class. 
            We add the client address (as a string) to the log record so the 
            logger can print that """
        # e = kwargs.copy()
        # e["address_string"] = self.address_string()
        format = '%s ' + format
        self.logger.log(lvl, format, self.address_string(), *args, **kwargs)

    # These three methods overwrite the corresponding methods from 
    # http.server.BaseHTTPRequestHandler
    # They just call self.log() with a numerical logging level added
    def log_message(self, format, *args, **kwargs):
        self.log(logging.INFO, format, *args, **kwargs)

    def log_request(self, code='-', size='-'):
        self.log(logging.DEBUG, '"%s" %s %s', 
                 self.requestline, code, size)

    def log_error(self, format, *args):
        # Log errors with WARNING level, except messages about no work being
        # available, as those are frequent and kinda spammy
        level = logging.DEBUG if self.no_work_available else logging.WARNING
        self.log(level, format, *args)

    def send_body(self, body):
        self.wfile.write(body)
        self.wfile.flush()
    
    def translate_path(self, path):
        """ Translate path in request URL to local file system, taking into 
        account registered file names.
        Overrides SimpleHTTPRequestHandler.translate_path(); paths that are not
        in registered_filenames are delegated to super().translate_path()
        """
        # Path in url always starts with '/'
        relpath = self.path.lstrip('/')
        if relpath in self.registered_filenames:
            self.log(logging.DEBUG, "Translated path %s to %s", relpath, 
                     self.registered_filenames[relpath])
            return self.registered_filenames[relpath]
        self.log(logging.DEBUG, "Not translating path %s ", relpath)
        return super().translate_path(path)
    
    def do_GET(self):
        """Generates a work unit if request is cgi-bin/getwu, otherwise calls
           parent class' do_GET()"""
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
        if self.is_upload():
            super().do_POST()
        else:
            self.send_error(501, "POST request allowed only for uploads")
        sys.stdout.flush()

    def is_upload(self):
        """Test whether request is a file upload."""
        splitpath = urllib.parse.urlsplit(self.path)
        if self.command == 'POST' and self.is_cgi() and \
                splitpath.path == self.upload_path:
            return True
        return False

    def is_getwu(self):
        """Test whether request is for a new WU."""
        splitpath = urllib.parse.urlsplit(self.path)
        return self.command == 'GET' and splitpath.path in ['/getwu']

    def is_getstatus(self):
        """Test whether request is for a a status page."""
        splitpath = urllib.parse.urlsplit(self.path)
        return self.command == 'GET' and splitpath.path in ['/status']

    def guess_type(self, path):
        type = super().guess_type(path)
        # Use text/plain for files in upload, unless the type was properly 
        # identified
        # FIXME: make path identification more robust
        cwd = os.getcwd().rstrip(os.sep)
        if type == "application/octet-stream" and path.startswith(cwd + os.sep + 'upload' + os.sep):
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
        if not self.clientid_pattern.match(clientid):
            return self.send_error(400, "Malformed client id specified")
        
        if self.db_pool:
            wu_text = self.db_pool.assign(clientid)
        else:
            wu_text = wudb.WuAccess(self.dbfilename).assign(clientid)
        if not wu_text:
            # This flag is to downgrade the logging level. Ugly.
            self.no_work_available = True
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
        if self.db_pool:
            wus = self.db_pool.query(**conditions)
        else:
            wus = wudb.WuAccess(self.dbfilename).query(**conditions)

        body = HtmlGen()

        body.append('<a href="/index.html">Back to index</a>')
        body.append("<p>Query for conditions = " + str(conditions) + "</p>")

        if not wus is None and len(wus) > 0:
            cwd = os.getcwd()
            body.append(str(len(wus)) + " records match.")
            keys = wus[0].keys()
            body.start_table(keys)
            for wu in wus:
                body.wu_row(wu, keys, cwd)
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

class ServerLauncher(object):
    def __init__(self, address, port, threaded, dbfilename,
                registered_filenames, uploaddir, bg = False,
                use_db_pool = True, scriptdir = None):
        
        self.name = "HTTP server"
        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(logging.NOTSET)
        # formatter = logging.Formatter(
        #    fmt='%(address_string)s - - [%(asctime)s] %(message)s')
        #self.ch = logging.StreamHandler()
        #self.ch.setFormatter(formatter)
        #self.logger.addHandler(self.ch)
        #self.logger.propagate = False
        
        self.bg = bg
        if threaded:
            logging.info("Using threaded server")
            ServerClass = ThreadedHTTPServer
        else:
            logging.info("Not using threaded server")
            ServerClass = http.server.HTTPServer
        if use_db_pool:
            self.db_pool = wudb.DbThreadPool(dbfilename, 1)
        else:
            self.db_pool = None
        # Generate a class with parameters stored in class variables
        upload_url_path = "/cgi-bin/upload.py"
        handler_params = {
            "registered_filenames": registered_filenames,
            "logger": self.logger,
            "dbfilename": dbfilename,
            "db_pool": self.db_pool, 
            "uploaddir": uploaddir,
            "cgi_directories" : ['/cgi-bin'],
            "upload_path": upload_url_path
        }
        MyHandlerWithParams = type("MyHandlerWithParams", (MyHandler, ), handler_params)
        
        # Find the upload.py script
        upload_path = self.findscript("upload.py", scriptdir)
        # Always register the upload script
        if upload_path is None:
            raise IOError("upload.py script not found")
        self.logger.debug("Found upload.py at %s" % upload_path)
        registered_filenames[upload_url_path.lstrip('/')] = upload_path

        # Set shell environment variables which the upload.py script needs if
        # spawned as subprocess
        os.environ[upload.DBFILENAMEKEY] = dbfilename
        os.environ[upload.UPLOADDIRKEY] = uploaddir
        if not os.path.isdir(uploaddir):
            os.mkdir(uploaddir)
        
        self.httpd = ServerClass((address, port), MyHandlerWithParams, )
        self.httpd.server_name = self.name
    
    def serve(self):
        logging.info("serving at http://%s:%d",
                     *self.httpd.server_address)
        
        if self.bg:
            from threading import Thread
            self.thread = Thread(target=self.httpd.serve_forever,
                                 name=self.name)
            self.thread.daemon = True
            self.thread.start()
        else:
            self.httpd.serve_forever()
    
    def shutdown(self):
        logging.info("Shutting down HTTP server")
        self.httpd.shutdown()
        if self.bg:
            self.thread.join()
        if self.db_pool:
            self.db_pool.terminate()
        #self.logger.removeHandler(self.ch)
    
    @staticmethod
    def findscript(scriptname, scriptdir = None):
        # If scriptdir is specified, use that unconditionally
        if not scriptdir is None:
            return scriptdir + os.sep + scriptname
        # Try the CWD
        if os.path.isfile(scriptname):
            return scriptname
        # Try the directory that contains wuserver.py
        if not __file__ is None:
            dirname = os.path.dirname(os.path.realpath(__file__))
            if os.path.isfile(dirname + os.sep + scriptname):
                return dirname + os.sep + scriptname
        # Not found
        return None

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-address", help="Listen address", default="localhost")
    parser.add_argument("-port", help="Listen port", default="8001")
    parser.add_argument("-uploaddir", help="Upload directory", 
                        default="upload/")
    parser.add_argument("-dbfile", help="Database file name", required=True)
    parser.add_argument("-threaded", help="Use threaded server", 
                        action="store_true", default=False)
    args = parser.parse_args()

    PORT = int(args.port)
    HTTP = args.address
    dbfilename = args.dbfile
    registered_filenames = {}

    logger = logging.getLogger()
    logger.setLevel(logging.NOTSET)

    httpd = ServerLauncher(HTTP, PORT, args.threaded, dbfilename, 
                           registered_filenames, args.uploaddir)
    
    try:
        httpd.serve()
    except KeyboardInterrupt:
        pass
    else:
        raise
    httpd.shutdown()
