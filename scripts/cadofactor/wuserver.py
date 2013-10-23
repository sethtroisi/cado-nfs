#!/usr/bin/env python3

import http.server
import socketserver
from socket import error as socket_error
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
import select
import errno

# Import upload to get the shell environment variable name in which we should
# store the path to the upload directory

class FixedHTTPServer(http.server.HTTPServer):
    """ Work-around class for http.server.HTTPServer that handles EINTR """
    def serve_forever(self, *args, **kwargs):
        """ Wrapper around http.server.HTTPServer.serve_forever() that
        restarts in case of EINTR.

        See http://bugs.python.org/issue7978
        """
        while True:
            try:
                return super().serve_forever(*args, **kwargs)
            except OSError as e:
                if e.errno != errno.EINTR:
                    raise
            except select.error as e:
                if e.args[0] != errno.EINTR:
                    raise


class ThreadedHTTPServer(socketserver.ThreadingMixIn, FixedHTTPServer):
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
    clientid_pattern = re.compile("^clientid=([\w.-]*)$")
    
    def __init__(self, *args, **kwargs):
        self.no_work_available = False
        super().__init__(*args, **kwargs)
    
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
        self.log(logging.DEBUG, '"%s" %s %s', self.requestline, code, size)

    def log_error(self, format, *args, **kwargs):
        # Log errors with WARNING level, except messages about no work being
        # available, as those are frequent and kinda spammy
        level = logging.DEBUG if self.no_work_available else logging.WARNING
        self.log(level, format, *args, **kwargs)

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
        # print("translate_path(%s)" % path)
        relpath = self.path.lstrip('/')
        if relpath in self.registered_filenames:
            self.log(logging.DEBUG, "Translated path %s to %s", self.path, 
                     self.registered_filenames[relpath])
            return self.registered_filenames[relpath]
        self.log(logging.DEBUG, "Not translating path %s ", self.path)
        if self.only_registered:
            return None
        else:
            return super().translate_path(path)
    
    def do_GET(self):
        """Generates a workunit if request is cgi-bin/getwu, generates a status
        page is requested, otherwise calls parent class' do_GET()"""
        if self.is_cgi():
            if self.is_getwu():
                self.send_WU()
            elif self.is_getstatus():
                self.send_status()
            else:
                self.send_error(404, "GET for CGI scripts allowed only "
                                "for workunit or status page request")
        elif self.only_registered and \
                not self.path.lstrip('/') in self.registered_filenames:
                self.send_error(404, "Access restricted to registered file "
                                "names, %s is not registered" % self.path)
        else:
            try:
                super().do_GET()
            except socket_error as e:
                self.log_error("Connection error: %s", str(e))
        sys.stdout.flush()

    def do_POST(self):
        """Set environment variable telling the upload directory 
           and call CGI handler to run upload CGI script"""
        if self.is_upload():
            super().do_POST()
        else:
            self.send_error(501, "POST request allowed only for uploads")
        sys.stdout.flush()
        sys.stderr.flush()

    def is_upload(self):
        """Test whether request is a file upload."""
        splitpath = urllib.parse.urlsplit(self.path)
        if self.command == 'POST' and self.is_cgi() and \
                splitpath.path == self.upload_path:
            return True
        return False

    def is_getwu(self):
        """Test whether request is for a new WU."""
        if not self.command == 'GET':
            return False
        splitpath = urllib.parse.urlsplit(self.path)
        # print("is_getwu(): path = %s, splitpath.path = %s" %
        #       (self.path, splitpath.path))
        return splitpath.path in [directory + '/getwu'
                                  for directory in self.cgi_directories]

    def is_getstatus(self):
        """Test whether request is for a a status page."""
        if not self.command == 'GET':
            return False
        splitpath = urllib.parse.urlsplit(self.path)
        # print("is_getwu(): path = %s, splitpath.path = %s" %
        #       (self.path, splitpath.path))
        return splitpath.path in [directory + '/status'
                                  for directory in self.cgi_directories]

    def guess_type(self, path):
        type = super().guess_type(path)
        # Use text/plain for files in upload, unless the type was properly 
        # identified
        # FIXME: testing the CWD here is wrong. If we want to expose files in
        # the upload directory, we either need to register all of them, or
        # register the upload directory and add directory name translation.
        cwd = os.getcwd().rstrip(os.sep)
        if type == "application/octet-stream" and \
                 path.startswith(cwd + os.sep + 'upload' + os.sep):
            return "text/plain"
        return type

    def send_WU(self):
        splitpath = urllib.parse.urlsplit(self.path)
        # print(str(splitpath))
        clientid_match = self.clientid_pattern.match(splitpath.query)
        if not clientid_match:
            return self.send_error(400, "No or malformed client id specified")
        clientid = clientid_match.group(1)
        
        if self.db_pool:
            wu_text = self.db_pool.assign(clientid)
        else:
            wu_text = wudb.WuAccess(self.dbfilename).assign(clientid)
        if not wu_text:
            # This flag is to downgrade the logging level. Ugly.
            self.no_work_available = True
            try:
                return self.send_error(404, "No work available")
            except socket_error as e:
                self.log_error("Connection error: %s", str(e))
                return
        
        self.log_message("Sending workunit " + Workunit(wu_text).get_id() + 
                         " to client " + clientid)
        # wu_text = wu.get_wu()
        try:
            self.send_response(200)
            self.send_header("Content-Type", "text/plain")
            self.send_header("Content-Length", len(wu_text))
            self.send_header("Cache-Control", "no-cache")
            self.end_headers()
            # FIXME: is ASCII enough for workunits? Is there any shell syntax
            # that needs more, or should we allow non-ASCII workunit names?
            self.send_body(bytes(wu_text, "ascii"))
        except socket_error as e:
            self.log_error("Connection error: %s", str(e))

    def send_status(self):
        self.send_query()
    
    def send_query(self):
        splitpath = urllib.parse.urlsplit(self.path)
        conditions = {}
        logging.debug("send_query(): query: %s", splitpath.query)
        if splitpath.query:
            # Parse query part into SELECT conditions
            # Now look at individual key=value pairs
            for query in splitpath.query.split("&"):
                query = urllib.parse.unquote_plus(query)
                # logging.debug("Processing token " + str(query))
                for (name, op) in wudb.MyCursor.name_to_operator.items():
                    if op in query:
                        (key, value) = query.split(op, 1)
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
        logging.debug("send_query(): conditions: %s", conditions)
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
                registered_filenames, uploaddir, *, bg = False,
                use_db_pool = True, scriptdir = None, only_registered=False):
        
        self.name = "HTTP server"
        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(logging.NOTSET)
        self.address = address
        self.port = port
        upload_scriptname = "upload.py"
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
            ServerClass = FixedHTTPServer
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
            "upload_path": upload_url_path,
            "only_registered": only_registered
        }
        MyHandlerWithParams = type("MyHandlerWithParams", (MyHandler, ), handler_params)
        
        # Find the upload.py script
        upload_path = self.findscript(upload_scriptname, scriptdir)
        # Always register the upload script
        if upload_path is None:
            raise IOError("%s script not found" % upload_scriptname)
        self.logger.debug("Found %s at %s" % (upload_scriptname, upload_path))
        registered_filenames[upload_url_path.lstrip('/')] = upload_path

        # Set shell environment variables which the upload.py script needs if
        # spawned as subprocess
        os.environ[upload.DBFILENAMEKEY] = dbfilename
        os.environ[upload.UPLOADDIRKEY] = uploaddir
        if not os.path.isdir(uploaddir):
            os.mkdir(uploaddir)
        
        try:
            self.httpd = ServerClass((address, port), MyHandlerWithParams, )
        except socket_error as e:
            if e.errno == errno.EADDRINUSE:
                self.logger.critical("Address %s:%d is already in use (maybe "
                        "another cadofactor running?)", address, port)
                self.logger.critical("You can choose a different port with "
                        "server.port=<integer>.")
                sys.exit(1)
        self.httpd.server_name = self.name
    
    def serve(self):
        logging.info("serving at http://%s:%d (%s)",
                     self.address, self.port, self.httpd.server_address[0])
        
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
    parser.add_argument("-onlyreg", help="Allow access only to registered files", 
                        action="store_true", default=False)
    args = parser.parse_args()

    PORT = int(args.port)
    HTTP = args.address
    dbfilename = args.dbfile
    registered_filenames = {}

    logger = logging.getLogger()
    logger.setLevel(logging.NOTSET)

    httpd = ServerLauncher(HTTP, PORT, args.threaded, dbfilename, 
                           registered_filenames, args.uploaddir,
                           only_registered=args.onlyreg)
    
    try:
        httpd.serve()
    except KeyboardInterrupt:
        pass
    else:
        raise
    httpd.shutdown()
