#!/usr/bin/env python3

import sys
import sqlite3
from datetime import datetime
from workunit import Workunit

debug = 1

def diag(level, text, var = None):
    if debug > level:
        if var == None:
            print (text, file=sys.stderr)
        else:
            print (text + str(var), file=sys.stderr)

# A DB entry for a WU uses the WUid as the key, and stores:
# - status (available, assigned, received, received with error, verified:correct, verified:error
# - Time created
# - Content of WU file, as dictionary of WU keys:values
# - Time assigned
# - Client-id assigned
# - Time result was received
# - Error code for result
# - Filenames uploaded with result
# - Time verified

# Dummy class for defining "constants"
class WuStatus:
    AVAILABLE = 0
    ASSIGNED = 1
    RECEIVED_OK = 2
    RECEIVED_ERROR = 3
    VERIFIED_OK = 4
    VERIFIED_ERROR = 5
    CANCELLED = 6

    def check(status):
        assert status in (WuStatus.AVAILABLE, WuStatus.ASSIGNED, WuStatus.RECEIVED_OK, 
            WuStatus.RECEIVED_ERROR, WuStatus.VERIFIED_OK, WuStatus.VERIFIED_ERROR, 
            WuStatus.CANCELLED)

# If we try to update the status in any way other than progressive 
# (AVAILABLE -> ASSIGNED -> ...), we raise this exception
class StatusUpdateError(Exception):
    pass

# This class represents a DB connection and provides wrappers around SQL queries
class WuDb:
    def __init__(self, filename):
        """ Open a connection to a sqlite database with the specified filename 
            and create a cursor. We also create the required tables if they do 
            not exits """
        self.db = sqlite3.connect(filename)

    def __del__(self):
        self.db.close()

    def _fieldlist(l):
        """ For a list ('a', 'b', 'c') returns 'a=?, b=?, c=?' """
        return ", ".join([k + "=?" for k in l])

    def _exec(cursor, command, values, name):
        """ Wrapper around self.cursor.execute() that prints arguments 
            for debugging """
        # Could use inspect module to remove name parameter
        diag (1, "WuDb." + name + "(): command = " + command);
        diag (1, "WuDb." + name + "(): values = ", values)
        cursor.execute(command, values)

    def create_table(self, table, layout):
        command = "CREATE TABLE IF NOT EXISTS " + table + \
            "( id INTEGER PRIMARY KEY ASC, " + \
            ", ".join([" ".join(col) for col in layout]) + " );"
        cursor = self.db.cursor()            
        WuDb._exec (cursor, command, (), "create_table")
        self.db.commit()
        cursor.close()
    
    def insert(self, table, d):
        """ Insert a new entry, where d is a dictionary containing the 
            field:value pairs. Returns the id of the newly created entry """
        # INSERT INTO WORKUNITS (field_1, field_2, ..., field_n) 
        # 	VALUES (value_1, value_2, ..., value_n)
        sqlformat = ", ".join(("?",) * len(d)) # sqlformat = "?, ?, ?, " ... "?"
        command = "INSERT INTO " + table + \
            " (" + ", ".join(d.keys()) + ") VALUES (" + sqlformat + ");"
        values = list(d.values())
        cursor = self.db.cursor()            
        WuDb._exec(cursor, command, values, "insert")
        id = cursor.lastrowid
        self.db.commit()
        cursor.close()
        return id

    def update(self, table, id, d):
        """ Update fields of an existing entry. id is the row id of the 
            entry to update, d is the dictionary of fields and their values 
            to update """
        # UPDATE WORKUNITS SET column_1=value1, column2=value_2, ..., 
        # column_n=value_n WHERE id="id"
        assert id is not None
        command = "UPDATE " + table + " SET " + WuDb._fieldlist(d.keys()) + \
            " WHERE id=?;"
        values = list(d.values()) + [id, ]
        cursor = self.db.cursor()            
        WuDb._exec(cursor, command, values, "update")
        self.db.commit()
        cursor.close()
    
    def read(self, table, id):
        """ Read an entry by id, return None if it does not exist """
        assert id is not None
        command = "SELECT * FROM " + table + " WHERE id=?;"
        values = (id,)
        cursor = self.db.cursor()            
        WuDb._exec(cursor, command, values, "read")
        row = cursor.fetchone()
        cursor.close()
        return row
    
    def where_eq(self, table, d = None, limit=None):
        """ Get a up to "limit" table rows (limit == 0: no limit) where 
            the key:value pairs of the dictionary d are set to the same 
            value in the database table "table" """
        result = []

        # Table/Column names cannot be substituted, so include in query directly.
        if d is None or len(d) == 0:
            WHERE = ""
            values = ()
        else:
            WHERE = " WHERE " + WuDb._fieldlist(d.keys())
            values = list(d.values())

        if limit is None:
            LIMIT = ""
        else:
            LIMIT = " LIMIT " + str(int(limit))
        command = "SELECT * FROM " + table + WHERE + LIMIT + ";"
        cursor = self.db.cursor()            
        WuDb._exec(cursor, command, values, "where_eq")
        
        desc = [k[0] for k in cursor.description]
        row = cursor.fetchone()
        while row is not None:
            diag(1, "WuDb.where_eq(): row = ", row)
            result.append(dict(zip(desc, row)))
            row = cursor.fetchone()
        cursor.close()
        return result


# One entry in the WU DB, including the text with the WU contents 
# (FILEs, COMMANDs, etc.) and info about the progress on this WU (when and 
# to whom assigned, received, etc.)

class DbWuEntry:
    # wuid is the unique wuid of the workunit
    # status is a status code as defined in WuStatus
    # data is the str containing the text of the workunit
    # timecreated is the string containing the date and time of when the WU was added to the db
    # timeassigned is the ... of when the WU was assigned to a client
    # assignedclient is the clientid of the client to which the WU was assigned
    # timeresult is the ... of when a result for this WU was received
    # resultclient is the clientid of the client that uploaded a result for this WU
    # errorcode is the exit status code of the first failed command, or 0 if none failed
    # timeverified is the ... of when the result was marked as verified
    table = "workunits"
    # Each entry specifies a sqlite column-def: column-name, type-name, column-constraint
    fields = (
        ("wuid", "TEXT", "UNIQUE NOT NULL"), 
        ("status", "INTEGER", "NOT NULL"), 
        ("wu", "TEXT", "NOT NULL"), 
        ("timecreated", "TEXT", ""), 
        ("timeassigned", "TEXT", ""), 
        ("assignedclient", "TEXT", ""), 
        ("timeresult", "TEXT", ""), 
        ("resultclient", "TEXT", ""), 
        ("errorcode", "INTEGER", ""), 
        ("timeverified", "TEXT", "")
    )

    _is_created = False

    def __init__(self, db):
        self.db = db
        self.id = None # The unique DB key (an integer) for this table row
        self.data = {key[0]: None for key in DbWuEntry.fields}
        if not DbWuEntry._is_created:
            self.db.create_table(DbWuEntry.table, DbWuEntry.fields)
            DbWuEntry._is_created = True
    
    def __str__(self):
        return str(self.data)

    def update(self, d):
        """ Assign the key:value pairs in d to self.data, and call 
            db.update() method to write these updates to the DB """
        self.data.update(d) # Python built-in dict.update() method
        self.db.update(DbWuEntry.table, self.id, d)

    def from_dict(self, d):
        self.id = d["id"]
        for k in DbWuEntry.fields:
            self.data[k[0]] = d[k[0]]
    
    def check(self):
        status = self.data["status"]
        WuStatus.check(status)
        wu = Workunit(self.data["wu"])
        assert wu.get_id() == self.data["wuid"]
        if status == WuStatus.AVAILABLE:
            assert self.data["timeassigned"] == None
            assert self.data["assignedclient"] == None
            return
        if status == WuStatus.ASSIGNED:
            assert self.data["timeresult"] == None
            assert self.data["resultclient"] == None
            return
        if status == WuStatus.RECEIVED_OK:
            assert self.data["errorcode"] == 0
        # etc.
    
    def get_wuid(self):
        return self.data["wuid"]
    
    def get_wu(self):
        return self.data["wu"]
    
    def _checkstatus(self, status):
        diag (1, "DbWuEntry._checkstatus(" + str(self) + ", " + str(status) + ")")
        if self.data["status"] != status:
            wuid = str(self.data["wuid"])
            wustatus = str(self.data["status"])
            msg = "WU " + wuid + " has status " + wustatus + ", expected " + str(status)
            diag (0, "DbWuEntry._checkstatus(): " + msg)
            # raise wudb.StatusUpdateError(msg)
            raise Exception(msg)

    def find_available(self):
        row = self.db.where_eq(DbWuEntry.table, {"status" : WuStatus.AVAILABLE}, limit=1)
        if len(row) == 0:
            return
        self.from_dict(row[0])

    def find_wuid(self, wuid):
        """ Find a row by its wuid """
        row = self.db.where_eq(DbWuEntry.table, {"wuid" : wuid}, limit=1)
        if len(row) == 0:
            return
        self.from_dict(row[0])

    def add(self, wu):
        self.data["wuid"] = Workunit(wu).get_id()
        self.data["wu"] = wu
        self.data["status"] = WuStatus.AVAILABLE
        self.data["timecreated"] = str(datetime.now())
        self.id = self.db.insert(DbWuEntry.table, self.data)

    def assign(self, clientid):
        self._checkstatus(WuStatus.AVAILABLE)
        if debug > 0:
            self.check()
        d = {"status": WuStatus.ASSIGNED, 
             "assignedclient": clientid,
             "timeassigned": str(datetime.now())}
        self.update(d)

    def result(self, clientid, errorcode, files):
        self._checkstatus(WuStatus.ASSIGNED)
        if debug > 0:
            self.check()
        d = {"resultclient": clientid,
             "errorcode": errorcode,
             "timeresult": str(datetime.now())}
        if errorcode == 0:
           d["status"] = WuStatus.RECEIVED_OK
        else:
            d["status"] = WuStatus.RECEIVED_ERROR
        self.update(d)

    def verification(self, ok):
        self._checkstatus(WuStatus.RECEIVED_OK)
        if debug > 0:
            self.check()
        d = {["timeverified"]: str(datetime.now())}
        if ok:
            d["status"] = WuStatus.VERIFIED_OK
        else:
            d["status"] = WuStatus.VERIFIED_ERROR
        self.update(d)


# A table storing names of uploaded files that belong to work units.
# The table stored that wuid, the filename under which the client 
# created and uploaded the file, and the file path and name under which 
# the server stored it.
class DbFilesEntry:
    
    table = "files"
    fields = (
        ("wuid", "TEXT", ""), 
        ("filename", "TEXT", ""), 
        ("path", "TEXT", "")
    )
    
    _is_created = False

    def __init__(self, db):
        self.db = db
        if not DbFilesEntry._is_created:
            self.db.create_table(DbFilesEntry.table, DbFilesEntry.fields)
            DbFilesEntry._is_created = True
        

if __name__ == '__main__':
    import sys
    dbname = "wudb"
    for arg in sys.argv:
        args = arg.split("=")
        if args[0] == "-dbname":
            dbname = args[1]
        if args[0] == "-add":
            db = WuDb(dbname)
            wutext = sys.stdin.read()
            wu = DbWuEntry(db)
            wu.add(wutext)
        if args[0] == "-avail":
            db = WuDb(dbname)
            available = db.where_eq("workunits", {"status": WuStatus.AVAILABLE})
            print("Available workunits: ")
            for wu in available:
                print (str(wu))
        if args[0] == "-assigned":
            db = WuDb(dbname)
            assigned = db.where_eq("workunits", {"status": WuStatus.ASSIGNED})
            print("Assigned workunits: ")
            for wu in available:
                print (str(wu))
        if args[0] == "-all":
            db = WuDb(dbname)
            all = db.where_eq("workunits")
            print("Existing workunits: ")
            for wu in all:
                print (str(wu))
        if args[0] == "-find_avail":
            db = WuDb(dbname)
            wu = DbWuEntry(db)
            wu.find_available()
            print("One aviailable workunit: ")
            print (str(wu))
