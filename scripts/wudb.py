#!/usr/bin/env python3

import sys
import sqlite3
from datetime import datetime
from workunit import Workunit

debug = 1

def diag(level, text, var = None):
    if debug > level:
        if var is None:
            print (text, file=sys.stderr)
        else:
            print (text + str(var), file=sys.stderr)

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
        """ Check whether status is equal to one of the constants """
        assert status in (WuStatus.AVAILABLE, WuStatus.ASSIGNED, WuStatus.RECEIVED_OK, 
            WuStatus.RECEIVED_ERROR, WuStatus.VERIFIED_OK, WuStatus.VERIFIED_ERROR, 
            WuStatus.CANCELLED)

# If we try to update the status in any way other than progressive 
# (AVAILABLE -> ASSIGNED -> ...), we raise this exception
class StatusUpdateError(Exception):
    pass

class WuDb: # {
    """ This class represents a DB connection and provides wrappers around 
        SQL queries to individual specified tables.
        That is, it acts as a Gateway between the SQL syntax and the table/row/column 
        structure of a relational database table on one hand, and the dictionary data 
        structure of Python on the other hand, where one row of one table maps to one 
        dictoriary """

    def __init__(self, filename):
        """ Open a connection to a sqlite database with the specified filename """
        self.db = sqlite3.connect(filename)

    def __del__(self):
        self.close()

    def close(self):
        self.db.close()

    @staticmethod
    def _fieldlist(l):
        """ For a list l = ('a', 'b', 'c') returns the string 'a=?, b=?, c=?' """
        return ", ".join([k + "=?" for k in l])

    @staticmethod
    def _without_None(d):
        """ Return a copy of the dictionary d, but without entries whose values 
            are None """
        return {k[0]:k[1] for k in d.items() if k[1] is not None}
    
    @staticmethod
    def _without_id(d):
        """ Return a copy of the dictionary d, but without the "id" entry """
        return {k[0]:k[1] for k in d.items() if k[0] != "id"}
    
    @staticmethod
    def _exec(cursor, command, values, name):
        """ Wrapper around self.cursor.execute() that prints arguments 
            for debugging """
        # Could use inspect module to remove name parameter
        diag (1, "WuDb." + name + "(): command = " + command);
        diag (1, "WuDb." + name + "(): values = ", values)
        cursor.execute(command, values)

    def create_table(self, table, layout):
        """ Creates a table with fields as described in the layout parameter """
        command = "CREATE TABLE IF NOT EXISTS " + table + \
            "( " + ", ".join([" ".join(col) for col in layout]) + " );"
        cursor = self.db.cursor()
        WuDb._exec (cursor, command, (), "create_table")
        self.db.commit()
        cursor.close()
    
    def insert(self, table, d):
        """ Insert a new entry, where d is a dictionary containing the 
            field:value pairs. Returns the row id of the newly created entry """
        # INSERT INTO table (field_1, field_2, ..., field_n) 
        # 	VALUES (value_1, value_2, ..., value_n)

        # Fields is a copy of d but with entries removed that have value None.
        # This is done primarily to avoid having "id" listed explicitly in the 
        # INSERT statement, because the DB fills in a new value automatically 
        # if "id" is the primary key. But I guess not listing field:NULL items 
        # explicitly in an INSERT is a good thing in general
        fields = self._without_None(d)

        sqlformat = ", ".join(("?",) * len(fields)) # sqlformat = "?, ?, ?, " ... "?"
        command = "INSERT INTO " + table + \
            " (" + ", ".join(fields.keys()) + ") VALUES (" + sqlformat + ");"
        values = list(fields.values())
        cursor = self.db.cursor()            
        WuDb._exec(cursor, command, values, "insert")
        id = cursor.lastrowid
        self.db.commit()
        cursor.close()
        return id

    def update(self, table, id, d):
        """ Update fields of an existing entry. id is the row id of the 
            entry to update, other entries in the dictionary d are the fields 
            and their values to update """
        # UPDATE table SET column_1=value1, column2=value_2, ..., 
        # column_n=value_n WHERE id="id"
        # FIXME: can generalize this a bit by passing a dict for the where clause
        command = "UPDATE " + table + " SET " + WuDb._fieldlist(d.keys()) + \
            " WHERE id=?;"
        values = list(d.values()) + [id, ]
        cursor = self.db.cursor()            
        WuDb._exec(cursor, command, values, "update")
        self.db.commit()
        cursor.close()
    
    def where_eq(self, table, d = None, limit=None):
        """ Get a up to "limit" table rows (limit == 0: no limit) where 
            the key:value pairs of the dictionary d are set to the same 
            value in the database table """
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
        
        # cursor.description is a list of lists, where the first element of 
        # each inner list is the column name
        desc = [k[0] for k in cursor.description]
        row = cursor.fetchone()
        while row is not None:
            diag(1, "WuDb.where_eq(): row = ", row)
            result.append(dict(zip(desc, row)))
            row = cursor.fetchone()
        cursor.close()
        return result
# }

class DbTable: # {
    """ A class template defining access methods to a database table """
    def __init__(self, db, tablename, fields):
        self.db = db
        self.tablename = tablename
        self.fields = fields

    @staticmethod
    def _subdict(d, l):
        """ Returns a dictionary of those key:value pairs of d for which key 
            exists l """
        if d is None:
            return None
        return {k:d[k] for k in d.keys() if k in l}

    def _get_colnames(self):
        return [k[0] for k in self.fields]

    def dictextract(self, d):
        """ Return a dictionary with all those key:value pairs of d
            for which key is in self._get_colnames() """
        return self._subdict(d, self._get_colnames())

    def create(self):
        db.create_table(self.tablename, self.fields)

    def insert(self, d):
        """ Insert a new row into this table. The column:value pairs are 
            specified key:value pairs of the dictionary d. 
            The database's row id for the new entry is returned """
        return self.db.insert(self.tablename, self.dictextract(d))

    def update(self, id, d):
        """ Update an existing row in this table. The column:value pairs to 
            be written are specified key:value pairs of the dictionary d """
        self.db.update(self.tablename, id, d)

    def where_eq(self, d = None, limit = None):
        return self.db.where_eq(self.tablename, d, limit)
# }

class WuTable(DbTable):
    name = "workunits"
    fields = (
        ("id", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"), 
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
    def __init__(self, db):
        DbTable.__init__(self, db, WuTable.name, WuTable.fields)


class FilesTable(DbTable):
    name = "files"
    fields = (
        ("id", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"), 
        ("wuid", "TEXT", ""), 
        ("filename", "TEXT", ""), 
        ("path", "TEXT", "UNIQUE NOT NULL")
    )
    def __init__(self, db):
        DbTable.__init__(self, db, FilesTable.name, FilesTable.fields)


class WuActiveRecord(): # {
    """ This class maps between the WORKUNIT and FILES tables 
        and a dictionary 
        {"wuid": string, ..., "timeverified": string, "files": list}
        where list is None or a list of dictionaries of the from
        {"id": int, "wuid": string, "filename1": string, "path": string
        Operations on instances of WuActiveRecord are directly carried 
        out on the database persistent storage, i.e., they behave kind 
        of as if the WuActiveRecord instance were itself a persistent 
        storage device """
    
    def __init__(self, db):
        self.wutable = WuTable(db)
        self.filestable = FilesTable(db)

    def __str__(self):
        s = "Workunit " + str(self.data["wuid"]) + ":\n"
        for (k,v) in self.data.items():
            if k != "wuid" and k != "files":
                s = s + "  " + k + ": " + repr(v) + "\n"
        if "files" in self.data:
            s = s + "  Associated files:\n"
            if self.data["files"] is None:
                s = s + "    None\n"
            else:
                for f in self.data["files"]:
                    s = s + "    " + str(f) + "\n"
        return s

    def create_tables(self):
        self.wutable.create()
        self.filestable.create()

    def add_files(self, files):
        if len(files) > 0 and self.data["files"] is None:
            self.data["files"] = []
        for f in files:
            d = {"filename": f[0], "path": f[1]}
            d["wuid"] = self.data["wuid"]
            d["id"] = self.filestable.insert(d)
            del d["wuid"]
            self.data["files"].append(d)

    def _from_db_row(self, wu_row, files_rows):
        """ Return a WuActiveRecord instance with the data from the workunit table row
            wu_row, and files table rows files_rows """
        # Create an instance
        wu = type(self)(self.wutable.db)
        # Fill in data from the workunit row
        wu.data = self.wutable.dictextract(wu_row)
        assert not "files" in wu.data
        if len(files_rows) > 0:
            wu.data["files"] = []
        else:
            wu.data["files"] = None
        for f in files_rows:
            fileentry = self.filestable.dictextract(f)
            assert fileentry["wuid"] == wu.data["wuid"]
            del(fileentry["wuid"])
            wu.data["files"].append(fileentry)
        return wu

    def where_eq(self, where = None, limit = None):
        result = []
        wu_d = self.wutable.dictextract(where)
        wu_rows = self.wutable.where_eq(wu_d, limit)
        # print ("* wu_rows = " + str(wu_rows))
        for wu_row in wu_rows:	
            # print ("* wu_row = " + str(wu_row))
            files_rows = self.filestable.where_eq({"wuid" : wu_row["wuid"]})
            wu = self._from_db_row(wu_row, files_rows)
            # print ("* wu = " + str(wu))
            result.append(wu)
        return result
        
    def update_wu(self, d):
        """ Assign the key:value pairs in d to self.data, and call 
            db.update() method to write these updates to the DB """
        self.data.update(d) # Python built-in dict.update() method
        self.wutable.update(self.data["id"], d)
    
    def get_wuid(self):
        return self.data["wuid"]
    
    def get_wu(self):
        return self.data["wu"]
    
    def find_by_wuid(self, wuid):
        r = self.where_eq({"wuid": wuid}, limit = 1)
        if len(r) == 0:
            raise Exception
        self.data = r[0].data

    def _checkstatus(self, status):
        diag (2, "WuActiveRecord._checkstatus(" + str(self) + ", " + str(status) + ")")
        if self.data["status"] != status:
            wuid = self.get_wuid()
            wustatus = str(self.data["status"])
            msg = "WU " + wuid + " has status " + wustatus + ", expected " + str(status)
            diag (0, "WuActiveRecord._checkstatus(): " + msg)
            # FIXME: this raise has no effect other than returning from the method. Why?
            # raise wudb.StatusUpdateError(msg)
            raise Exception(msg)

    def check(self):
        status = self.data["status"]
        WuStatus.check(status)
        wu = Workunit(self.get_wu())
        assert wu.get_id() == self.get_wuid()
        if status > WuStatus.RECEIVED_ERROR:
            return
        if status == WuStatus.RECEIVED_ERROR:
            assert self.data["errorcode"] != 0
            return
        if status == WuStatus.RECEIVED_OK:
            assert self.data["errorcode"] == 0
            return
        assert self.data["errorcode"] is None
        assert self.data["timeresult"] is None
        assert self.data["resultclient"] is None
        if status == WuStatus.ASSIGNED:
            return
        assert self.data["timeassigned"] is None
        assert self.data["assignedclient"] is None
        if status == WuStatus.AVAILABLE:
            return
        assert self.data["timecreated"] is None
        # etc.
    
    # Here come the application-visible functions that implement the 
    # "business logic": creating a new workunit from the text of a WU file,
    # assigning it to a client, receiving a result for the WU, marking it as
    # verified, or marking it as cancelled
    def create(self, wu):
        """ Create a new workunit from wu which contains the text of the 
            workunit file """
        self.data = {}
        self.data["wuid"] = Workunit(wu).get_id()
        self.data["wu"] = wu
        self.data["status"] = WuStatus.AVAILABLE
        self.data["timecreated"] = str(datetime.now())
        self.data["files"] = []
        self.data["id"] = self.wutable.insert(self.data)

    def assign(self, clientid):
        """ Finds an available workunit and assigns it to clientid.
            Returns None of no available workunit exists """
        r = self.where_eq({"status": WuStatus.AVAILABLE}, limit = 1)
        if len(r) == 0:
            return None
        self.data = r[0].data
        self._checkstatus(WuStatus.AVAILABLE)
        if debug > 0:
            self.check()
        d = {"status": WuStatus.ASSIGNED, 
             "assignedclient": clientid,
             "timeassigned": str(datetime.now())}
        self.update_wu(d)

    def result(self, wuid, clientid, errorcode, files):
        self.find_by_wuid(wuid)
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
        self.update_wu(d)
        self.add_files(files)

    def verification(self, wuid, ok):
        self.find_by_wuid(wuid)
        self._checkstatus(WuStatus.RECEIVED_OK)
        if debug > 0:
            self.check()
        d = {["timeverified"]: str(datetime.now())}
        if ok:
            d["status"] = WuStatus.VERIFIED_OK
        else:
            d["status"] = WuStatus.VERIFIED_ERROR
        self.update_wu(d)

    def cancel(self):
        pass
# }

# One entry in the WU DB, including the text with the WU contents 
# (FILEs, COMMANDs, etc.) and info about the progress on this WU (when and 
# to whom assigned, received, etc.)

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


        


if __name__ == '__main__':
    import sys
    dbname = "wudb"
    for arg in sys.argv:
        args = arg.split("=")
        if args[0] == "-dbname":
            dbname = args[1]
        if args[0] == "-create":
            db = WuDb(dbname)
            wu = WuActiveRecord(db)
            wu .create_tables()
            db.close()
        if args[0] == "-add":
            db = WuDb(dbname)
            wutext = sys.stdin.read()
            wu = WuActiveRecord(db)
            wu.create(wutext)
            db.close()
        # Functions for testing
        if args[0] == "-assign":
            db = WuDb(dbname)
            wu = WuActiveRecord(db)
            clientid = "client1"
            wu.assign(clientid)
            db.close()
        if args[0] == "-result":
            db = WuDb(dbname)
            wu = WuActiveRecord(db)
            wuid = "a"
            clientid = "client1"
            files = (("output", "/foo/bar"),)
            errorcode = 0
            wu.result(wuid, clientid, errorcode, files)
            db.close()
        # Functions for queries
        if args[0] == "-avail":
            db = WuDb(dbname)
            wu = WuActiveRecord(db)
            wus = wu.where_eq({"status": WuStatus.AVAILABLE})
            print("Available workunits: ")
            for wu in wus:
                print (str(wu))
            db.close()
        if args[0] == "-assigned":
            db = WuDb(dbname)
            wu = WuActiveRecord(db)
            wus = wu.where_eq({"status": WuStatus.ASSIGNED})
            print("Assigned workunits: ")
            for wu in wus:
                print (str(wu))
            db.close()
        if args[0] == "-receivedok":
            db = WuDb(dbname)
            wu = WuActiveRecord(db)
            wus = wu.where_eq({"status": WuStatus.RECEIVED_OK})
            print("Received ok workunits: ")
            for wu in wus:
                print (str(wu))
            db.close()
        if args[0] == "-receivederr":
            db = WuDb(dbname)
            wu = WuActiveRecord(db)
            wus = wu.where_eq({"status": WuStatus.RECEIVED_ERROR})
            print("Received with error workunits: ")
            for wu in wus:
                print (str(wu))
            db.close()
        if args[0] == "-all":
            db = WuDb(dbname)
            wu = WuActiveRecord(db)
            wus = wu.where_eq()
            print("Existing workunits: ")
            for wu in wus:
                print (str(wu))
            db.close()
