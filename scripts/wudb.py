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

    # This is used in where queries; it converts from named arguments such as 
    # "eq" to a binary operator such as "="
    name_to_operator = {"lt": "<", "le": "<=", "eq": "=", "ge": ">=", "gt" : ">", "ne": "!="}
    
    def __init__(self, filename):
        """ Open a connection to a sqlite database with the specified filename """
        # DEFERRED causes sqlite to create a SHARED lock after a read access, 
        # which (hopefully) prevents, e.g., race conditions between two threads
        # looking up an available workunit and assigning it to a client
        self.db = sqlite3.connect(filename, isolation_level="DEFERRED")

    def __del__(self):
        self.close()

    def cursor(self):
        return self.db.cursor()

    def commit(self):
        return self.db.commit()

    def close(self):
        self.db.close()

    @staticmethod
    def _fieldlist(l, r = "="):
        """ For a list l = ('a', 'b', 'c') returns the string 'a = ?, b = ?, c = ?',
            or with a different string r in place of the "=" """
        return ", ".join([k + " " + r + " ?" for k in l])

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

    @classmethod
    def where_str(cls, name, **args):
        where = ""
        values = []
        for opname in args:
            if args[opname] is None:
                continue
            if where == "":
                where = " " + name + " "
            else:
                where = where + ", "
            where = where + cls._fieldlist(args[opname].keys(), cls.name_to_operator[opname])
            values = values + list(args[opname].values())
        return (where, values)

    def create_table(self, cursor, table, layout):
        """ Creates a table with fields as described in the layout parameter """
        command = "CREATE TABLE IF NOT EXISTS " + table + \
            "( " + ", ".join([" ".join(col) for col in layout]) + " );"
        self.__class__._exec (cursor, command, (), "create_table")
    
    def insert(self, cursor, table, d):
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
        self.__class__._exec(cursor, command, values, "insert")
        id = cursor.lastrowid
        return id

    def update(self, cursor, table, d, **conditions):
        """ Update fields of an existing entry. conditions specifies the where 
            clause to use for to update, entries in the dictionary d are the 
            fields and their values to update """
        # UPDATE table SET column_1=value1, column2=value_2, ..., 
        # column_n=value_n WHERE column_n+1=value_n+1, ...,
        setstr = " SET " + self.__class__._fieldlist(d.keys())
        setvalues = d.values()
        (wherestr, wherevalues) = self.__class__.where_str("WHERE", **conditions)
        command = "UPDATE " + table + setstr + wherestr
        values = list(setvalues) + wherevalues
        self.__class__._exec(cursor, command, values, "update")
    
    def where(self, cursor, table, limit = None, order = None, **conditions):
        """ Get a up to "limit" table rows (limit == 0: no limit) where 
            the key:value pairs of the dictionary d are set to the same 
            value in the database table """
        result = []

        # Table/Column names cannot be substituted, so include in query directly.
        (WHERE, values) = self.__class__.where_str("WHERE", **conditions)

        if order is None:
            ORDER = ""
        else:
            if not order[1] in ("ASC", "DESC"):
                raise Exception
            ORDER = " ORDER BY " + str(order[0]) + " " + order[1]

        if limit is None:
            LIMIT = ""
        else:
            LIMIT = " LIMIT " + str(int(limit))

        command = "SELECT * FROM " + table + WHERE + ORDER + LIMIT + ";"
        self.__class__._exec(cursor, command, values, "where")
        
        # cursor.description is a list of lists, where the first element of 
        # each inner list is the column name
        desc = [k[0] for k in cursor.description]
        row = cursor.fetchone()
        while row is not None:
            diag(1, "WuDb.where(): row = ", row)
            result.append(dict(zip(desc, row)))
            row = cursor.fetchone()
        return result
# }

class DbTable: # {
    """ A class template defining access methods to a database table """
    def __init__(self, db):
        self.db = db
        self.tablename = type(self).name
        self.fields = type(self).fields

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

    def create(self, cursor):
        db.create_table(cursor, self.tablename, self.fields)

    def insert(self, cursor, d):
        """ Insert a new row into this table. The column:value pairs are 
            specified key:value pairs of the dictionary d. 
            The database's row id for the new entry is returned """
        return self.db.insert(cursor, self.tablename, self.dictextract(d))

    def update(self, cursor, d, **conditions):
        """ Update an existing row in this table. The column:value pairs to 
            be written are specified key:value pairs of the dictionary d """
        self.db.update(cursor, self.tablename, d, **conditions)

    def where(self, cursor, limit = None, order = None, **conditions):
        assert order is None or order[0] in self._get_colnames()
        return self.db.where(cursor, self.tablename, limit=limit, order=order, **conditions)
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
        ("failedcommand", "INTEGER", ""), 
        ("timeverified", "TEXT", ""),
        ("retryof", "TEXT", ""),
        ("priority", "INTEGER", "")
    )

class FilesTable(DbTable):
    name = "files"
    fields = (
        ("id", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"), 
        ("wurowid", "INTEGER", "REFERENCES " + WuTable.name + " (id)"), 
        ("filename", "TEXT", ""), 
        ("path", "TEXT", "UNIQUE NOT NULL")
    )

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
        self.db = db
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

    def add_files(self, cursor, files):
        if len(files) > 0 and self.data["files"] is None:
            self.data["files"] = []
        for f in files:
            d = {"filename": f[0], "path": f[1]}
            d["wurowid"] = self.data["id"]
            d["id"] = self.filestable.insert(cursor, d)
            del d["wurowid"]
            self.data["files"].append(d)

    def _from_db_row(self, wu_row, files_rows):
        """ Return a WuActiveRecord instance with the data from the workunit table row
            wu_row, and files table rows files_rows """
        # Create an instance
        wu = type(self)(self.db)
        # Fill in data from the workunit row
        wu.data = self.wutable.dictextract(wu_row)
        assert not "files" in wu.data
        if len(files_rows) > 0:
            wu.data["files"] = []
        else:
            wu.data["files"] = None
        for f in files_rows:
            fileentry = self.filestable.dictextract(f)
            assert fileentry["wurowid"] == wu.data["id"]
            del(fileentry["wurowid"])
            wu.data["files"].append(fileentry)
        return wu

    def where(self, cursor, limit = None, order = None, **conditions):
        result = []
        wu_rows = self.wutable.where(cursor, limit=limit, order=order, **conditions)
        # print ("* wu_rows = " + str(wu_rows))
        for wu_row in wu_rows:	
            # print ("* wu_row = " + str(wu_row))
            files_rows = self.filestable.where(cursor, eq={"wurowid" : wu_row["id"]})
            wu = self._from_db_row(wu_row, files_rows)
            # print ("* wu = " + str(wu))
            result.append(wu)
        return result
    
    def get_by_wuid(self, cursor, wuid):
        r = self.where(cursor, eq={"wuid": wuid}, limit = 1)
        if len(r) == 0:
            return False
        self.data = r[0].data
        return True

    def update_wu(self, cursor, d):
        """ Assign the key:value pairs in d to self.data, and call 
            db.update() method to write these updates to the DB """
        self.data.update(d) # Python built-in dict.update() method
        self.wutable.update(cursor, d, eq={"id": self.data["id"]})
    
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

    def create_tables(self):
        cursor = self.db.cursor()
        self.wutable.create(cursor)
        self.filestable.create(cursor)
        self.db.commit()
        cursor.close()

    def get_wuid(self):
        return self.data["wuid"]
    
    def get_wu(self):
        return self.data["wu"]
    
    def create1(self, cursor, wu, priority = None):
        self.data = {}
        self.data["wuid"] = Workunit(wu).get_id()
        self.data["wu"] = wu
        self.data["status"] = WuStatus.AVAILABLE
        self.data["timecreated"] = str(datetime.now())
        self.data["files"] = []
        if not priority is None:
            self.data["priority"] = priority
        self.data["id"] = self.wutable.insert(cursor, self.data)

    def create(self, wus, priority = None):
        """ Create a new workunit from wu which contains the text of the 
            workunit file """
        cursor = self.db.cursor()
        if isinstance(wus, str):
            self.create1(cursor, wus, priority)
        elif isinstance(wus, tuple) or isinstance(wus, list):
            for wu in wus:
                self.create1(cursor, wu, priority)
        else:
            cursor.close()
            raise Exception
        self.db.commit()
        cursor.close()

    def assign(self, clientid):
        """ Finds an available workunit and assigns it to clientid.
            Returns False of no available workunit exists """
        cursor = self.db.cursor()
        r = self.where(cursor, limit = 1, order=("priority", "DESC"), eq={"status": WuStatus.AVAILABLE})
        if len(r) == 1:
            self.data = r[0].data
            self._checkstatus(WuStatus.AVAILABLE)
            if debug > 0:
                self.check()
            d = {"status": WuStatus.ASSIGNED, 
                 "assignedclient": clientid,
                 "timeassigned": str(datetime.now())}
            self.update_wu(cursor, d)
            self.db.commit()
        cursor.close()
        return len(r) == 1

    def result(self, wuid, clientid, files, errorcode = None, failedcommand = None):
        cursor = self.db.cursor()
        if not self.get_by_wuid(cursor, wuid):
            cursor.close()
            raise Exception
        self._checkstatus(WuStatus.ASSIGNED)
        if debug > 0:
            self.check()
        d = {"resultclient": clientid,
             "errorcode": errorcode,
             "failedcommand": failedcommand, 
             "timeresult": str(datetime.now())}
        if errorcode == 0:
           d["status"] = WuStatus.RECEIVED_OK
        else:
            d["status"] = WuStatus.RECEIVED_ERROR
        self.update_wu(cursor, d)
        self.add_files(cursor, files)
        self.db.commit()
        cursor.close()

    def verification(self, wuid, ok):
        cursor = self.db.cursor()
        if not self.get_by_wuid(cursor, wuid):
            cursor.close()
            raise Exception
        self.data = r[0].data
        self._checkstatus(WuStatus.RECEIVED_OK)
        if debug > 0:
            self.check()
        d = {["timeverified"]: str(datetime.now())}
        if ok:
            d["status"] = WuStatus.VERIFIED_OK
        else:
            d["status"] = WuStatus.VERIFIED_ERROR
        self.update_wu(cursor, d)
        self.db.commit()
        cursor.close()

    def cancel(self, wuid):
        cursor = self.db.cursor()
        if not self.get_by_wuid(cursor, wuid):
            cursor.close()
            raise Exception
        self.data = r[0].data
        d = {"status": WuStatus.CANCELLED}
        self.update_wu(cursor, d)
        self.db.commit()
        cursor.close()

    def query(self, limit = None, **conditions):
        cursor = self.db.cursor()
        r = self.where(cursor, limit=limit, **conditions)
        cursor.close()
        return r
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


        


if __name__ == '__main__': # {
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-create', action="store_true", required=False, help='Create the database tables if they do not exist')
    parser.add_argument('-add', action="store_true", required=False, help='Add new work units. Contents of WU(s) are read from stdin, separated by blank line')
    parser.add_argument('-assign', required = False, nargs = 1, metavar = 'clientid', help = 'Assign an available WU to clientid')
    parser.add_argument('-prio', required = False, nargs = 1, metavar = 'N', help = 'If used with -add, newly added WUs receive priority N')

    for arg in ("avail", "assigned", "receivedok", "receivederr", "all"):
        parser.add_argument('-' + arg, action="store_true", required = False)
    for arg in ("dbname", "debug"):
        parser.add_argument('-' + arg, required = False, nargs = 1)
    # Parse command line, store as dictionary
    args = vars(parser.parse_args())
    # print(args)

    dbname = "wudb"
    if args["dbname"]:
        dbname = args["dbname"]

    if args["debug"]:
        debug = int(args["debug"][0])
    prio = 0
    if args["prio"]:
        debug = int(args["prio"][0])

    db = WuDb(dbname)
    wu = WuActiveRecord(db)
    
    if args["create"]:
        wu.create_tables()
    if args["add"]:
        s = ""
        wus = []
        for line in sys.stdin:
            if line == "\n":
                wus.append(s)
                s = ""
            else:
                s = s + line
        if s != "":
            wus.append(s)
        wu.create(wus, priority=prio)
    # Functions for queries
    if args["avail"]:
        wus = wu.query(eq={"status": WuStatus.AVAILABLE})
        print("Available workunits: ")
        for wu in wus:
            print (str(wu))
    if args["assigned"]:
        wus = wu.query(eq={"status": WuStatus.ASSIGNED})
        print("Assigned workunits: ")
        for wu in wus:
            print (str(wu))
    if args["receivedok"]:
        wus = wu.query(eq={"status": WuStatus.RECEIVED_OK})
        print("Received ok workunits: ")
        for wu in wus:
            print (str(wu))
    if args["receivederr"]:
        wus = wu.query(eq={"status": WuStatus.RECEIVED_ERROR})
        print("Received with error workunits: ")
        for wu in wus:
            print (str(wu))
    if args["all"]:
        wus = wu.query()
        print("Existing workunits: ")
        for wu in wus:
            print (str(wu))
    # Functions for testing
    if args["assign"]:
        clientid = args["assign"][0]
        wu.assign(clientid)
    
    db.close()
# }
