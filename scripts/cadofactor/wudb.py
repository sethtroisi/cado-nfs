#!/usr/bin/env python3

# TODO: 
# FILES table: OBSOLETE column
#     OBSOLETE says that this file was replaced by a newer version, for example checkrels may want to create a new file with only part of the output. 
#     Should this be a pointer to the file that replaced the obsolete one? How to signal a file that is obsolete, but not replaced by anything?
#     If one file is replaced by several (say, due to a data corruption in the middle), we need a 1:n relationship. If several files are replaced by one (merge), 
#     we need n:1. What do? Do we really want an n:n relationship here? Disallow fragmenting files, or maybe simply not track it in the DB if we do?
# FILES table: CHECKSUM column
#     We need a fast check that the information stored in the DB still accurately reflects the file system contents. The test should also warn about files in upload/ which are not listed in DB


# Make Python 2.7 use the print() syntax from Python 3
from __future__ import print_function

import sys
import sqlite3
import threading
import traceback
import collections
import abc
from datetime import datetime
from workunit import Workunit
if sys.version_info.major == 3:
    from queue import Queue
else:
    from Queue import Queue
import patterns
import signalhandler

debug = 1

def diag(level, text, var = None):
    if debug > level:
        if var is None:
            print (text, file=sys.stderr)
        else:
            print (text + str(var), file=sys.stderr)

def join3(l, pre = None, post = None, sep = ", "):
    """ 
    If any parameter is None, it is interpreted as the empty string 
    >>> join3 ( ('a'), pre = "+", post = "-", sep = ", ")
    '+a-'
    >>> join3 ( ('a', 'b'), pre = "+", post = "-", sep = ", ")
    '+a-, +b-'
    >>> join3 ( ('a', 'b'))
    'a, b'
    >>> join3 ( ('a', 'b', 'c'), pre = "+", post = "-", sep = ", ")
    '+a-, +b-, +c-'
    """
    if pre is None:
        pre = ""
    if post is None:
        post = ""
    if sep is None:
        sep = "";
    return sep.join([pre + k + post for k in l])

def dict_join3(d, sep=None, op=None, pre=None, post=None):
    """ 
    If any parameter is None, it is interpreted as the empty string 
    >>> dict_join3 ( {"a": "1", "b": "2"}, sep = ",", op = "=", pre="-", post="+")
    '-a=1+,-b=2+'
    """
    if pre is None:
        pre = ""
    if post is None:
        post = ""
    if sep is None:
        sep = "";
    if op is None:
        op = ""
    return sep.join([pre + op.join(k) + post for k in d.items()])

# Dummy class for defining "constants"
class WuStatus:
    AVAILABLE = 0
    ASSIGNED = 1
    RECEIVED_OK = 2
    RECEIVED_ERROR = 3
    VERIFIED_OK = 4
    VERIFIED_ERROR = 5
    CANCELLED = 6

    @classmethod
    def check(cls, status):
        """ Check whether status is equal to one of the constants """
        assert status in (cls.AVAILABLE, cls.ASSIGNED, cls.RECEIVED_OK, 
            cls.RECEIVED_ERROR, cls.VERIFIED_OK, cls.VERIFIED_ERROR, 
            cls.CANCELLED)

# If we try to update the status in any way other than progressive 
# (AVAILABLE -> ASSIGNED -> ...), we raise this exception
class StatusUpdateError(Exception):
    pass

class MyCursor(sqlite3.Cursor):
    """ This class represents a DB cursor and provides convenience functions 
        around SQL queries. In particular it is meant to provide an  
        (1) an interface to SQL functionality via method calls with parameters, 
        and 
        (2) hiding some particularities of the SQL variant of the underlying 
            DBMS as far as possible """
        
    # This is used in where queries; it converts from named arguments such as 
    # "eq" to a binary operator such as "="
    name_to_operator = {"lt": "<", "le": "<=", "eq": "=", "ge": ">=", "gt" : ">", "ne": "!="}
    
    def __init__(self, conn):
        # Enable foreign key support
        super(MyCursor, self).__init__(conn)
        self._exec("PRAGMA foreign_keys = ON;")

    @staticmethod
    def _without_None(d):
        """ Return a copy of the dictionary d, but without entries whose values 
            are None """
        return {k[0]:k[1] for k in d.items() if k[1] is not None}

    @staticmethod
    def as_string(d):
        if d is None:
            return ""
        else:
            return ", " + dict_join3(d, sep=", ", op=" AS ")
    
    @classmethod
    def _where_str(cls, name, **args):
        where = ""
        values = []
        for opname in args:
            if args[opname] is None:
                continue
            if where == "":
                where = " " + name + " "
            else:
                where = where + " AND "
            where = where + join3(args[opname].keys(), post = " " + cls.name_to_operator[opname] + " ?", sep = " AND ")
            values = values + list(args[opname].values())
        return (where, values)

    def _exec(self, command, values = None):
        """ Wrapper around self.execute() that prints arguments 
            for debugging and retries in case of "database locked" exception """
        if debug > 1:
            # FIXME: should be the caller's class name, as _exec could be 
            # called from outside this class
            classname = self.__class__.__name__
            parent = sys._getframe(1).f_code.co_name
            diag (1, classname + "." + parent + "(): command = " + command)
            if not values is None:
                diag (1, classname + "." + parent + "(): values = ", values)
        i = 0
        while True:
            try:
                if values is None:
                    self.execute(command)
                else:
                    self.execute(command, values)
                break
            except sqlite3.OperationalError as e:
                if i == 10 or str(e) != "database is locked":
                    raise

    def create_table(self, table, layout):
        """ Creates a table with fields as described in the layout parameter """
        command = "CREATE TABLE IF NOT EXISTS " + table + \
            "( " + ", ".join(" ".join(k) for k in layout) + " );"
        self._exec (command)
    
    def create_index(self, name, table, columns):
        """ Creates an index with fields as described in the columns list """
        command = "CREATE INDEX IF NOT EXISTS " + name + " ON " + \
            table + "( " + ", ".join(columns) + " );"
        self._exec (command)
    
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
        self._exec(command, values)
        rowid = self.lastrowid
        return rowid

    def update(self, table, d, **conditions):
        """ Update fields of an existing entry. conditions specifies the where 
            clause to use for to update, entries in the dictionary d are the 
            fields and their values to update """
        # UPDATE table SET column_1=value1, column2=value_2, ..., 
        # column_n=value_n WHERE column_n+1=value_n+1, ...,
        setstr = " SET " + join3(d.keys(), post = " = ?", sep = ", ")
        (wherestr, wherevalues) = self._where_str("WHERE", **conditions)
        command = "UPDATE " + table + setstr + wherestr
        values = list(d.values()) + wherevalues
        self._exec(command, values)
    
    def where(self, joinsource, col_alias = None, limit = None, order = None, 
              **conditions):
        """ Get a up to "limit" table rows (limit == 0: no limit) where 
            the key:value pairs of the dictionary "conditions" are set to the 
            same value in the database table """

        # Table/Column names cannot be substituted, so include in query directly.
        (WHERE, values) = self._where_str("WHERE", **conditions)

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

        AS = self.as_string(col_alias);

        command = "SELECT *" + AS + " FROM " + joinsource + WHERE + \
            ORDER + LIMIT + ";"
        self._exec(command, values)
        
    def count(self, joinsource, **conditions):
        """ Count rows where the key:value pairs of the dictionary "conditions" are 
            set to the same value in the database table """

        # Table/Column names cannot be substituted, so include in query directly.
        (WHERE, values) = self._where_str("WHERE", **conditions)

        command = "SELECT COUNT(*)" + " FROM " + joinsource + WHERE + ";"
        self._exec(command, values)
        r = self.fetchone()
        return int(r[0])
        
    def delete(self, table, **conditions):
        """ Delete the rows specified by conditions """
        (WHERE, values) = self._where_str("WHERE", **conditions)
        command = "DELETE FROM " + table + WHERE + ";"
        self._exec(command, values)

    def where_as_dict(self, joinsource, col_alias = None, limit = None, 
                      order = None, **conditions):
        self.where(joinsource, col_alias=col_alias, limit=limit, 
                      order=order, **conditions)
        # cursor.description is a list of lists, where the first element of 
        # each inner list is the column name
        result = []
        desc = [k[0] for k in self.description]
        row = self.fetchone()
        while row is not None:
            diag (2, "MyCursor.where_as_dict(): row = ", row)
            result.append(dict(zip(desc, row)))
            row = self.fetchone()
        return result


class DbTable(object):
    """ A class template defining access methods to a database table """
    
    @staticmethod
    def _subdict(d, l):
        """ Returns a dictionary of those key:value pairs of d for which key 
            exists l """
        if d is None:
            return None
        return {k:d[k] for k in d.keys() if k in l}

    def _get_colnames(self):
        return [k[0] for k in self.fields]

    def getname(self):
        return self.tablename

    def getpk(self):
        return self.primarykey

    def dictextract(self, d):
        """ Return a dictionary with all those key:value pairs of d
            for which key is in self._get_colnames() """
        return self._subdict(d, self._get_colnames())

    def create(self, cursor):
        fields = list(self.fields)
        if self.references:
            # If this table references another table, we use the primary
            # key of the referenced table as the foreign key name
            r = self.references # referenced table
            fk = (r.getpk(), "INTEGER", "REFERENCES " + 
                  r.getname() + " (" + r.getpk() + ") ")
            fields.append(fk)
        cursor.create_table(self.tablename, fields)
        if self.references:
            # We always create an index on the foreign key
            cursor.create_index(self.tablename + "_pkindex", r.tablename, 
                                (fk[0], ))
        for indexname in self.index:
            cursor.create_index(self.tablename + "_" + indexname, 
                                self.tablename, self.index[indexname])

    def insert(self, cursor, values, foreign=None):
        """ Insert a new row into this table. The column:value pairs are 
            specified key:value pairs of the dictionary d. 
            The database's row id for the new entry is stored in 
            d[primarykey] """
        d = self.dictextract(values)
        assert self.primarykey not in d or d[self.primarykey] is None
        # If a foreign key is specified in foreign, add it to the column
        # that is marked as being a foreign key
        if foreign:
            r = self.references.primarykey
            assert not r in d or d[r] is None
            d[r] = foreign
        values[self.primarykey] = cursor.insert(self.tablename, d)

    def insert_list(self, cursor, values, foreign=None):
        for v in values:
            self.insert(cursor, v, foreign)

    def update(self, cursor, d, **conditions):
        """ Update an existing row in this table. The column:value pairs to 
            be written are specified key:value pairs of the dictionary d """
        cursor.update(self.tablename, d, **conditions)

    def delete(self, cursor, **conditions):
        """ Delete an existing row in this table """
        cursor.delete(self.tablename, **conditions)

    def where(self, cursor, limit = None, order = None, **conditions):
        assert order is None or order[0] in self._get_colnames()
        return cursor.where_as_dict(self.tablename, limit=limit, 
                                    order=order, **conditions)


class WuTable(DbTable):
    tablename = "workunits"
    fields = (
        ("wurowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"), 
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
    primarykey = fields[0][0]
    references = None
    index = {"wuidindex": (fields[1][0],), "statusindex" : (fields[2][0],)}

class FilesTable(DbTable):
    tablename = "files"
    fields = (
        ("filesrowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"), 
        ("filename", "TEXT", ""), 
        ("path", "TEXT", "UNIQUE NOT NULL")
    )
    primarykey = fields[0][0]
    references = WuTable()
    index = {}


class DictDbTable(DbTable):
    fields = (
        ("rowid", "INTEGER PRIMARY KEY ASC", "UNIQUE NOT NULL"),
        ("key", "TEXT", "UNIQUE NOT NULL"),
        ("type", "INTEGER", "NOT NULL"),
        ("value", "TEXT", "")
        )
    primarykey = fields[0][0]
    references = None
    index = {"keyindex": ("key",)}
    def __init__(self, *args, name = None, **kwargs):
        self.tablename = name
        super().__init__(*args, **kwargs)


class DictDbAccess(collections.MutableMapping):
    """ A DB-backed flat dictionary.
    
    Flat means that the value of each dictionary entry must be a type that
    the underlying DB understands, like integers, strings, etc., but not
    collections or other complex types.
    
    A copy of all the data in the table is kept in memory; read accesses 
    are always served from the in-memory dict. Write accesses write through
    to the DB.
    
    >>> conn = sqlite3.connect(':memory:')
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {}
    True
    >>> d['a'] = '1'
    >>> d == {'a': '1'}
    True
    >>> d['a'] = 2
    >>> d == {'a': 2}
    True
    >>> d['b'] = '3'
    >>> d == {'a': 2, 'b': '3'}
    True
    >>> del(d)
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {'a': 2, 'b': '3'}
    True
    >>> del(d['b'])
    >>> d == {'a': 2}
    True
    >>> d.setdefault('a', '3')
    2
    >>> d == {'a': 2}
    True
    >>> d.setdefault('b', 3.0)
    3.0
    >>> d == {'a': 2, 'b': 3.0}
    True
    >>> d.setdefault(None, {'a': '3', 'c': '4'})
    >>> d == {'a': 2, 'b': 3.0, 'c': '4'}
    True
    >>> d.update({'a': '3', 'd': True})
    >>> d == {'a': '3', 'b': 3.0, 'c': '4', 'd': True}
    True
    >>> del(d)
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {'a': '3', 'b': 3.0, 'c': '4', 'd': True}
    True
    >>> d.clear('a', 'd')
    >>> d == {'b': 3.0, 'c': '4'}
    True
    >>> del(d)
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {'b': 3.0, 'c': '4'}
    True
    >>> d.clear()
    >>> d == {}
    True
    >>> del(d)
    >>> d = DictDbAccess(conn, 'test')
    >>> d == {}
    True
    """
    
    types = (str, int, float, bool)
    
    def __init__(self, db, name):
        ''' Attaches to a DB table and reads values stored therein. 
        
        db can be a string giving the file name for the DB (same as for 
        sqlite3.connect()), or an open DB connection. The latter is allowed 
        primarily for making the doctest work, so we can reuse the same 
        memory-backed DB connection, but it may be useful in other contexts.
        Note that writes to the dict cause a commit() on the DB connection, 
        so sharing a DB connection may be ill-advised if you want control 
        over when commits happen.
        '''
        
        if isinstance(db, str):
            self._conn = sqlite3.connect(db)
            self._ownconn = True
        else:
            self._conn = db
            self._ownconn = False
        self._table = DictDbTable(name = name)
        # Create an empty table if none exists
        cursor = self._conn.cursor(MyCursor)
        self._table.create(cursor);
        # Get the entries currently stored in the DB
        self._data = self._getall()
        cursor.close()
    
    # Implement the abstract methods defined by collections.MutableMapping
    # All but __del__ and __setitem__ are simply passed through to the self._data
    # dictionary
    
    def __getitem__(self, key):
        return self._data.__getitem__(key)
    
    def __iter__(self):
        return self._data.__iter__()
    
    def  __len__(self):
        return self._data.__len__()
    
    def __str__(self):
        return self._data.__str__()
    
    def __del__(self):
        """ Close the DB connection and delete the dictionary """
        if self._ownconn:
            self._conn.close()
        # http://docs.python.org/2/reference/datamodel.html#object.__del__
        # MutableMapping does not have __del__, but in a complex class 
        # hierarchy, it may not be next in the MRO
        if hasattr(super(), "__del__"):
            super().__del__()
    
    def __convert_value(self, row):
        valuestr = row["value"]
        valuetype = row["type"]
        # Look up constructor for this type
        typecon = self.types[int(valuetype)]
        # Bool is handled separately as bool("False") == True
        if typecon == bool:
            if valuestr == "True":
                return True
            elif valuestr == "False":
                return False
            else:
                raise ValueError("Value %s invalid for Bool type", valuestr)
        return typecon(valuestr)
    
    def __get_type_idx(self, value):
        valuetype = type(value)
        for (idx, t) in enumerate(self.types):
            if valuetype == t:
                return idx
        raise TypeError("Type %s not supported" % str(valuetype))
    
    def _getall(self):
        """ Reads the whole table and returns it as a dict """
        cursor = self._conn.cursor(MyCursor)
        rows = self._table.where(cursor)
        cursor.close()
        dict = {r["key"]: self.__convert_value(r) for r in rows}
        return dict
    
    def __setitem_nocommit(self, cursor, key, value):
        """ Set dictioary key to value and update/insert into table,
        but don't commit
        """
        update = {"value": str(value), "type": self.__get_type_idx(value)}
        if key in self._data:
            # Update the table row where column "key" equals key
            self._table.update(cursor, update, eq={"key": key})
        else:
            # Insert a new row
            update["key"] = key
            self._table.insert(cursor, update)
        # Update the in-memory dict
        self._data[key] = value
    
    def __setitem__(self, key, value):
        """ Set a dict entry to a value and update the DB """
        cursor = self._conn.cursor(MyCursor)
        self.__setitem_nocommit(cursor, key, value)
        self._conn.commit()
        cursor.close()
    
    def __delitem__(self, key):
        """ Delete a key from the dictionary """
        del(self._data[key])
        cursor = self._conn.cursor(MyCursor)
        self._table.delete(cursor, eq={"key": key})
        self._conn.commit()
        cursor.close()
    
    def setdefault(self, key, default = None):
        ''' Setdefault function that allows a mapping as input
        
        Values from default dict are merged into self, *not* overwriting
        existing values in self '''
        if key is None and isinstance(default, collections.Mapping):
            cursor = self._conn.cursor(MyCursor)
            for (key, value) in default.items():
                if not key in self:
                    self.__setitem_nocommit(cursor, key, value)
            self._conn.commit()
            cursor.close()
            return None
        elif not key in self:
            self.__setitem__(key, default)
        return self._data[key]
    
    def update(self, other):
        cursor = self._conn.cursor(MyCursor)
        for (key, value) in other.items():
            self.__setitem_nocommit(cursor, key, value)
        self._conn.commit()
        cursor.close()
    
    def clear(self, *args):
        """ Overriden clear that allows removing several keys atomically """
        cursor = self._conn.cursor(MyCursor)
        if not args:
            self._data.clear()
            self._table.delete(cursor)
        else:
            for key in args:
                del(self._data[key])
                self._table.delete(cursor, eq={"key": key})
        self._conn.commit()
        cursor.close()


class Mapper(object):
    """ This class translates between application objects, i.e., Python 
        directories, and the relational data layout in an SQL DB, i.e.,
        one or more tables which possibly have foreign key relationships 
        that map to hierarchical data structures. For now, only one 
        foreign key / subdirectory."""

    def __init__(self, table, subtables = None):
        self.table = table
        self.subtables = {}
        if subtables:
            for s in subtables.keys():
                self.subtables[s] = Mapper(subtables[s])

    def __sub_dict(self, d):
        """ For each key "k" that has a subtable assigned in "self.subtables",
        pop the entry with key "k" from "d", and store it in a new directory
        which is returned. I.e., the directory d is separated into 
        two parts: the part which corresponds to subtables and is the return 
        value, and the rest which is left in the input dictionary. """
        sub_dict = {}
        for s in self.subtables.keys():
            # Don't store s:None entries even if they exist in d
            t = d.pop(s, None)
            if not t is None:
                sub_dict[s] = t
        return sub_dict

    def getname(self):
        return self.table.getname()

    def getpk(self):
        return self.table.getpk()

    def create(self, cursor):
        self.table.create(cursor)
        for t in self.subtables.values():
            t.create(cursor)

    def insert(self, cursor, wus, foreign=None):
        pk = self.getpk()
        for wu in wus:
            # Make copy so sub_dict does not change caller's data
            wuc = wu.copy()
            # Split off entries that refer to subtables
            sub_dict = self.__sub_dict(wuc)
            # We add the entries in wuc only if it does not have a primary 
            # key yet. If it does have a primary key, we add only the data 
            # for the subtables
            if not pk in wuc:
                self.table.insert(cursor, wuc, foreign=foreign)
                # Copy primary key into caller's data
                wu[pk] = wuc[pk]
            for subtable_name in sub_dict.keys():
                self.subtables[subtable_name].insert(
                    cursor, sub_dict[subtable_name], foreign=wu[pk])

    def update(self, cursor, wus):
        pk = self.getpk()
        for wu in wus:
            assert not wu[pk] is None
            wuc = wu.copy()
            sub_dict = self.__sub_dict(wuc)
            rowid = wuc.pop(pk, None)
            if rowid:
                self.table.update(cursor, wuc, {wp: rowid})
            for s in sub.keys:
                self.subtables[s].update(cursor, sub_dict[s])
    
    def count(self, cursor, **cond):
        joinsource = self.table.tablename
        return cursor.count(joinsource, **cond)
    
    def where(self, cursor, limit = None, order = None, **cond):
        pk = self.getpk()
        joinsource = self.table.tablename
        for s in self.subtables.keys():
            # FIXME: this probably breaks with more than 2 tables
            joinsource = joinsource + " LEFT JOIN " + \
                self.subtables[s].getname() + \
                " USING ( " + pk + " )"
        # FIXME: don't get result rows as dict! Leave as tuple and
        # take them apart positionally
        rows = cursor.where_as_dict(joinsource, limit=limit, order=order, 
                                    **cond)
        wus = []
        for r in rows:
            # Collapse rows with identical primary key
            if len(wus) == 0 or r[pk] != wus[-1][pk]:
                wus.append(self.table.dictextract(r))
                for s in self.subtables.keys():
                    wus[-1][s] = None

            for (sn, sm) in self.subtables.items():
                spk = sm.getpk()
                if spk in r and not r[spk] is None:
                    if wus[-1][sn] == None:
                        # If this sub-array is empty, init it
                        wus[-1][sn] = [sm.table.dictextract(r)]
                    elif r[spk] != wus[-1][sn][-1][spk]:
                        # If not empty, and primary key of sub-table is not
                        # same as in previous entry, add it
                        wus[-1][sn].append(sm.table.dictextract(r))
        return wus

class WuAccess(object): # {
    """ This class maps between the WORKUNIT and FILES tables 
        and a dictionary 
        {"wuid": string, ..., "timeverified": string, "files": list}
        where list is None or a list of dictionaries of the from
        {"id": int, "wuid": string, "filename": string, "path": string
        Operations on instances of WuAcccess are directly carried 
        out on the database persistent storage, i.e., they behave kind 
        of as if the WuAccess instance were itself a persistent 
        storage device """
    
    def __init__(self, dbfilename):
        self.conn = sqlite3.connect(dbfilename)
        self.mapper = Mapper(WuTable(), {"files": FilesTable()})
    
    def __del__(self):
        self.conn.close()
    
    @staticmethod
    def to_str(wus):
        r = []
        for wu in wus:
            s = "Workunit " + str(wu["wuid"]) + ":\n"
            for (k,v) in wu.items():
                if k != "wuid" and k != "files":
                    s = s + "  " + k + ": " + repr(v) + "\n"
            if "files" in wu:
                s = s + "  Associated files:\n"
                if wu["files"] is None:
                    s = s + "    None\n"
                else:
                    for f in wu["files"]:
                        s = s + "    " + str(f) + "\n"
            r.append(s)
        return '\n'.join(r)

    @staticmethod
    def _checkstatus(wu, status):
        diag (2, "WuAccess._checkstatus(" + str(wu) + ", " + str(status) + ")")
        if wu["status"] != status:
            msg = "Workunit " + str(wu["wuid"]) + " has status " + str(wu["status"]) \
                + ", expected " + str(status)
            diag (0, "WuAccess._checkstatus(): " + msg)
            raise StatusUpdateError(msg)

    def check(self, data):
        status = data["status"]
        WuStatus.check(status)
        wu = Workunit(data["wu"])
        assert wu.get_id() == data["wuid"]
        if status > WuStatus.RECEIVED_ERROR:
            return
        if status == WuStatus.RECEIVED_ERROR:
            assert data["errorcode"] != 0
            return
        if status == WuStatus.RECEIVED_OK:
            assert data["errorcode"] is None or data["errorcode"] == 0
            return
        assert data["errorcode"] is None
        assert data["timeresult"] is None
        assert data["resultclient"] is None
        if status == WuStatus.ASSIGNED:
            return
        assert data["timeassigned"] is None
        assert data["assignedclient"] is None
        if status == WuStatus.AVAILABLE:
            return
        assert data["timecreated"] is None
        # etc.
    
    # Here come the application-visible functions that implement the 
    # "business logic": creating a new workunit from the text of a WU file,
    # assigning it to a client, receiving a result for the WU, marking it as
    # verified, or marking it as cancelled

    def add_files(self, cursor, files, wuid = None, rowid = None):
        # Exactly one must be given
        assert not wuid is None or not rowid is None
        assert wuid is None or rowid is None
        # FIXME: allow selecting row to update directly via wuid, without 
        # doing query for rowid first
        pk = self.mapper.getpk()
        if rowid is None:
            wu = get_by_wuid(cursor, wuid)
            if wu:
                rowid = wu[pk]
            else:
                return False
        d = ({"filename": f[0], "path": f[1]} for f in files)
        # These two should behave identically
        if True:
            self.mapper.insert(cursor, [{pk:rowid, "files": d},])
        else:
            self.mapper.subtables["files"].insert(cursor, d, foreign=rowid)

    def create_tables(self):
        cursor = self.conn.cursor(MyCursor)
        cursor._exec("PRAGMA journal_mode=WAL;")
        self.mapper.create(cursor)
        self.conn.commit()
        cursor.close()

    def create1(self, cursor, wutext, priority = None):
        d = {
            "wuid": Workunit(wutext).get_id(),
            "wu": wutext,
            "status": WuStatus.AVAILABLE,
            "timecreated": str(datetime.now())
            }
        if not priority is None:
            d["priority"] = priority
        # Insert directly into wu table
        self.mapper.table.insert(cursor, d)

    def create(self, wus, priority = None):
        """ Create new workunits from wus which contains the texts of the 
            workunit files """
        cursor = self.conn.cursor(MyCursor)
        if isinstance(wus, str):
            self.create1(cursor, wus, priority)
        else:
            for wu in wus:
                self.create1(cursor, wu, priority)
        self.conn.commit()
        cursor.close()

    def assign(self, clientid):
        """ Finds an available workunit and assigns it to clientid.
            Returns the text of the workunit, or None if no available 
            workunit exists """
        cursor = self.conn.cursor(MyCursor)
        r = self.mapper.table.where(cursor, limit = 1, 
                                    order=("priority", "DESC"), 
                                    eq={"status": WuStatus.AVAILABLE})
        assert len(r) <= 1
        if len(r) == 1:
            self._checkstatus(r[0], WuStatus.AVAILABLE)
            if debug > 0:
                self.check(r[0])
            d = {"status": WuStatus.ASSIGNED, 
                 "assignedclient": clientid,
                 "timeassigned": str(datetime.now())}
            pk = self.mapper.getpk()
            self.mapper.table.update(cursor, d, eq={pk:r[0][pk]})
            self.conn.commit()
            cursor.close()
            return r[0]["wu"]
        else:
            cursor.close()
            return None

    def get_by_wuid(self, cursor, wuid):
        r = self.mapper.where(cursor, eq={"wuid": wuid})
        assert len(r) <= 1
        if len(r) == 1:
            return r[0]
        else:
            return None

    def result(self, wuid, clientid, files, errorcode = None, 
               failedcommand = None):
        cursor = self.conn.cursor(MyCursor)
        data = self.get_by_wuid(cursor, wuid)
        if data is None:
            cursor.close()
            return False
        self._checkstatus(data, WuStatus.ASSIGNED)
        if debug > 0:
            self.check(data)
        d = {"resultclient": clientid,
             "errorcode": errorcode,
             "failedcommand": failedcommand, 
             "timeresult": str(datetime.now())}
        if errorcode is None or errorcode == 0:
           d["status"] = WuStatus.RECEIVED_OK
        else:
            d["status"] = WuStatus.RECEIVED_ERROR
        pk = self.mapper.getpk()
        self.mapper.table.update(cursor, d, eq={pk:data[pk]})
        self.add_files(cursor, files, rowid = data[pk])
        self.conn.commit()
        cursor.close()

    def verification(self, wuid, ok):
        cursor = self.conn.cursor(MyCursor)
        data = self.get_by_wuid(cursor, wuid)
        if data is None:
            cursor.close()
            return False
        # FIXME: should we do the update by wuid and skip these checks?
        self._checkstatus(data, WuStatus.RECEIVED_OK)
        if debug > 0:
            self.check(data)
        d = {"timeverified": str(datetime.now())}
        if ok:
            d["status"] = WuStatus.VERIFIED_OK
        else:
            d["status"] = WuStatus.VERIFIED_ERROR
        pk = self.mapper.getpk()
        self.mapper.table.update(cursor, d, eq={pk:data[pk]})
        self.conn.commit()
        cursor.close()

    def cancel(self, wuid):
        self.cancel_by_condition(eq={"wuid":wuid})
    
    def cancel_all_available(self):
        self.cancel_by_condition(eq={"status": 0})
    
    def cancel_all_assigned(self):
        self.cancel_by_condition(eq={"status": 1})
    
    def cancel_by_condition(self, **conditions):
        cursor = self.conn.cursor(MyCursor)
        d = {"status": WuStatus.CANCELLED}
        self.mapper.table.update(cursor, d, **conditions)
        self.conn.commit()
        cursor.close()

    def query(self, limit = None, **conditions):
        cursor = self.conn.cursor(MyCursor)
        r = self.mapper.where(cursor, limit=limit, **conditions)
        cursor.close()
        return r

    def count(self, **cond):
        cursor = self.conn.cursor(MyCursor)
        count = self.mapper.count(cursor, **cond)
        cursor.close()
        return count
    
    def count_available(self):
        return self.count(eq={"status": WuStatus.AVAILABLE})
    
    def get_one_result(self):
        r = self.query(limit = 1, eq={"status": WuStatus.RECEIVED_OK})
        if not r:
            r = self.query(limit = 1, eq={"status": WuStatus.RECEIVED_ERROR})
        if not r:
            return None
        else:
            return r[0]
#}

class WuResultMessage(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def get_wu_id(self):
        pass
    @abc.abstractmethod
    def get_output_files(self):
        pass
    @abc.abstractmethod
    def get_stdout(self, command_nr):
        pass
    @abc.abstractmethod
    def get_stderr(self, command_nr):
        pass
    @abc.abstractmethod
    def get_exitcode(self, command_nr):
        pass

class ResultInfo(WuResultMessage):
    def __init__(self, record):
        # record looks like this:
        # {'status': 0, 'errorcode': None, 'timeresult': None, 'wuid': 'testrun_polyselect_0-5000', 
        #  'wurowid': 1, 'timecreated': '2013-05-23 22:28:08.333310', 'timeverified': None, 
        #  'failedcommand': None, 'priority': None, 'wu': "WORKUNIT [..rest of workunit text...] \n", 
        #  'assignedclient': None, 'retryof': None, 'timeassigned': None, 'resultclient': None, 
        #  'files': None}
        self.record = record
    
    def get_wu_id(self):
        return self.record["wuid"]
    
    def get_output_files(self):
        """ Returns the list of output files of this workunit
        
        Only files that were specified in RESULT lines appear here;
        automatically captured stdout and stderr does not. 
        """
        files = []
        for f in self.record["files"]:
            if not f["filename"].startswith("stdout") and not f["filename"].startswith("stderr"):
                files.append(f["path"])
        return files
    
    def get_file(self, filename):
        """ Get the file location of the file that was uploaded under the
        name filename. Used internally.
        """
        for f in self.record["files"]:
            if f["filename"] == "filename":
                return f["path"]
        return None
    
    def get_stdout(self, command_nr):
        """ Return the path to the file that captured stdout of the 
        command_nr-th COMMAND in the workunit, or None if there was no stdout 
        output. Note that explicitly redirected stdout that was uploaded via 
        RESULT does not appear here, but in get_files()
        """
        return self.get_file("stdout%d" % command_nr)
    
    def get_stderr(self, command_nr):
        """ Like get_stdout(), but for stderr """
        return self.get_file("stderr%d" % command_nr)
    
    def get_exitcode(self, command_nr):
        """ Return the exit code of the command_nr-th command """
        if not self.record["failedcommand"] is None \
                and command_nr == int(self.record["failedcommand"]):
            return int(self.record["errorcode"])
        else:
            return 0


class DbListener(patterns.Observable):
    """ Class that queries the Workunit database for available results
    and sends them to its Observers.
    
    The query is triggered by receiving a SIGUSR1 (the instance subscribes to
    the signal handler relay), or by calling send_result().
    """
    def __init__(self, *args, db, **kwargs):
        super().__init__(*args, **kwargs)
        self.wuar = WuAccess(db)
    
    def send_result(self):
        # Check for results
        r = self.wuar.get_one_result()
        if not r:
            return False
        message = ResultInfo(r)
        self.notifyObservers(message)
        return True


class IdMap(object):
    """ Identity map. Ensures that DB-backed dictionaries on the same DB and
    of the same table name are instanciated only once. """
    def __init__(self):
        self.db_dicts = {}
    
    def make_db_dict(self, db, name):
        assert isinstance(db, str)
        key = (db, name)
        if not key in self.db_dicts:
            self.db_dicts[key] = DictDbAccess(db, name)
        return self.db_dicts[key]

# Singleton instance of IdMap
idmap = IdMap()

class DbAccess(object):
    """ Base class that lets subclasses create DB-backed dictionaries or 
    WuAccess instances on a database whose file name is specified in the db 
    parameter to __init__.
    Meant to be used as a cooperative class; it strips the db parameter from
    the parameter list and remembers it in a private variable so that it can 
    later be used to open DB connections.
    """
    
    def __init__(self, *args, db, **kwargs):
        super().__init__(*args, **kwargs)
        self.__db = db
    
    def get_db_filename(self):
        return self.__db
    
    def make_db_dict(self, name):
        return idmap.make_db_dict(self.__db, name)
    
    def make_wu_access(self):
        return WuAccess(self.__db)
    
    def make_db_listener(self):
        return DbListener(db = self.__db)


class DbWorker(DbAccess, threading.Thread):
    """Thread executing WuAccess requests from a given tasks queue"""
    def __init__(self, taskqueue, *args, daemon=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.taskqueue = taskqueue
        if not daemon is None:
            self.daemon = daemon
        self.start()
    
    def run(self):
        # One DB connection per thread. Created inside the new thread to make
        # sqlite happy
        wuar = self.make_wu_access()
        while True:
            # We expect a 4-tuple in the task queue. The elements of the tuple:
            # a 2-array, where element [0] receives the result of the DB call, 
            #  and [1] is an Event variable to notify the caller when the 
            #  result is available
            # fn_name, the name (as a string) of the WuAccess method to call
            # args, a tuple of positional arguments
            # kargs, a dictionary of keyword arguments
            (result_tuple, fn_name, args, kargs) = self.taskqueue.get()
            if fn_name == "terminate":
                break
            ev = result_tuple[1]
            # Assign to tuple in-place, so result is visible to caller. 
            # No slice etc. here which would create a copy of the array
            try: result_tuple[0] = getattr(wuar, fn_name)(*args, **kargs)
            except Exception as e: 
                traceback.print_exc()
            ev.set()
            self.taskqueue.task_done()

class DbRequest(object):
    """ Class that represents a request to a given WuAccess function.
        Used mostly so that DbThreadPool's __getattr__ can return a callable 
        that knows which of WuAccess's methods should be called by the 
        worker thread """
    def __init__(self, taskqueue, func):
        self.taskqueue = taskqueue
        self.func = func
    
    def do_task(self, *args, **kargs):
        """Add a task to the queue, wait for its completion, and return the result"""
        ev = threading.Event()
        result = [None, ev]
        self.taskqueue.put((result, self.func, args, kargs))
        ev.wait()
        return result[0]

class DbThreadPool(object):
    """Pool of threads consuming tasks from a queue"""
    def __init__(self, dbfilename, num_threads = 1):
        self.taskqueue = Queue(num_threads)
        self.pool = []
        for _ in range(num_threads):
            worker = DbWorker(self.taskqueue, daemon=True, db=dbfilename)
            self.pool.append(worker)

    def terminate(self):
        for t in self.pool:
            self.taskqueue.put((None, "terminate", None, None))
        self.wait_completion

    def wait_completion(self):
        """Wait for completion of all the tasks in the queue"""
        self.taskqueue.join()
    
    def __getattr__(self, name):
        """ Delegate calls to methods of WuAccess to a worker thread.
            If the called method exists in WuAccess, creates a new 
            DbRequest instance that remembers the name of the method that we 
            tried to call, and returns the DbRequest instance's do_task 
            method which will process the method call via the thread pool. 
            We need to go through a new object's method since we cannot make 
            the caller pass the name of the method to call to the thread pool 
            otherwise """
        if hasattr(WuAccess, name):
            task = DbRequest(self.taskqueue, name)
            return task.do_task
        else:
            raise AttributeError(name)


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

    queries = {"avail" : ("Available workunits: ", {"eq": {"status": WuStatus.AVAILABLE}}), 
               "assigned": ("Assigned workunits: ", {"eq": {"status": WuStatus.ASSIGNED}}), 
               "receivedok": ("Received ok workunits: ", {"eq":{"status": WuStatus.RECEIVED_OK}}), 
               "receivederr": ("Received with error workunits: ", {"eq": {"status": WuStatus.RECEIVED_ERROR}}), 
               "verifiedok": ("Verified ok workunits: ", {"eq": {"status": WuStatus.VERIFIED_OK}}), 
               "verifiederr": ("Verified with error workunits: ", {"eq": {"status": WuStatus.VERIFIED_ERROR}}), 
               "cancelled": ("Cancelled workunits: ", {"eq": {"status": WuStatus.CANCELLED}}), 
               "all": ("All existing workunits: ", {})
              }

    use_pool = False

    parser = argparse.ArgumentParser()
    parser.add_argument('-create', action="store_true", 
                        help='Create the database tables if they do not exist')
    parser.add_argument('-add', action="store_true", 
                        help='Add new workunits. Contents of WU(s) are ' 
                        'read from stdin, separated by blank line')
    parser.add_argument('-assign', nargs = 1, metavar = 'clientid', 
                        help = 'Assign an available WU to clientid')
    parser.add_argument('-cancel', nargs = 1, metavar = 'wuid', 
                        help = 'Cancel a WU with the given id')
    parser.add_argument('-prio', metavar = 'N', 
                        help = 'If used with -add, newly added WUs ' 
                        'receive priority N')
    parser.add_argument('-result', nargs = 4, 
                        metavar = ('wuid', 'clientid', 'filename', 'filepath'), 
                        help = 'Return a result for wu from client')
    parser.add_argument('-test', action="store_true", 
                        help='Run some self tests')

    for arg in queries:
        parser.add_argument('-' + arg, action="store_true", required = False)
    parser.add_argument('-dump', nargs='?', default = None, const = "all", 
                        metavar = "FIELD", 
                        help='Dump WU contents, optionally a single field')
    for arg in ("dbname", "debug"):
        parser.add_argument('-' + arg)
    # Parse command line, store as dictionary
    args = vars(parser.parse_args())
    # print(args)

    dbname = "wudb"
    if args["dbname"]:
        dbname = args["dbname"]

    if args["test"]:
        import doctest
        doctest.testmod()

    if args["debug"]:
        debug = int(args["debug"])
    prio = 0
    if args["prio"]:
        prio = int(args["prio"][0])
    
    if use_pool:
        db_pool = DbThreadPool(dbname)
    else:
        db_pool = WuAccess(dbname)
    
    if args["create"]:
        db_pool.create_tables()
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
        db_pool.create(wus, priority=prio)

    # Functions for queries
    for (arg, (msg, condition)) in queries.items():
        if not args[arg]:
            continue
        print(msg)
        if not args["dump"]:
            count = db_pool.count(**condition)
            print (count)
        elif args["dump"]:
            wus = db_pool.query(**condition)
            if wus is None:
                print("None")
                continue
            print (len(wus))
            if args["dump"] == "all":
                print(WuAccess.to_str(wus))
            else:
                field = args["dump"]
                for wu in wus:
                    print(wu[field])


    # Functions for testing
    if args["assign"]:
        clientid = args["assign"][0]
        wus = db_pool.assign(clientid)

    if args["cancel"]:
        wuid = args["cancel"][0]
        if wuid == "available":
            wus = db_pool.cancel_all_available()
        if wuid == "assigned":
            wus = db_pool.cancel_all_assigned()
        else:
            wus = db_pool.cancel(wuid)

    if args["result"]:
        (wuid, clientid, filename, filepath) = args["result"]
        db_pool.result(wuid, clientid, [(filename, filepath)])
    
    if use_pool:
        db_pool.terminate()
# }

# Local Variables:
# version-control: t
# End:
