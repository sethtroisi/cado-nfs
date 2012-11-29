#!/usr/bin/env python3

import dbm
import sys
import pickle
import pickletools
from datetime import datetime
from workunit import Workunit

debug = 1


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

class StatusUpdateError(Exception):
    pass

# One entry in the WU DB, including the dictionary with the WU contents 
# (FILEs, COMMANDs, etc.) and info about the progress on this WU (when and 
# to whom assigned, received, etc.)

class DbWuEntry:
    def __init__(self, wu_data):
        self.data = {}
        self.data["status"] = WuStatus.AVAILABLE
        self.data["wu_data"] = wu_data
        self.data["timecreated"] = str(datetime.now())
        self.data["timeassigned"] = ""
        self.data["assignedclient"] = ""
        self.data["timeresult"] = ""
        self.data["resultclient"] = ""
        self.data["errorcode"] = ""
        self.data["filenames"] = []
        self.data["timeverified"] = ""

    def __str__(self):
        return str(self.data)
    
    def get_WU(self):
        if debug > 1:
            print ("DbWuEntry.get_WU(): self.data = " + str(self.data))
            print ('DbWuEntry.get_WU(): self.data["wu_data"] = ' + str(self.data["wu_data"]))
        return self.data["wu_data"]
    
    def assign(self, clientid):
        if not self.data["status"] == WuStatus.AVAILABLE:
            raise StatusUpdateError("WU " + self.data["wu_data"].get_id() + 
                                    " has status " + self.data["status"])
        self.data["status"] = WuStatus.ASSIGNED
        self.data["assignedclient"] = clientid
        self.data["timeassigned"] = str(datetime.now())

    def result(self, clientid, errorcode, filenames):
        if not self.data["status"] == WuStatus.ASSIGNED:
            raise StatusUpdateError("WU " + self.data["wu_data"].get_id() + 
                                    " has status " + str(self.data["status"]))
        if errorcode == 0:
            self.data["status"] = WuStatus.RECEIVED_OK
        else:
            self.data["status"] = WuStatus.RECEIVED_ERROR
        self.data["clientid_res"] = clientid
        self.data["errorcode"] = errorcode
        self.data["filenames"] = filenames
        self.data["timeresult"] = str(datetime.now())

    def verification(self, ok):
        if not self.data["status"] == WuStatus.RECEIVED_OK:
            raise StatusUpdateError("WU " + self.data["wu_data"].get_id() + 
                                    " has status " + self.data["status"])
        if ok:
            self.data["status"] = WuStatus.VERIFIED_OK
        else:
            self.data["status"] = WuStatus.VERIFIED_ERROR
        self.data["timeverified"] = str(datetime.now())
    
# If we try to add a WU under a key (WUid) that is already in the database, 
# we raise an exception. Python's built-in KeyError is not a good choice,
# as it is specified as indicating a lookup with a key that does not exist.
class KeyCollisionError(LookupError):
    pass
    
class WuDb:
    def __init__(self, filename):
        self.db = dbm.open(filename, 'c', 0o644)

    # Pickle and write an entry
    def _write(self, key, entry):
        if debug > 1:
            print ("wudb._write(" + key + "): entry = " + str(entry))
            print ("wudb._write(" + key + "): entry.data = " + str(entry.data))
        pickledata = pickle.dumps(entry)
        if debug > 1:
            pickletools.dis(pickledata, annotate=1)
        self.db[key] = pickledata
        if hasattr(self.db, "sync"):
            self.db.sync()
    
    # Read an entry and unpickle it
    def _read(self, key):
        if key in self.db:
            data = self.db[key]
            if debug > 1:
                pickletools.dis(data, annotate=1)
            entry = pickle.loads(data)
            if debug > 1:
                print ("wudb._read(" + str(key) + "): entry.data = " + str(entry.data))
            return entry
        else:
            return None
    
    def add(self, WU):
        """ Add WU under its WUid. Set added-time, available flag """
        WUid = WU.get_id()
        if debug > 0:
            print ("WuDb.add(" + WUid + ")")
        if WUid in self.db:
            raise KeyCollisionError("Work unit id " + WUid + " already in database")
        self._write(WUid, DbWuEntry(WU))

    def get(self, WUid):
        entry = self._read(WUid)
        if entry is None:
            return None
        return entry
    
    def where_eq(self, key, value, limit=0):
        """ Get a WU where DbWuEntry[key] == value """
        result = []
        found = 0;
        for wuid in self.db.keys():
            entry = self._read(wuid)
            if entry.data[key] == value:
                result.append(entry.get_WU())
                found = found + 1
                if limit > 0 and found == limit:
                    break
        return result

    def assign(self, WUid, clientid):
        """ Assign a database WU entry to a client """
        if debug > 0:
            print ("WuDb.assign(" + WUid + ", " + clientid + ")")
        entry = self._read(WUid)
        entry.assign(clientid)
        self._write(WUid, entry)
        
    def result(self, WUid, clientid, errorcode, filenames):
        """ Mark a database WU entry as having received a result """
        if debug > 0:
            print ("WuDb.result(" + WUid + ", " + clientid + ", " + str(errorcode) + \
                   ", " + str(filenames) + ")")
        entry = self._read(WUid)
        entry.result(clientid, errorcode, filenames)
        self._write(WUid, entry)
    
    def verify(self, WUid, ok):
        """ Mark a database WU entry as having been verified, ok == True means 
            verification succeeded, ok = False means it failed """
        if debug > 0:
            print ("WuDb.verify(" + WUid + ", " + ok + ")")
        entry = self._read(WUid)
        entry.verification(ok)
        self._write(WUid, entry)
    
if __name__ == '__main__':
    dbname = "wudb"
    for arg in sys.argv:
        args = arg.split("=")
        if args[0] == "-dbname":
            dbname = args[1]
        if args[0] == "-add":
            db = WuDb(dbname)
            wutext = sys.stdin.read()
            wu = Workunit(wutext)
            db.add(wu)
        if args[0] == "-avail":
            db = WuDb(dbname)
            available = db.where_eq("status", WuStatus.AVAILABLE)
            print("Available workunits: ")
            for wu in available:
                print (str(wu))
        if args[0] == "-all":
            db = WuDb(dbname)
            for key in db.db.keys():
                print ("Workunit: " + key.decode())
                print(db.get(key))
