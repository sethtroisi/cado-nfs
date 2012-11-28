import dbm
import pickle
import pickletools
from datetime import datetime
from Workunit import Workunit

debug = 0


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
class WU_Status:
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

class db_WU_entry:
    def __init__(self, WU_data):
        self.data = {}
        self.data["status"] = WU_Status.AVAILABLE
        self.data["WU_data"] = WU_data
        self.data["timecreated"] = str(datetime.now())
        self.data["timeassigned"] = ""
        self.data["assignedclient"] = ""
        self.data["timeresult"] = ""
        self.data["resultclient"] = ""
        self.data["errorcode"] = ""
        self.data["filenames"] = []
        self.data["timeverified"] = ""
    
    def get_WU(self):
        if debug > 0:
            print ("db_WU_entry.get_WU(): self.data = " + str(self.data))
            print ('db_WU_entry.get_WU(): self.data["WU_data"] = ' + str(self.data["WU_data"]))
        return self.data["WU_data"]
    
    def assign(self, clientid):
        if not self.data["status"] == WU_Status.AVAILABLE:
            raise StatusUpdateError("WU " + self.data["WU_data"]["WORKUNIT"] + 
                                    " has status " + self.data["status"])
        self.data["status"] = WU_Status.ASSIGNED
        self.data["assignedclient"] = clientid
        self.data["timeassigned"] = str(datetime.now())

    def result(self, clientid, errorcode, filenames):
        if not self.data["status"] == WU_Status.ASSIGNED:
            raise StatusUpdateError("WU " + self.data["WU_data"]["WORKUNIT"] + 
                                    " has status " + self.data["status"])
        if errorcode == 0:
            self.data["status"] = WU_Status.RECEIVED_OK
        else:
            self.data["status"] = WU_Status.RECEIVED_ERROR
        self.data["clientid_res"] = clientid
        self.data["errorcode"] = errorcode
        self.data["filenames"] = filenames
        self.data["timeresult"] = str(datetime.now())

    def verification(self, ok):
        if not self.data["status"] == WU_Status.RECEIVED_OK:
            raise StatusUpdateError("WU " + self.data["WU_data"]["WORKUNIT"] + 
                                    " has status " + self.data["status"])
        if ok:
            self.data["status"] = WU_Status.VERIFIED_OK
        else:
            self.data["status"] = WU_Status.VERIFIED_ERROR
        self.data["timeverified"] = str(datetime.now())
    
# If we try to add a WU under a key (WUid) that is already in the database, 
# we raise an exception. Python's built-in KeyError is not a good choice,
# as it is specified as indicating a lookup with a key that does not exist.
class KeyCollisionError(LookupError):
    pass
    
class wudb:
    def __init__(self, filename):
        self.db = dbm.open(filename, 'c', 0o644)

    # Pickle and write an entry
    def _write(self, key, entry):
        if debug > 0:
            print ("wudb._write(" + key + "): entry = " + str(entry))
            print ("wudb._write(" + key + "): entry.data = " + str(entry.data))
        pickledata = pickle.dumps(entry)
        if debug > 1:
            pickletools.dis(pickledata, annotate=1)
        self.db[key] = pickledata
    
    # Read an entry and unpickle it
    def _read(self, key):
        if key in self.db:
            data = self.db[key]
            if debug > 1:
                pickletools.dis(data, annotate=1)
            entry = pickle.loads(data)
            if debug > 0:
                print ("wudb._read(" + str(key) + "): entry.data = " + str(entry.data))
            return entry
        else:
            return None
    
    def add(self, WU):
        """ Add WU under its WUid. Set added-time, available flag """
        WUid = WU.get_id()
        if debug > 0:
            print ("Adding WU " + WUid)
        if WUid in self.db:
            raise KeyCollisionError("Work unit id " + WUid + " already in database")
        self._write(WUid, db_WU_entry(WU))

    def get(self, WUid):
        entry = self._read(WUid)
        if entry is None:
            return None
        return entry.get_WU()
    
    def where_eq(self, key, value, limit=0):
        """ Get a WU where db_WU_entry[key] == value """
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
        entry = self._read(key)
        entry.assign(clientid)
        self._write(key, entry)
        
    def result(self, WUid, clientid, errorcode, filenames):
        """ Mark a database WU entry as having received a result """
        entry = self._read(key)
        entry.result(clientid, errorcode, filenames)
        self._write(key, entry)
    
    def verify(self, WUid, ok):
        """ Mark a database WU entry as having been verified, ok == True means 
            verification succeeded, ok = False means it failed """
        entry = self._read(key)
        entry.verification(ok)
        self._write(key, entry)
    
if __name__ == '__main__':
    db = wudb("wudb")
    wu = Workunit("WORKUNIT abc2")
    if db.get(wu.get_id()) == None:
        db.add(wu)
    available = db.where_eq("status", WU_Status.AVAILABLE)
    print("Available workunits: ")
    for wu in available:
        print (str(wu))
