#!/usr/bin/env python3

import toplevel
import re
import os
import sys
import subprocess
import locale

#if use_mysql:
#    import wudb_mysql
#    wudb_mysql.dbname = parameters.get_simple("name")
#    print("Loading Mysql")
#    import wudb_mysql
#    wudb_mysql.username = parameters.get_simple("mysql.username",'cado')
#    wudb_mysql.password = parameters.get_simple("mysql.password",'***REMOVED***')
#    from wudb_mysql import *
#else:
#    print("Loading Sqlite")
#    from wudb_sqlite import *

from wudb_sqlite import *
