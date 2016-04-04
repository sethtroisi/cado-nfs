#!/usr/bin/env python3

import toplevel
import re
import os
import sys
import subprocess
import locale

pathdict = dict()

if not re.search("^/", "@CMAKE_INSTALL_PREFIX@"):
    # We are not in the installed tree, but in the source tree. We need
    # to find the .py libraries.
    pathdict["source"] = os.path.dirname(sys.argv[0])
    pathdict["pylib"] = pathdict["source"] + "/scripts/cadofactor"
    pathdict["data"] = pathdict["source"] + "/parameters"
    # find out where the binaries are. We need to use
    # ./scripts/build_environment.sh for that.
    helper = pathdict["source"] + "/scripts/build_environment.sh"
    pipe = subprocess.Popen([helper, "--show"], stdout=subprocess.PIPE)
    loc = locale.getdefaultlocale()[1]
    if not loc:
        loc="ascii"
    output = pipe.communicate()[0].decode(loc)
    cado_bin_path = [x.split("=",2)[1] for x in output.split("\n") if re.match("^build_tree",x)][0]
    cado_bin_path = re.sub("^\"(.*)\"$", "\\1", cado_bin_path)
    pathdict["lib"] = cado_bin_path
    pathdict["bin"] = cado_bin_path
else:
    pathdict["pylib"] = "@CMAKE_INSTALL_PREFIX@/@LIBSUFFIX@/scripts/cadofactor"
    pathdict["data"] ="@CMAKE_INSTALL_PREFIX@/@DATASUFFIX@"
    # binaries are installed in subdirectories of $LIBDIR.
    pathdict["lib"] ="@CMAKE_INSTALL_PREFIX@/@LIBSUFFIX@"
    pathdict["bin"] ="@CMAKE_INSTALL_PREFIX@/@BINSUFFIX@"


toplevel_params = toplevel.Cado_NFS_toplevel()
for key, value in pathdict.items():
    toplevel_params.setpath(key, value)
parameters = toplevel_params.get_cooked_parameters()
use_mysql = parameters.get_simple("mysql",False)
if use_mysql:
    import wudb_mysql
    wudb_mysql.dbname = parameters.get_simple("name")
    print("Loading Mysql")
    from wudb_mysql import *
else:
    print("Loading Sqlite")
    from wudb_sqlite import *
