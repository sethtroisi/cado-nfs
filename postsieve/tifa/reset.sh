#!/bin/sh

#
# Reset the SCons environment to build the project from a clean state.
#
rm -fr .sconf_temp
rm -f  .sconsign.dblite
rm -f  config.log
