#!/bin/bash

# This installs cmake if it is not found in the current system

name=cmake
version=2.6.3
package=${name}-${version}.tar.gz
url=http://www.cmake.org/files/v2.6/${package}

me=`pwd`
tmpdir=`mktemp -d /tmp/${name}-build.XXXXXXXX`
cd $tmpdir
rm -f ${package}
wget="`which wget 2>/dev/null`"
if [ "$?" = "0" ] ; then
    wget $url
else
    curl="`which curl 2>/dev/null`"
    if [ "$?" = "0" ] ; then
        curl $url > $package
    else
        echo "Need either wget or curl to get $url" >&2
        exit 1
    fi
fi
tar xzf ${package}
# Better unset these, since we expect cmake to do the right thing without
# our custom flags !
unset CC
unset CXX
unset CFLAGS
unset CXXFLAGS
cd ${name}-${version}
./configure --prefix=$me/${name}
make -j 4
make install
rm -rf $tmpdir
