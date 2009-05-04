#!/bin/bash

# This installs cmake if it is not found in the current system

name=cmake
version=2.6.3
package=${name}-${version}.tar.gz
url=http://www.cmake.org/files/v2.6/${package}

me=`pwd`
tmpdir=/tmp/${name}-build-`id -un`
mkdir $tmpdir
cd $tmpdir
rm -f ${package}
wget $url
tar xzf ${package}
cd ${name}-${version}
./configure --prefix=$me/${name}
make -j 4
make install
