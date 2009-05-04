#!/bin/bash

# This script is responsible of handing over the build process, in a
# proper out of source build directory. It takes care of calling cmake
# first if needed, and then cd's into the proper sub-directory of the
# build tree, and runs make there. The intent is that this script is
# called from a Makefile within the source tree.

# The location of the build tree is, by default,
# $cado_source_tree/build/`hostname`, but customizing it is easy.

up_path=${0%%scripts/call_cmake.sh}
if [ "$up_path" = "$0" ] ; then
    echo "Error: call_cmake.sh must be called with a path ending in scripts/call_cmake.sh" >&2
    exit 1
elif [ "$up_path" = "" ] ; then
    up_path=./
fi
called_from=`pwd`
absolute_path_of_source=`cd $up_path ; pwd`
relative_path_of_cwd=${called_from##$absolute_path_of_source}

########################################################################
# Set some default variables which can be overridden from local.sh

# The default behaviour
# build_tree="/tmp/cado-build-`id -un`"
: ${build_tree:="${up_path}build/`hostname`"}

# By default, we also avoid /usr/local ; of course, it may be overridden.
: ${PREFIX:="$absolute_path_of_source/dist"}

# XXX XXX XXX LOOK LOOK LOOK: here you've got an entry point for customizing.

# The source directory may contain a hint script with several useful
# preferences. Its job is to put environment variables. Quite notably,
# $build_tree is amongst these.
if [ -f "${up_path}local.sh" ] ; then
    . "${up_path}local.sh"
fi

# If no CFLAGS have been set yet, set something sensible: get optimization by
# default, as well as asserts.  If you want to disable this, use either
# local.sh or the environment to set an environment variable CFLAGS to be
# something non-empty, like CFLAGS=-g or CFLAGS=-O0 ; Simply having CFLAGS=
# <nothing> won't do, because bash makes no distinction between null and unset
# here.
: ${CFLAGS:=-O2}
: ${CXXFLAGS:=-O2}


########################################################################
# Arrange so that relevant stuff is passed to cmake -- the other end of
# the magic is in CMakeLists.txt. The two lists must agree.

export PREFIX
export CFLAGS
export CXXFLAGS
export CC
export CXX
export CADO_VERSION

if [ "$1" = "tidy" ] ; then
    echo "Wiping out $build_tree"
    rm -rf $build_tree
    exit 0
fi

if [ "$1" = "show" ] ; then
    echo "build_tree=$build_tree"
    echo "up_path=$up_path"
    echo "called_from=$called_from"
    echo "absolute_path_of_source=$absolute_path_of_source"
    echo "relative_path_of_cwd=$relative_path_of_cwd"
    echo "CFLAGS=$CFLAGS"
    echo "CXXFLAGS=$CXXFLAGS"
    echo "CC=$CC"
    echo "CXX=$CXX"
    echo "PREFIX=$PREFIX"
    exit 0
fi

# Make sure we have cmake, by the way !
cmake_path="`which cmake 2>/dev/null`"
if [ "$?" != "0" ] || ! [ -x "$cmake_path" ] ; then 
    cmake_path="$absolute_path_of_source/cmake/bin/cmake"
    if ! [ -x "$cmake_path" ] ; then
        echo "Need to get cmake first -- this takes long !"
        cd $up_path
        if ! scripts/install-cmake.sh ; then
            echo "cmake install Failed, sorry" >&2
            exit 1
        fi
        cd $called_from
    fi
fi

if [ ! -f "$build_tree/Makefile" ] ; then
    mkdir -p $build_tree
    (cd $build_tree ; $cmake_path $absolute_path_of_source)
fi

# Now cd into the target directory, and build everything required.
# Note that it's useful to kill MAKELEVEL, or otherwise we'll get scores
# and scores of ``Entering directory'' messages (sure, there's the
# --no-print-directory option -- but it's not the right cure here).
# env | grep -i make
unset MAKELEVEL
(cd $build_tree$relative_path_of_cwd ; make "$@")
