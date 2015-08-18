#!/usr/bin/env bash

########################################################################
# This script is responsible of handing over the build process, in a
# proper out of source build directory. It takes care of calling cmake
# first if needed, and then cd's into the proper sub-directory of the
# build tree, and runs make there. The intent is that this script is
# called from a Makefile within the source tree.
# In particuar, the following tasks are done,
#  - check if the calling path is correct?
#  - if exists, parse the file ${up_path}local.sh
#  - check if cmake is installed, if not install it.
#  - "cmake" to generate Makefile
#  - "make"
#
# Type "make ?" for more options.
########################################################################

: ${MAKE=make}
export MAKE

if echo "${MAKEFLAGS}" | grep -q "jobserver-fds=0," ; then
    echo "# You are calling the top-level cado makefile with file descriptor 0 closed.">&2
    echo "# This is unsupported (at least for a parallel build), because in that">&2
    echo "# case GNU Make opens and uses a pipe on file descriptor 0, and we">&2
    echo "# suspect that cmake closes it right away, causing the compilation to">&2
    echo "# fail.">&2
    echo "#">&2
    echo "# Simple fix: make -j \$number_of_cpus < /dev/null">&2
    echo "#">&2
    exit 42
fi

########################################################################
# The location of the build tree is, by default,
# $cado_source_tree/build/`hostname`, but customizing it is easy.

up_path="${0%%scripts/call_cmake.sh}"
if [ "$up_path" = "$0" ] ; then
    echo "Error: call_cmake.sh must be called with a path ending in scripts/call_cmake.sh" >&2
    exit 1
elif [ "$up_path" = "" ] ; then
    up_path=./
fi
pwdP="pwd -P"
if ! $pwdP >/dev/null 2>&1 ; then
    pwdP=pwd
fi
called_from="`pwd`"
absolute_path_of_source="`cd "$up_path" ; pwd`"
relative_path_of_cwd="${called_from##$absolute_path_of_source}"

########################################################################
# Set some default variables which can be overridden from local.sh

# The default behaviour
# build_tree="/tmp/cado-build-`id -un`"
: ${build_tree:="${up_path}build/`hostname`"}

# Joe gets confused by " in prev line; this line fixes it

# By default, we also avoid /usr/local ; of course, it may be overridden.
# Note that cmake does not really have phony targets, so one must avoid
# basenames which correspond to targets !
: ${PREFIX:="$absolute_path_of_source/installed"}

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
: ${PTHREADS:=1}

########################################################################
# Arrange so that relevant stuff is passed to cmake -- the other end of
# the magic is in CMakeLists.txt. The two lists must agree.
# (no, it's not as simple. As long as the cmake checks care about
# *environment variables*, we are here at the right place for setting
# them. However, we might also be interested in having cmake export test
# results to scripts. This is done by cmake substitutions, but the
# corresponding names need not match the ones below).

export CC
export CXX
export CFLAGS
export CXXFLAGS
export FLAGS_SIZE
export LDFLAGS
export PREFIX
export CADO_DIST_ARCHIVE_NAME
export MPI
export GF2X_CONFIGURE_EXTRA
export GMP
export GMP_INCDIR
export GMP_LIBDIR
export MPIR
export MPIR_INCDIR
export MPIR_LIBDIR
export PTHREADS
export CURL
export CURL_INCDIR
export CURL_LIBDIR
export HWLOC
export HWLOC_INCDIR
export HWLOC_LIBDIR
export NUMA
export NUMA_INCDIR
export NUMA_LIBDIR
export GF2X_CONFIGURE_EXTRA_FLAGS
export CMAKE_DUMP_VARIABLES
export ENABLE_SHARED
export NO_PYTHON_CHECK
export NO_SSE
export NO_INLINE_ASSEMBLY
export CHECKS_EXPENSIVE
export BWC_GF2_MATMUL_BACKENDS
export BWC_GFP_MATMUL_BACKENDS
export BWC_GF2_ARITHMETIC_BACKENDS
export BWC_GFP_ARITHMETIC_BACKENDS
export BWC_EXTRA_BACKENDS

########################################################################
# "make ?" or "make help" (when the Makefile does not exist)
info=""
function make_usage {
    echo "-------------------------------------------------------------- "
    echo "[Options] (see $0 for details)"
    echo "-------------------------------------------------------------- "
    echo " $info \"make ?\"       -- this help"
    echo " $info \"make show\"    -- only show env variables"
    echo " $info \"make cmake\"   -- only run cmake to generate Makefile"
    echo " $info \"make\"         -- run cmake first and then make"
    echo " $info \"make tidy\"    -- delete folder $build_tree (dangerous)"
    echo " $info  Any other options will be passed to the actual make followed."
    echo "-------------------------------------------------------------- "
    exit 0
}
if [ "$1" == "?" ] ; then
    make_usage
fi
if [ "$1" == "help" ] && [ ! -f "$build_tree/Makefile" ] ; then
    make_usage
fi

########################################################################
# make tidy (warn, this delete the whole build folder)
wn="[warning]"
if [ "$1" == "tidy" ] ; then
    echo "$wn this deletes the whole folder $build_tree. Do you want to continue? (n/Y)"
    if [ -e "`tty`" ] ; then
        read TIDY_BUILD
    else
        echo "$wn no input terminal, assuming no"
        TIDY_BUILD=n
    fi
    if [ "$TIDY_BUILD" == "Y" ]; then
        echo "$wn wiping out $build_tree"
        rm -rf "$build_tree"
    else
        echo "$wn no action and quit now"
    fi
    exit 0
fi


########################################################################
# make show
if [ "$1" = "show" ] ; then
    echo "build_tree=\"$build_tree\""
    echo "up_path=\"$up_path\""
    echo "called_from=\"$called_from\""
    echo "absolute_path_of_source=\"$absolute_path_of_source\""
    echo "relative_path_of_cwd=\"$relative_path_of_cwd\""
    echo "CC=\"$CC\""
    echo "CXX=\"$CXX\""
    echo "CFLAGS=\"$CFLAGS\""
    echo "CXXFLAGS=\"$CXXFLAGS\""
    echo "FLAGS_SIZE=\"$FLAGS_SIZE\""
    echo "LDFLAGS=\"$CFLAGS\""
    echo "PREFIX=\"$PREFIX\""
    echo "CADO_DIST_ARCHIVE_NAME=\"$CADO_DIST_ARCHIVE_NAME\""
    echo "MPI=\"$MPI\""
    echo "GF2X_CONFIGURE_EXTRA=\"$GF2X_CONFIGURE_EXTRA\""
    echo "PTHREADS=\"$PTHREADS\""
    echo "CMAKE_GENERATOR=\"$CMAKE_GENERATOR\""
    echo "ENABLE_SHARED=\"$ENABLE_SHARED\""
    echo "CHECKS_EXPENSIVE=\"$CHECKS_EXPENSIVE\""
    echo "BWC_GF2_MATMUL_BACKENDS=\"$BWC_GF2_MATMUL_BACKENDS\""
    echo "BWC_GFP_MATMUL_BACKENDS=\"$BWC_GFP_MATMUL_BACKENDS\""
    echo "BWC_GF2_ARITHMETIC_BACKENDS=\"$BWC_GF2_ARITHMETIC_BACKENDS\""
    echo "BWC_GFP_ARITHMETIC_BACKENDS=\"$BWC_GFP_ARITHMETIC_BACKENDS\""
    echo "BWC_EXTRA_BACKENDS=\"$BWC_EXTRA_BACKENDS\""
    exit 0
fi

########################################################################
# Make sure we have cmake, by the way !
:  ${cmake_path:="`which cmake 2>/dev/null`"}
cmake_companion_install_location="$absolute_path_of_source/cmake-installed"
if [ "$?" != "0" ] || ! [ -x "$cmake_path" ] ; then
    echo "CMake not found" >&2
    cmake_path=
# Recall that (some versions of) bash do not want quoting for regex patterns.
elif [[ "`"$cmake_path" --version`" =~ ^cmake\ version\ [012] ]] && ! [[ "`"$cmake_path" --version`" =~ ^cmake\ version\ 2.[89] ]] ; then
    echo "CMake found, but not with version 2.8 or newer" >&2
    cmake_path=
fi
if ! [ "$cmake_path" ] ; then
    cmake_path="$cmake_companion_install_location/bin/cmake"
    if [ -x "$cmake_path" ] ; then
        echo "Using custom cmake in $cmake_companion_install_location" >&2
    else
        echo "I am about to download and compile a compatible version of Cmake."
        echo "Do you want to continue ? (y/n)"
        if [ -e "`tty`" ] ; then
            read INSTALL_CMAKE
        else
            echo "No input terminal, assuming yes"
            INSTALL_CMAKE=y
        fi
        if [ ! "$INSTALL_CMAKE" = "y" ]; then
            echo "Please install a compatible version of Cmake."
            exit 1
        fi
        echo "Need to get cmake first -- this takes long !"
        cd "$up_path"
        if ! scripts/install-cmake.sh "$cmake_companion_install_location" ; then
            echo "cmake install Failed, sorry" >&2
            exit 1
        fi
        cd "$called_from"
    fi
fi


########################################################################
# handle "make clean"
if [ "$1" == "clean" ] && [ ! -f "$build_tree/Makefile" ] ; then
    echo "There is no $build_tree/Makefile. Nothing to clean."
    exit 0
fi

########################################################################
# call cmake (if Makefile does not exist)
if [ "$1" = "cmake" ] || [ ! -f "$build_tree/Makefile" ] ; then
    mkdir -p "$build_tree"
    absolute_path_of_build_tree="`cd "$build_tree" ; $pwdP`"
    if [ ! "x$CMAKE_GENERATOR" == "x" ] ; then
      CMAKE_GENERATOR_OPT="-G$CMAKE_GENERATOR"
    else
      unset CMAKE_GENERATOR_OPT
    fi
    (cd "$absolute_path_of_build_tree" ; "$cmake_path" "$CMAKE_GENERATOR_OPT" "$absolute_path_of_source")
fi

if [ "$1" = "cmake" ] ; then
    exit 0
fi

########################################################################
# Now cd into the target directory, and build everything required.
# Note that it's useful to kill MAKELEVEL, or otherwise we'll get scores
# and scores of ``Entering directory'' messages (sure, there's the
# --no-print-directory option -- but it's not the right cure here).
# env | grep -i make
unset MAKELEVEL
absolute_path_of_build_tree="`cd "$build_tree" ; $pwdP`"
(cd "$absolute_path_of_build_tree$relative_path_of_cwd" ; ${MAKE} "$@")
