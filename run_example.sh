#!/bin/bash

# Run cadofactor.pl with 2 jobs (or 1 job with 2 threads) on the current
# machine, to factor the number given by the specified example file. One may
# use a abbreviation like c59 or c79, provided the corresponding params.xxx
# file is found in @CMAKE_INSTALL_PREFIX@/@example_subdir@
#
# this prints OK or FAILED, and the exit status corresponds.
#
# All partial data are put in a temp dir that is removed at the end in
# case of success. The temp dir is left here in case of failure. If
# CADO_DEBUG is set, the temp dir is never removed.

# Mac OS X's mktemp requires a prefix.
# And we need $t to be an absolute pathname, or else ./new_cadofactor.pl
# will fail.


: ${t:=`mktemp -d /tmp/cado.XXXXXXXXXX`}

chmod 755 $t

file=$1
shift

if ! [ "$file" ] ; then
    echo "Usage: $0 <file>" >&2
    echo "(in some cases, provided params.XXX can be found, XXX is an acceptable argument)" >&2
    exit 1
fi

#########################################################################
# Set paths properly.
cado_prefix="@CMAKE_INSTALL_PREFIX@"
example_subdir="@example_subdir@"
mpiexec="@MPIEXEC@"

if [ -d "$cado_prefix/$example_subdir" ] ; then
    # We're called in the install tree.
    if [ -f "$cado_prefix/$example_subdir/params.$file" ] ; then
        file="$cado_prefix/$example_subdir/params.$file" 
    fi
    bindir="$cado_prefix/bin"
    cadofactor="$bindir/cadofactor.pl"
elif [ -x "@CADO_NFS_SOURCE_DIR@/cadofactor.pl" ] ; then
    # Otherwise we're called from the source tree (or we hope so)
    if [ -f "@CADO_NFS_SOURCE_DIR@/params/params.$file" ] ; then
        file="@CADO_NFS_SOURCE_DIR@/params/params.$file"
    fi
    cadofactor="@CADO_NFS_SOURCE_DIR@/cadofactor.pl"
    # Make the path absolute.
    bindir=$(cd "`dirname $0`" ; pwd)
else
    # Otherwise we're called from the source tree (or we hope so)
    call_cmake="`dirname $0`/scripts/call_cmake.sh"
    if ! [ -x "$call_cmake" ] ; then
        echo "I don't know where I am !" >&2
    fi
    eval `$call_cmake show`
    if [ -f "${up_path}params/params.$file" ] ; then
        file="${up_path}params/params.$file"
    fi
    cadofactor="${up_path}cadofactor.pl"
    # Make the path absolute.
    build_tree=`cd "$build_tree" ; pwd`
    bindir="$build_tree"
fi

if [ ! -f $file ] ; then
    echo "$file not found" >&2
    exit 1
fi

#########################################################################
echo "Testing factorization as given by $file in $t"

[ "$CADO_DEBUG" ] && echo "(debug mode, temporary files will be kept in $t)"

cp $file $t/param

$bindir/fixup_params_file.sh $t/param

# Sets the machine description file
if [ -z "$CADO_USEHOST" ] ; then
   host=localhost
else
   host=`hostname --short`
fi

cat > $t/mach_desc <<EOF
[local]
tmpdir=$t/tmp
bindir=$bindir
$host cores=2
EOF

$cadofactor $t/param bindir=$bindir machines=$t/mach_desc wdir=$t delay=60 sievenice=0 selectnice=0 logfile=$t/out "$@"

rc=$?

#########################################################################
# Check result, clean up the mess afterwards.
if [ "$rc" = 0 ] ; then
    echo OK
    if ! [ "$CADO_DEBUG" ] ; then
        rm -rf $t
    fi
else
    echo "FAILED ; data left in $t"
fi
exit $rc
