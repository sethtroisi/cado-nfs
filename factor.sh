#!/usr/bin/env bash

# Run cadofactor.pl on the current machine. This script is recommended only
# for small factorizations (less than 100 digits). For larger numbers, you
# should use the cadofactor.pl script (see the documentation in that script).
#
# All partial data are put in a temp dir that is removed at the end in
# case of success. The temp dir is left here in case of failure. If
# CADO_DEBUG is set, the temp dir is never removed.

# Mac OS X's mktemp requires a prefix.
# And we need $t to be an absolute pathname, or else ./new_cadofactor.pl
# will fail.

usage() {
    cat >&2 <<EOF
Usage: $0 <integer> [options] [arguments passed to cadofactor.pl]
    <integer>     - integer to be factored. must be at least 57 digits,
                    and free of small prime factors [parameters are
                    optimized only for 85 digits and above].
options:
    -t <integer>  - numbers of cores to be used.
    -ssh          - use ssh (see README) for distributing the polynomial
                    selection and sieve steps on localhost. For broader
                    use on several machines, use the advanced script
                    cadofactor.pl
EOF
}

: ${t:=`mktemp -d /tmp/cado.XXXXXXXXXX`}
chmod 755 $t

n=$1
shift
if [ ! "$(grep "^[ [:digit:] ]*$" <<< $n)" ]; then
  usage
  exit 1
fi

ssh=0
cores=1
for ((i=1; i<=2; i=i+1)) ; do
  if [ "$1" = "-t" ]; then
    cores=$2
    if [ ! "$(grep "^[ [:digit:] ]*$" <<< $cores)" ]; then
      usage
      exit 1
    fi
    shift 2
  fi
  if [ "$1" = "-ssh" ]; then
    ssh=1
    shift
  fi
done


#########################################################################
# Set paths properly.
cado_prefix="@CMAKE_INSTALL_PREFIX@"
example_subdir="@example_subdir@"
mpiexec="@MPIEXEC@"

if [ "$0" -ef "$cado_prefix/bin/factor.sh" ] ; then
    # We're called in the install tree.
    paramdir="$cado_prefix/$example_subdir"
    bindir="$cado_prefix/bin"
    cadofactor="$bindir/cadofactor.pl"
    cputime="$bindir/cpu_time.sh"
elif [ "$0" -ef "@CADO_NFS_BINARY_DIR@/factor.sh" ] ; then
    # We're called in the build tree.
    paramdir="@CADO_NFS_SOURCE_DIR@/params"
    cputime="@CADO_NFS_SOURCE_DIR@/scripts/cpu_time.sh"
    cadofactor="@CADO_NFS_SOURCE_DIR@/cadofactor.pl"
    # Make the path absolute.
    bindir="@CADO_NFS_BINARY_DIR@"
elif [ -f "`dirname $0`/cado_config_h.in" ] ; then
    # Otherwise we're called from the source tree (or we hope so)
    srcdir=$(cd "`dirname $0`" ; pwd)
    call_cmake="`dirname $0`/scripts/call_cmake.sh"
    if ! [ -x "$call_cmake" ] ; then
        echo "I don't know where I am !" >&2
    fi
    eval `cd $srcdir ; $call_cmake show`
    paramdir="${srcdir}/params/"
    cputime="$srcdir/scripts/cpu_time.sh"
    cadofactor="${srcdir}/cadofactor.pl"
    # Make the path absolute.
    bindir=`cd "$build_tree" ; pwd`
else
    echo "I don't know where I am !" >&2
    # but I do care
fi

if ! [ -d "$paramdir" ] ; then
    echo "Parameter dir $paramdir not found." >&2 ; exit 1
elif ! [ -d "$bindir" ] ; then
    echo "Binary dir $bindir not found." >&2 ; exit 1
elif ! [ -x "$cadofactor" ] ; then
    echo "Script $cadofactor not found." >&2 ; exit 1
elif ! [ -x "$cputime" ] ; then
    echo "Script $cputime not found." >&2 ; exit 1
else
    # Ok, everything looks good.
    :
fi

size=${#n}

# Try n, n+1, n-1, n+2...
for ((i=1; i<=4; i=i+1)) ; do
  file="$paramdir/params.c$size"
  if [ -f $file ] ; then
    break
  fi
  size=`expr $size + \( 2 \* \( $i % 2 \) - 1 \) \* $i`
done

if [ ! -f $file ] ; then
    echo "no parameter file found for c${#n} (last tried was $file)" >&2
    usage
    exit 1
fi

#########################################################################
[ "$CADO_DEBUG" ] && echo "(debug mode, temporary files will be kept in $t)"

cp $file $t/param

# Sets the machine description file
if [ -z "$CADO_USEHOST" ] ; then
   host=localhost
else
   host=`hostname --short`
fi

mkdir $t/tmp

cat > $t/mach_desc <<EOF
[local]
tmpdir=$t/tmp
bindir=$bindir
$host cores=$cores
EOF

if [ $ssh -eq 0 ]; then
  $cadofactor params=$t/param n=$n bindir=$bindir parallel=0 \
  sieve_max_threads=$cores nthchar=$cores\
  bwmt=$cores wdir=$t sievenice=0 polsel_nice=0 logfile=$t/out "$@"
else
  $cadofactor params=$t/param n=$n bindir=$bindir parallel=1 \
  machines=$t/mach_desc nthchar=$cores bwmt=$bwmt wdir=$t \
  sievenice=0 polsel_nice=0 logfile=$t/out "$@"
fi

rc=$?

#########################################################################
# Check result, clean up the mess afterwards.
if [ "$rc" = 0 ] ; then
    echo OK
    # Print timings
    cd $t
    $cputime
    if ! [ "$CADO_DEBUG" ] ; then
        rm -rf $t
    fi
else
    echo "FAILED ; data left in $t"
fi
exit $rc
