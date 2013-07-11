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
    -py           - use Python cadofactor script as worker
EOF
}

# Set "t" to a temp directory, if not already set.
# The ":" shell built-in causes variable expansion but nothing else to happen
: ${t:=`mktemp -d /tmp/cado.XXXXXXXXXX`}
# Make temp dir writable by us and readable by everyone
chmod 755 $t

n=$1
shift
if [ ! "$(grep "^[ [:digit:] ]*$" <<< $n)" ]; then
  echo "Error, first parameter must be the number to factor" >&2
  usage
  exit 1
fi

ssh=false
python=false
cores=1
for ((i=1; i<=3; i=i+1)) ; do
  if [ "$1" = "-t" ]; then
    cores=$2
    if [ ! "$(grep "^[ [:digit:] ]*$" <<< $cores)" ]; then
      echo "Error, number of cores '$cores' is not an integer" >&2
      usage
      exit 1
    fi
    shift 2
  elif [ "$1" = "-ssh" ]; then
    ssh=true
    shift
  elif [ "$1" = "-py" ]; then
    python=true
    shift
  else
    break
  fi
done


#########################################################################
# Set paths properly.

# Try to find out from which location this script is running:
# the install directory (i.e., where 'make install' put it),
# the build directory (where 'cmake' put it, substituting the 
# CMake variables below), or in the source directory tree.
cado_prefix="@CMAKE_INSTALL_PREFIX@"
cado_source_dir="@CADO_NFS_SOURCE_DIR@"
cado_build_dir="@CADO_NFS_BINARY_DIR@"
example_subdir="@example_subdir@"
mpiexec="@MPIEXEC@"

if [ "$0" -ef "$cado_prefix/bin/factor.sh" ] ; then
    # We're called in the install tree.
    paramdir="$cado_prefix/$example_subdir"
    bindir="$cado_prefix/bin"
    if $python
    then
    	scriptpath="$bindir"
    	cadofactor="$scriptpath/cadofactor.py"
    else
    	cadofactor="$bindir/cadofactor.pl"
    fi
    cputime="$bindir/cpu_time.sh"
elif [ "$0" -ef "$cado_build_dir/factor.sh" ] ; then
    # We're called in the build tree.
    paramdir="$cado_source_dir/params"
    cputime="$cado_source_dir/scripts/cpu_time.sh"
    if $python
    then
    	scriptpath="$cado_source_dir/scripts/cadofactor"
        cadofactor="$scriptpath/cadofactor.py"
    else
        cadofactor="$cado_source_dir/cadofactor.pl"
    fi
    # Make the path absolute.
    bindir="$cado_build_dir"
elif [ -f "`dirname $0`/cado_config_h.in" ] ; then
    # Otherwise we're called from the source tree (or we hope so)
    srcdir=$(cd "`dirname $0`" ; pwd)
    call_cmake="`dirname $0`/scripts/call_cmake.sh"
    unset build_tree
    if [ -x "$call_cmake" ] ; then
      # This should set $build_tree:
      eval `cd $srcdir ; $call_cmake show`
    fi
    if ! [ -z "$build_tree" ] ; then
      paramdir="${srcdir}/params/"
      cputime="$srcdir/scripts/cpu_time.sh"
      if $python
      then
          scriptpath="${srcdir}/scripts/cadofactor"
          cadofactor="${scriptpath}/cadofactor.py"
      else
          cadofactor="${srcdir}/cadofactor.pl"
      fi
      # Make the path absolute.
      bindir=`cd "$build_tree" ; pwd`
    fi
fi
if [ -z "$cadofactor" ] ; then
    echo "I don't know where I am !" >&2
    # but I do care!
    exit 1
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

# This copy is not strictly necessary (it used to be when we used the
# fixup_params_file.sh script), but it might not be a bad idea to have the
# parameter file together with the rest in the working directory, which also
# allows for easy modification of parameters for the specific factorization.
cp $file $t/param

# Sets the machine description file
if [ -z "$CADO_USEHOST" ] ; then
   host=localhost
else
   host=`hostname --short`
fi

mkdir $t/tmp


if $python; then
  $cadofactor --old "$t/param" n=$n bindir="$bindir" \
  sieve_max_threads=$cores poly_max_threads=$cores nthchar=$cores \
  bwmt=$cores wdir="$t" slaves="$host" scriptpath="$scriptpath" \
  serveraddress=localhost "$@"
else
  cat > $t/mach_desc <<EOF
[local]
tmpdir=$t/tmp
bindir=$bindir
$host cores=$cores
EOF
  if ! $ssh; then
    $cadofactor params=$t/param n=$n bindir=$bindir parallel=0 \
    sieve_max_threads=$cores nthchar=$cores\
    bwmt=$cores wdir=$t sievenice=0 polsel_nice=0 logfile=$t/out "$@"
  else
    $cadofactor params=$t/param n=$n bindir=$bindir parallel=1 \
    machines=$t/mach_desc nthchar=$cores bwmt=$bwmt wdir=$t \
    sievenice=0 polsel_nice=0 logfile=$t/out "$@"
  fi
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
