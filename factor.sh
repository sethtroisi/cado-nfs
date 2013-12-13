#!/usr/bin/env bash

# Run cadofactor.pl or cadofactor.py on the current machine. This script is 
# recommended only for small factorizations (less than 100 digits). For 
# larger numbers, you should use the cadofactor.pl or cadofactor.py script
# (see the documentation in that script).
#
# All partial data are put in a temp dir that is removed at the end in
# case of success. The temp dir is left here in case of failure. If
# CADO_DEBUG is set, the temp dir is never removed.

# Mac OS X's mktemp requires a prefix.
# And we need $t to be an absolute pathname, or else ./cadofactor.pl
# will fail.

usage() {
    cat >&2 <<EOF
Usage: $0 <integer> [options] [arguments passed to cadofactor.py]
    <integer>     - integer to be factored. must be at least 57 digits,
                    and free of small prime factors [parameters are
                    optimized only for 85 digits and above]
options:
    -t <integer>  - numbers of cores to be used
    -pl           - use Perl cadofactor.pl script as worker
    -py           - use Python cadofactor.py script as worker (default)
    -s <integer>  - numbers of slave scripts to be used (only with -py)
    -h <host1>,<host2>,...,<hostn>
                  - comma-separated list of hosts to use for factorization.
                    With -s <n>, run <n> many scripts per host.
                    Only for -py. Default: localhost
    -ssh          - When using cadofactor.pl, use ssh (see README) for
                    distributing the polynomial selection and sieve steps 
                    on localhost. For broader use on several machines, use
                    the advanced script cadofactor.pl (or use cadofactor.py)
    -timeout <integer>
                  - abort computation after this many seconds

    If the shell environment variable CADO_DEBUG is set to any non-empty
    string, the working directory is not deleted even if the factorization
    succeeds.

    If CADO_USEHOST is set to any non-empty string, the Python script will
    listen on "*" even if all clients given with -h are on localhost.
EOF
}

# If a duration is specified, find the "timeout" binary and store its
# path and the duration in the TIMEOUT array. If no duration is specified,
# or if the timeout binary is not found, leave TIMEOUT unset.
# Returns 1 (false) if a timeout was requested and the binary was not found.
# Otherwise returns 0 (true).
# We use an array for TIMEOUT so that we can expand it to exactly 3 words on
# the resulting command line.
function find_timeout() {
    local TIMEOUT_DURATION="$1"
    test -z "$TIMEOUT_DURATION" && return 0
    # Maybe it's not a GNU system, but with GNU coreutils installed
    # with "g" prefix to all binaries
    for TIMEOUT_NAME in timeout gtimeout
    do
        TIMEOUT_BIN="`which timeout 2>/dev/null`"
        if [ -n "$TIMEOUT_BIN" ]
        then
            TIMEOUT[0]="$TIMEOUT_BIN"
            TIMEOUT[1]="--signal=SIGINT"
            TIMEOUT[2]="$TIMEOUT_DURATION"
            return 0
        fi
    done
    return 1
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
python=true
cores=1
slaves=1
verbose=false
while [ -n "$1" ] ; do
  if [ "$1" = "-t" ]; then
    cores=$2
    if [ ! "$(grep "^[ [:digit:] ]*$" <<< $cores)" ]; then
      echo "Error, number of cores '$cores' is not an integer" >&2
      usage
      exit 1
    fi
    shift 2
  elif [ "$1" = "-s" ]; then
    slaves=$2
    if [ ! "$(grep "^[ [:digit:] ]*$" <<< $slaves)" ]; then
      echo "Error, number of slaves '$slaves' is not an integer" >&2
      usage
      exit 1
    fi
    shift 2
  elif [ "$1" = "-h" ]; then
    hostnames="$2"
    shift 2
  elif [ "$1" = "-timeout" ]; then
    timeout="$2"
    if [ ! "$(grep "^[ [:digit:] ]*$" <<< $timeout)" ]; then
      echo "Error, number of seconds '$timeout' is not an integer" >&2
      usage
      exit 1
    fi
    if ! find_timeout $2
    then
      echo "timeout binary not found. Timeout disabled"
    fi
    shift 2
  elif [ "$1" = "-ssh" ]; then
    ssh=true
    shift
  elif [ "$1" = "-py" ]; then
    python=true
    shift
  elif [ "$1" = "-pl" ]; then
    python=false
    shift
  elif [ "$1" = "-v" ]; then
    verbose=true
    shift
  else
    break
  fi
done

if [ $slaves -ne 1 ] && ! $python
then
  echo "The -s parameter is used only for the Python script. Ignoring it."
fi

if [ -n "$hostnames" ] && ! $python
then
  echo "The -h parameter is used only for the Python script. Ignoring it."
fi

# Set default value
: ${hostnames:=localhost}

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
example_subdir_py="@example_subdir_py@"
mpiexec="@MPIEXEC@"

if [ "$0" -ef "$cado_prefix/bin/factor.sh" ] ; then
    # We're called in the install tree.
    bindir="$cado_prefix/bin"
    if $python
    then
    	scriptpath="$bindir"
    	cadofactor="$scriptpath/cadofactor.py"
    	paramdir="$cado_prefix/$example_subdir_py"
    else
    	cadofactor="$bindir/cadofactor.pl"
    	paramdir="$cado_prefix/$example_subdir"
    fi
    cputime="$bindir/cpu_time.sh"
elif [ "$0" -ef "$cado_build_dir/factor.sh" ] ; then
    # We're called in the build tree.
    cputime="$cado_source_dir/scripts/cpu_time.sh"
    if $python
    then
    	scriptpath="$cado_source_dir/scripts/cadofactor"
        paramdir="$cado_source_dir/params_py"
        cadofactor="$scriptpath/cadofactor.py"
    else
        cadofactor="$cado_source_dir/cadofactor.pl"
        paramdir="$cado_source_dir/params"
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
      cputime="$srcdir/scripts/cpu_time.sh"
      if $python
      then
          scriptpath="${srcdir}/scripts/cadofactor"
          cadofactor="${scriptpath}/cadofactor.py"
          paramdir="${srcdir}/params_py/"
      else
          cadofactor="${srcdir}/cadofactor.pl"
          paramdir="${srcdir}/params/"
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
elif ! $python && ! [ -x "$cputime" ] ; then
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
if [ -z "$CADO_USEHOST" ] && echo "$hostnames" | grep -E "^(localhost, *)*localhost$" > /dev/null
then
   host=localhost # For the Perl script
   server_address="server.address=localhost" # For the Python script
else
   host=`hostname --short`
   server_address=""
fi

mkdir $t/tmp


if $python; then
    # $PYTHON is there to expand a shell variable having that name, if
    # provided, in the case there is no python3 script in the path, or if
    # one which is named otherwise, or placed in a non-prority location
    # is the path, is preferred. If $PYTHON is empty, this is a no-op
  "${TIMEOUT[@]}" $PYTHON $cadofactor "$t/param" N=$n tasks.execpath="$bindir" \
  tasks.threads=$cores tasks.workdir="$t" slaves.hostnames="$hostnames" \
  slaves.nrclients=$slaves \
  slaves.scriptpath="$scriptpath" "$server_address" \
  slaves.basepath="$t/client/" "$@"
else
  cat > $t/mach_desc <<EOF
[local]
tmpdir=$t/tmp
bindir=$bindir
$host cores=$cores
EOF
  if $verbose; then
    VERBOSE="-v"
  fi
  if ! $ssh; then
    "${TIMEOUT[@]}" $cadofactor params=$t/param n=$n bindir=$bindir parallel=0 \
    sieve_max_threads=$cores nthchar=$cores\
    bwmt=$cores wdir=$t sievenice=0 polsel_nice=0 logfile=$t/out $VERBOSE "$@"
  else
    "${TIMEOUT[@]}" $cadofactor params=$t/param n=$n bindir=$bindir parallel=1 \
    machines=$t/mach_desc nthchar=$cores bwmt=$cores wdir=$t \
    sievenice=0 polsel_nice=0 logfile=$t/out $VERBOSE "$@"
  fi
fi

rc=$?

#########################################################################
# Check result, clean up the mess afterwards.
if [ "$rc" = 0 ] ; then
    echo OK
    if ! $python
    then
        # Print timings
        cd $t
        $cputime
    fi
    if ! [ "$CADO_DEBUG" ] ; then
        rm -rf $t
    fi
else
    echo "FAILED ; data left in $t"
fi
exit $rc
