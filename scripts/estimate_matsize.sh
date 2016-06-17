#!/bin/bash
#
# The purpose of this script is to help predict the size of the matrix
# that one would obtain with a given set of parameters.
# It takes as input a polynomial file, and sieving / filtering parameters
# given as env variables.
# The CADO_BUILD variable must point to the directory where cado-nfs was
# compiled.
#
# Typical usage will be something like:
#   wdir=/tmp/test I=12 lim0=... ./estimate_matsize.sh toto.poly
# The root files and the renumber file are cached, so if you change lim or
# lpb without changing wdir, you'll have to remove those cache-files
# first in wdir.
#
# WARNING: this works only if we sieve only on one side.

if [ $# != 1 ]; then
    echo "Usage: $0 <polyfile>"
    exit 1
fi

## default parameters: can be overriden using env variables
: ${I=15}
: ${lim0=30000000}
: ${lim1=30000000}
: ${lpb0=27}
: ${lpb1=27}
: ${mfb0=54}
: ${mfb1=54}
: ${qmax=100000000}
: ${sqside=1}
: ${maxlevel=25}
: ${target_density=150}

## read poly file on command line
polyfile=$1
if [ ! -e $1 ]; then
    echo "Error: file $1 does not exist?"
    exit 1
fi

## if wdir is not set, build one in /tmp
if [ -z ${wdir+x} ]; then
    wdir=`mktemp -d /tmp/est_mat.XXXXXXX`;
else
    mkdir -p $wdir || (echo "mkdir -p $wdir failed" && exit 1)
fi
echo "Working directory is $wdir"

## if wdir does not contain a rootfile, build it
rootfile0="$wdir/roots0.gz"
rootfile1="$wdir/roots1.gz"
if [ ! -e $rootfile0 ]; then
    echo "Root file for side 0 not in wdir: building it..."
    $CADO_BUILD/sieve/makefb -poly $polyfile -lim $lim0 -side 0 -t 2 -out $rootfile0
fi
if [ ! -e $rootfile1 ]; then
    echo "Root file for side 1 not in wdir: building it..."
    $CADO_BUILD/sieve/makefb -poly $polyfile -lim $lim1 -side 1 -t 2 -out $rootfile1
fi

## if wdir does not contain a renumber table, build it
renumberfile=$wdir/renumber.gz
if [ ! -e $renumberfile ]; then
    echo "Renumber file not in wdir: building it..."
    $CADO_BUILD/sieve/freerel -poly $polyfile -renumber $renumberfile \
      -out /dev/null -pmax 1 -lpb0 $lpb0 -lpb1 $lpb1 -t 2
fi

## split the qrange into 5 chunks
NCHUNKS=5
if [ $sqside == "0" ]; then
    qmin=$lim0
    qqmax=`echo "2 ^ $lpb0" | bc`
else
    qmin=$lim1
    qqmax=`echo "2 ^ $lpb1" | bc`
fi
if [ $qqmax -le $qmax ]; then
    echo "qmax should be less then lpb"
    exit 1
fi
qrange=$(((qmax-qmin)/NCHUNKS))

## sample with real sieving and build fake rels
NBSAMPLE=50
fakefiles=()
for i in `seq 1 $NCHUNKS`; do
    q0=$((qmin + (i-1)*qrange))
    q1=$((qmin + i*qrange))
    echo "Dealing with qrange=[$q0,$q1]"
    echo "  Sampling..."
    $CADO_BUILD/sieve/las -I $I -poly $polyfile -q0 $q0 -q1 $q1 \
      -lim0 $lim0 -lim1 $lim1 -lpb0 $lpb0 -lpb1 $lpb1 -sqside $sqside \
      -mfb0 $mfb0 -mfb1 $mfb1 \
      -fb0 $rootfile0 -fb1 $rootfile1 -random-sample $NBSAMPLE \
      -t 1 -dup > $wdir/sample.${q0}-${q1}
    echo "  Building fake relations..."
    $CADO_BUILD/sieve/fake_rels -poly $polyfile -lpb0 $lpb0 -lpb1 $lpb1 \
      -q0 $q0 -q1 $q1 -sqside $sqside -sample $wdir/sample.${q0}-${q1} \
      -renumber $renumberfile > $wdir/fakerels.${q0}-${q1}
    fakefiles+=("$wdir/fakerels.${q0}-${q1}")
done
nrels=`cat ${fakefiles[@]} | grep -v "^#" | wc -l`
echo "We have $nrels fake relations"

## Filtering
# take a huge upper bound for the number of primes...
nprimes=`echo "2*2^$lpb0/l(2^$lpb0) + 2*2^$lpb1/l(2^$lpb1)" | bc -l | cut -d "." -f 1`

# purge
$CADO_BUILD/filter/purge -out $wdir/purged.gz -nrels $nrels -keep 3 \
    -col-min-index 10000 -col-max-index $nprimes -t 2 ${fakefiles[@]} \
    2>&1 | tee $wdir/purge.log

# Did we get a positive excess?
if (grep "number of rows < number of columns + keep" $wdir/purge.log > /dev/null); then
    echo "Negative excess: no way to build a matrix!"
    echo "You'll have to change the parameters"
    exit 1
fi

# merge
$CADO_BUILD/filter/merge-dl -mat $wdir/purged.gz -out $wdir/history.gz \
    -maxlevel $maxlevel -keep 3 -skip 0 -target_density $target_density
