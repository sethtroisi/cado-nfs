#!/usr/bin/env bash

N=100
dens=10
seed=1
bindir=

set -e

usage() {
    echo "Usage: $0 <N>" >&2
    exit 1
}

while [ $# -gt 0 ] ; do
    if [ "$1" = "--matrix-size" ] ; then
        shift
        N=$1
        shift
    elif [ "$1" = "--density" ] ; then
        shift
        dens=$1
        shift
    elif [ "$1" = "--seed" ] ; then
        shift
        seed=$1
        shift
    elif [ "$1" = "--bindir" ] ; then
        shift
        bindir=$1
        shift
    else
        usage
    fi
done

if ! [ "$N" ] ; then usage ; fi

wdir=$(mktemp -d  /tmp/cado.XXXXXXXX)

cleanup() { rm -f $wdir ; }
trap cleanup EXIT

$bindir/random_matrix $N -d 10  --binary -o $wdir/mat.bin -s $seed
$bindir/mf_scan mfile=$wdir/mat.bin --binary-in --freq

eval $($bindir/mf_bal  1 1 $wdir/mat.bin skip_decorrelating_permutation=true 2>&1 | tee /dev/stderr | perl -ne 's/Writing.*to\s/B=/ && print;')

"`dirname $0`"/perlrandom.pl $((8*N)) $seed  > $wdir/Y.0

$bindir/spmv_test wdir=$wdir mn=64 prime=2 balancing=$B matrix=$wdir/mat.bin nullspace=left
mv $wdir/MY.0 $wdir/YM.0

$bindir/short_matmul  $wdir/mat.bin  $wdir/Y.0  > $wdir/sMY.0

$bindir/short_matmul -t $wdir/mat.bin  $wdir/Y.0  > $wdir/sYM.0

$bindir/spmv_test wdir=$wdir mn=64 prime=2 balancing=$B matrix=$wdir/mat.bin nullspace=RIGHT

diff -q $wdir/YM.0 $wdir/sYM.0
diff -q $wdir/MY.0 $wdir/sMY.0

rm -rf $wdir


