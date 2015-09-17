#!/usr/bin/env bash

nrows=100
ncols=100
density=10
seed=1
bindir=
bwc_extra=()
mf_bal_extra=()
nh=1
nv=1

set -e

if [ "$CADO_DEBUG" ] ; then set -x ; fi

usage() {
    echo "Usage: $0 <N>" >&2
    exit 1
}

while [ $# -gt 0 ] ; do
    if [ "$1" = "--matrix-size" ] ; then
        shift
        nrows=`echo $1 | cut -d, -f1`
        ncols=`echo $1 | cut -d, -f2`   # = nrows if no comma
        shift
    elif [ "$1" = "--density" ] ; then
        shift
        density=$1
        shift
    elif [ "$1" = "--seed" ] ; then
        shift
        seed=$1
        shift
    elif [ "$1" = "--bindir" ] ; then
        shift
        bindir=$1
        shift
    elif [ "$1" = "--bindir" ] ; then
        shift
        bindir=$1
        shift
    else
        case "$1" in
            thr=*|mpi=*)
                x=`echo $1 | cut -d= -f2 | cut -dx -f1`
                nh=$((x*nh))
                x=`echo $1 | cut -d= -f2 | cut -dx -f2`
                nv=$((x*nv))
                bwc_extra=("${bwc_extra[@]}" "$1")
                shift
                ;;
            *) usage;;
        esac
    fi
done

redirect_unless_debug() {
    file="$1"
    shift
    if [ "$CADO_DEBUG" ] ; then
        "$@" 2>&1 | tee "$file"
    else
        "$@" > $file 2>&1
    fi
}

if ! [ "$nrows" ] ; then usage ; fi

wdir=$(mktemp -d  /tmp/cado.XXXXXXXX)

cleanup() { if ! [ "$CADO_DEBUG" ] ; then rm -rf $wdir ; fi ; }
trap cleanup EXIT

redirect_unless_debug $wdir/random_matrix.out $bindir/linalg/bwc/random_matrix $nrows $ncols -d $density  --binary -o $wdir/mat.bin -s $seed
redirect_unless_debug $wdir/scan.out $bindir/linalg/bwc/mf_scan mfile=$wdir/mat.bin --binary-in --freq

if [ $nrows != $ncols ] ; then
    mf_bal_extra=(--rectangular --reorder both)
fi
redirect_unless_debug $wdir/bal.out $bindir/linalg/bwc/mf_bal $nh $nv $wdir/mat.bin skip_decorrelating_permutation=true out=$wdir/bal.bin "${mf_bal_extra[@]}"
B=$wdir/bal.bin

"`dirname $0`"/perlrandom.pl $((8*nrows)) $seed  > $wdir/Y.0
$bindir/tests/linalg/bwc/short_matmul -t $wdir/mat.bin  $wdir/Y.0  $wdir/sYM.0 > /dev/null 2>&1
for impl in basic sliced bucket threaded ; do
    redirect_unless_debug $wdir/spmv-$impl-left.out $bindir/linalg/bwc/bwc.pl :mpirun ${bwc_extra} -- $bindir/tests/linalg/bwc/spmv_test wdir=$wdir mn=64 prime=2 balancing=$B matrix=$wdir/mat.bin nullspace=left mm_impl=$impl no_save_cache=1 "${bwc_extra[@]}"
    mv $wdir/MY.0 $wdir/YM.0
    diff -q $wdir/YM.0 $wdir/sYM.0
    echo "spmv ${nrows}x${ncols} $impl left ok ${bwc_extra[@]}"
done

"`dirname $0`"/perlrandom.pl $((8*ncols)) $seed  > $wdir/Y.0
$bindir/tests/linalg/bwc/short_matmul  $wdir/mat.bin  $wdir/Y.0  $wdir/sMY.0 >/dev/null 2>&1
for impl in basic sliced bucket threaded ; do
    redirect_unless_debug $wdir/spmv-$impl-right $bindir/linalg/bwc/bwc.pl :mpirun ${bwc_extra} -- $bindir/tests/linalg/bwc/spmv_test wdir=$wdir mn=64 prime=2 balancing=$B matrix=$wdir/mat.bin nullspace=RIGHT mm_impl=$impl no_save_cache=1 "${bwc_extra[@]}"
    diff -q $wdir/MY.0 $wdir/sMY.0
    echo "spmv ${nrows}x${ncols} $impl right ok ${bwc_extra[@]}"
done

