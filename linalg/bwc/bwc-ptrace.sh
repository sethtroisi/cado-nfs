#!/bin/bash

set -e
set -x

top=`dirname $0`/../..
export DEBUG=1
# export MPI=1
make -s -C $top -j 4
eval `make -s -C $top show`
bins=$top/$build_tree/linalg/bwc
mats=$HOME/Local/mats
if ! [ -d $mats ] ; then
    mats=/local/rsa768/mats
fi

wdir=/tmp/bwc
if [ -d $wdir ] ; then rm -rf $wdir 2>/dev/null ; fi
mkdir $wdir

mn=1

Mh=1; Mv=1;
Th=2; Tv=2;
Nh=$((Mh*Th))
Nv=$((Mv*Tv))

mpi=${Mh}x${Mv}
thr=${Th}x${Tv}

# The test matrix may be created by:
# $bins/random  100 -c 10 --kleft 10 > $mats/t100p.txt
# $bins/mf_scan  --ascii-in --withcoeffs --mfile $mats/t100p.txt  --freq --binary-out --ofile $mats/t100p.bin

matrix=$mats/t100p.bin
nullspace=left

# I'm avoiding the balancing completely here.

# and also, I acknowledge there is a bug presently.
set +e
$bins/bench  -t --nmax 1000 --prime 65521 --nbys 1 --impl basicp  -- $matrix

mdir=$wdir

cmd=`dirname $0`/convert_magma.pl

rwfile=${matrix%%bin}rw.bin
cwfile=${matrix%%bin}cw.bin
$cmd weights < $rwfile > $mdir/rw.m
$cmd weights < $cwfile > $mdir/cw.m

$cmd bpmatrix < $matrix > $mdir/t.m

$cmd pvector32 < /tmp/Lsrc > $mdir/lsrc.m
$cmd pvector32 < /tmp/Lmul > $mdir/lmul.m
$cmd pvector32 < /tmp/Rsrc > $mdir/rsrc.m
$cmd pvector32 < /tmp/Rmul > $mdir/rmul.m
