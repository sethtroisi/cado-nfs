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
#
# Pay attention to the fact that when the implementation layers expect the
# matrix in column major order, we must not use --kleft, but rather
# --kright, even though we use -t in the end.
# $bins/random  90 100 -c 10 --kright 10 > $mats/t100p.txt
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
echo "var:=Transpose(var);" >> $mdir/t.m


$cmd spvector32 < /tmp/rowvec1 > $mdir/rowvec1.m
$cmd spvector32 < /tmp/rowvec2 > $mdir/rowvec2.m
$cmd spvector32 < /tmp/colvec1 > $mdir/colvec1.m
$cmd spvector32 < /tmp/colvec2 > $mdir/colvec2.m
