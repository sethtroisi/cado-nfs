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

wdir=/tmp/bwcp
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

# -c 10 imposes a bound on the coefficients.

$bins/random  90 100 -c 10 --kright 10 > $mats/t100p.txt
$bins/mf_scan  --ascii-in --withcoeffs --mfile $mats/t100p.txt  --freq --binary-out --ofile $mats/t100p.bin

matrix=$mats/t100p.bin
nullspace=right



# Note that it's better to look for a kernel which is not trivial. Thus
# specifying --kright for random generation is a good move prior to
# running this script for nullspace=right

if [ "$shuffle" = 1 ] ; then
    shuffle_option=--shuffled-product
fi

$bins/mf_bal $shuffle_option mfile=$matrix $Nh $Nv out=$wdir/ --withcoeffs

if [ $(ls $wdir/`basename $matrix .bin`.${Nh}x${Nv}.*.bin | wc -l) != 1 ] ; then
    echo "Weird -- should have only one balancing file as output." >&2
    exit 1
fi
bfile=$(ls $wdir/`basename $matrix .bin`.${Nh}x${Nv}.????????.bin)
echo "Using balancing file $bfile"
checksum=${bfile#$wdir/`basename $matrix .bin`.${Nh}x${Nv}.}
checksum=`basename $checksum .bin`

common="matrix=$matrix mpi=$mpi thr=$thr balancing=$bfile mn=$mn wdir=$wdir"
if [ "$nullspace" = left ] ; then
    common="$common nullspace=left"
    transpose_if_left="Transpose"
else
    common="$common nullspace=right"
    transpose_if_left=""
fi

set +e

all_splits=0
j0=0
while [ $j0 -lt $mn ] ; do
    let j0=$j0+1
    all_splits=$all_splits,$j0
done


# ys=0..64 here is really a hack. It merely has to match the version
# which is used in production.
$bins/bwc.pl dispatch sanity_check_vector=H1   $common save_submatrices=1 ys=0..64

set +e
$bins/bench  --nmax 1000 --prime 65521 --nbys 1 --impl basicp  -- $matrix

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
