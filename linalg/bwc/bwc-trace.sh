#!/bin/bash

set -e
set -x

top=`dirname $0`/../..
MPI=1 make -s -C $top -j 32
eval `MPI=1 make -s -C $top show`
bins=$top/$build_tree/linalg/bwc
mats=$HOME/Local/mats

wdir=/tmp/bwc
if [ -d $wdir ] ; then rm -rf $wdir 2>/dev/null ; fi
mkdir $wdir


Mh=3; Mv=4;
Th=1; Tv=1;
Nh=$((Mh*Th))
Nv=$((Mv*Tv))

mpi=${Mh}x${Mv}
thr=${Th}x${Tv}

# The test matrix may be created by:
# linalg/bwc/random 100 100 15 --kleft 4 > $mats/t100.txt
# linalg/bwc/mf_scan $mats/t100.txt ofile=$mats/t100.bin --binary-out --freq
#
# For the purpose of this complete test, a very tiny matrix is a good
# choice. Note though that the code may encounter random failures due to
# singular inputs. In particular, a matrix rank below 64 won't work.
matrix=$mats/t100b.bin
nullspace=left
shuffle=1

# Note that it's better to look for a kernel which is not trivial. Thus
# specifying --kleft for random generation is a good move prior to
# running this script for nullspace=left

if [ "$shuffle" = 1 ] ; then
    shuffle_option=--shuffled-product
fi

$bins/mf_bal $shuffle_option mfile=$matrix $Nh $Nv out=$wdir/

if [ $(ls $wdir/`basename $matrix .bin`.${Nh}x${Nv}.*.bin | wc -l) != 1 ] ; then
    echo "Weird -- should have only one balancing file as output." >&2
    exit 1
fi
bfile=$(ls $wdir/`basename $matrix .bin`.${Nh}x${Nv}.????????.bin)
echo "Using balancing file $bfile"
checksum=${bfile#$wdir/`basename $matrix .bin`.${Nh}x${Nv}.}
checksum=`basename $checksum .bin`

# $bins/mf_dobal --matrix /local/rsa768/mats/c72.bin mpi=2x3 /local/rsa768/mats/c72.2x3.6e90c700.bin

common="matrix=$matrix mpi=$mpi thr=$thr balancing=$bfile mn=64 wdir=$wdir"
if [ "$nullspace" = left ] ; then
    common="$common nullspace=left"
    transpose_if_left="Transpose";
else
    common="$common nullspace=right"
    transpose_if_left="";
fi

set +e

$bins/bwc.pl u64_dispatch         $common save_submatrices=1
[ "$?" = 0 ] && $bins/bwc.pl u64n_prep   $common
[ "$?" = 0 ] && $bins/bwc.pl u64_secure  $common interval=10
[ "$?" = 0 ] && $bins/bwc.pl :ysplit     $common
[ "$?" = 0 ] && $bins/bwc.pl u64_krylov  $common interval=1 end=10 skip_online_checks=1
[ "$?" = 0 ] && rm -f $wdir/A*
[ "$?" = 0 ] && $bins/bwc.pl u64_krylov  $common interval=10
[ "$?" = 0 ] && $bins/bwc.pl acollect    $common -- --remove-old
[ "$?" = 0 ] && $bins/bwc.pl lingen      $common lingen_threshold=64
[ "$?" = 0 ] && $bins/bwc.pl :fsplit     $common
[ "$?" = 0 ] && $bins/bwc.pl u64_mksol   $common interval=10 ys=0..64
[ "$?" = 0 ] && $bins/bwc.pl u64n_gather $common interval=10
[ "$?" = 0 ] && $bins/mf_twistvec --nullspace $nullspace --untwist $bfile $wdir/W.$checksum >  $wdir/W

failure="$?"

split=${Nh}x${Nv}
b=`basename $matrix .bin`
c=$b.$split.$checksum

mdir=$wdir

cmd=`dirname $0`/convert_magma.pl

$cmd balancing < $wdir/$c.bin > $mdir/b.m

$cmd bmatrix < $matrix > $mdir/t.m

# For nullspace=left, the matrices are transposed.
(
echo "xtr:=func<x|$transpose_if_left(x)>;";
echo "nh:=$Nh;"
echo "nv:=$Nv;"
echo "nrp:=nv*Ceiling(nr/(nh*nv));"
echo "ncp:=nh*Ceiling(nc/(nh*nv));"
echo "Mt:=Matrix(GF(2),nh*nrp,nv*ncp,[]);"
for i in `seq 0 $((Nh-1))` ; do
    echo "nr$i:=nrp;" # nr div nh + ($i lt nr mod nh select 1 else 0);"
    echo "snr$i:=$i *nrp;" # $i*(nr div nh) + Min($i, nr mod nh);"
    for j in `seq 0 $((Nv-1))` ; do
        $cmd bmatrix < $wdir/$c.h$j.v$i> $mdir/t$j$i.m
        echo "nc$j:=ncp;" #  div nv + ($j lt nc mod nv select 1 else 0);"
        echo "snc$j:=$j*ncp;" # (nc div nv) + Min($j, nc mod nv);"
        echo "load \"$mdir/t$j$i.m\"; M$j$i:=Matrix(GF(2),Matrix(var));"
        echo "x:=RMatrixSpace(GF(2),nr$i,nc$j)!0;InsertBlock(~x,$transpose_if_left(M$j$i),1,1);M$j$i:=x;"
        echo "InsertBlock(~Mt,M$j$i,1+snr$i,1+snc$j);"
    done
done) > $mdir/placemats.m

$cmd x $wdir/X.$checksum > $mdir/x.m
# ok, strictly speaking the A files have nothing to do with vectors. But
# the procedure sort of works as is anyway... Same holds for F.
for f in $wdir/[AVCYFSW]*.$checksum ; do $cmd vector < $f > $mdir/`basename $f .$checksum`.m ; done
for f in $wdir/[AF]*[0-9] ; do $cmd vector < $f > $mdir/`basename $f .$checksum`.m ; done
$cmd vector < $wdir/W > $mdir/Wu.m

if [ "$failure" != 0 ] ; then
    echo "ERROR: partial failure encountered"
fi
