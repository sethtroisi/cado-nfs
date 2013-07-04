#!/bin/bash

set -e
set -x

top=`dirname $0`/../..
export DEBUG=1
# export MPI=1
make -s -C $top -j 4
make -s -C $top -j 4 plingen
eval `make -s -C $top show`
bins=$top/$build_tree/linalg/bwc
mats=$HOME/Local/mats
if ! [ -d $mats ] ; then
    mats=/local/rsa768/mats
fi

wdir=/tmp/bwcp
if [ -d $wdir ] ; then rm -rf $wdir 2>/dev/null ; fi
mkdir $wdir

m=8 n=4

Mh=1; Mv=1;
Th=1; Tv=1;
Nh=$((Mh*Th))
Nv=$((Mv*Tv))

mpi=${Mh}x${Mv}
thr=${Th}x${Tv}

prime=1766847064778384329583297500742918515827483896875618958121606201292619891
bits_per_coeff=256

if ! echo "$prime-2^$bits_per_coeff" | bc | grep -q '^-' ; then
    echo "2^bits_per_coeff must be above the prime" >&2
    echo "Please fix $0" >&2
    exit 1
elif echo "$prime-2^($bits_per_coeff-64)" | bc | grep -q '^-' ; then
    echo "2^(bits_per_coeff-64) must be below the prime" >&2
    echo "Please fix $0" >&2
    exit 1
fi

# The test matrix may be created by:
#
# Pay attention to the fact that when the implementation layers expect the
# matrix in column major order, we must not use --kleft, but rather
# --kright, even though we use -t in the end.

# -c 10 imposes a bound on the coefficients.

$bins/random_matrix  1000 -c 10 --kright 10 > $mats/t1000p.txt
matrix_txt=$mats/t1000p.txt
matrix=$mats/t1000p.bin

$bins/mf_scan  --ascii-in --withcoeffs --mfile $matrix_txt  --freq --binary-out --ofile $matrix

nullspace=right
interval=50



# Note that it's better to look for a kernel which is not trivial. Thus
# specifying --kright for random generation is a good move prior to
# running this script for nullspace=right

if [ "$shuffle" = 1 ] ; then
    shuffle_option=--shuffled-product
else
    # I have the impression that it's better to leave this here !
    shuffle_option="noshuffle=1"
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

common="matrix=$matrix mpi=$mpi thr=$thr balancing=$bfile m=$m n=$n wdir=$wdir prime=$prime mm_impl=basicp"
if [ "$nullspace" = left ] ; then
    common="$common nullspace=left"
    transpose_if_left="Transpose"
else
    common="$common nullspace=right"
    transpose_if_left=""
fi

set -e

all_splits=0
j0=0
while [ $j0 -lt $n ] ; do
    let j0=$j0+1
    all_splits=$all_splits,$j0
done


# ys=0..1 here is really a hack. It merely has to match the version
# which is used in production.
$bins/bwc.pl dispatch $common save_submatrices=1 ys=0..1

if [ $n -eq 1 ] ; then
    $bins/bwc.pl prep   $common
    ln -s Y.0 $wdir/V0-1.0
else
    # Otherwise prep won't work. Let's be stupid.
    set `head -1 $matrix_txt`
    ncols=$2
    nbytes=$((bits_per_coeff/8 * $2))
    j0=0
    while [ $j0 -lt $n ] ; do
        let j1=$j0+1
        echo "Generating $wdir/V${j0}-${j1}.0"
        dd if=/dev/urandom bs=1 count=$nbytes of=$wdir/V${j0}-${j1}.0
        j0=$j1
    done
    (echo 1 ; seq 0 $((m-1))) > $wdir/X
    # TODO: create Y, too.
fi

$bins/bwc.pl secure  $common interval=$interval
$bins/bwc.pl secure  $common interval=1
# [ "$?" = 0 ] && $bins/bwc.pl :ysplit $common splits=$all_splits
# [ "$?" = 0 ] && $bins/bwc.pl krylov  $common interval=$interval end=$interval ys=0..1 skip_online_checks=1
$bins/bwc.pl krylov  $common interval=1 end=$interval ys=0..1 skip_online_checks=1
rm -f $wdir/A*
j0=0
while [ $j0 -lt $n ] ; do
    let j1=$j0+1
    $bins/bwc.pl krylov  $common interval=$interval ys=$j0..$j1
    j0=$j1
done

afile=$($bins/acollect wdir=$wdir m=$m n=$n bits-per-coeff=$bits_per_coeff --remove-old | tail -1)
$bins/plingen lingen-mpi-threshold=10000 lingen-threshold=10 m=$m n=$n wdir=$wdir prime=$prime afile=$afile

ln $wdir/$afile.gen $wdir/F0-$n
j0=0
while [ $j0 -lt $n ] ; do
    let j1=$j0+1
    $bins/bwc.pl mksol  $common interval=$interval ys=$j0..$j1
    j0=$j1
done

$bins/bwc.pl gather  $common interval=$interval || :


# [ "$?" = 0 ] && $bins/bwc.pl acollect    $common -- --remove-old
# 
# 
# 
# # set +e
# # $bins/bench  --nmax 1000 --prime $prime --nbys 1 --impl basicp  -- $matrix
# 
mdir=$wdir

split=${Nh}x${Nv}
b=`basename $matrix .bin`
c=$b.$split.$checksum

cmd=`dirname $0`/convert_magma.pl

rwfile=${matrix%%bin}rw.bin
cwfile=${matrix%%bin}cw.bin
$cmd weights < $rwfile > $mdir/rw.m
$cmd weights < $cwfile > $mdir/cw.m

# This is for the **unbalanced** matrix !!

$cmd bpmatrix < $matrix > $mdir/t.m

$cmd balancing < $wdir/$c.bin > $mdir/b.m

(
echo "nullspace:=\"$nullspace\";"
echo "xtr:=func<x|$transpose_if_left(x)>;"
echo "nh:=$Nh;"
echo "nv:=$Nv;"
echo "nrp:=nv*Ceiling(nr/(nh*nv));"
echo "ncp:=nh*Ceiling(nc/(nh*nv));"
echo "Mt:=Matrix(GF($prime),nh*nrp,nv*ncp,[]);"
for i in `seq 0 $((Nh-1))` ; do
    echo "nr$i:=nrp;" # nr div nh + ($i lt nr mod nh select 1 else 0);"
    echo "snr$i:=$i *nrp;" # $i*(nr div nh) + Min($i, nr mod nh);"
    for j in `seq 0 $((Nv-1))` ; do
        $cmd bpmatrix < $wdir/$c.h$i.v$j> $mdir/t$i$j.m
        echo "nc$j:=ncp;" #  div nv + ($j lt nc mod nv select 1 else 0);"
        echo "snc$j:=$j*ncp;" # (nc div nv) + Min($j, nc mod nv);"
        echo "load \"$mdir/t$i$j.m\"; M$i$j:=Matrix(GF($prime),Matrix(var));"
        echo "x:=RMatrixSpace(GF($prime),nr$i,nc$j)!0;InsertBlock(~x,$transpose_if_left(M$i$j),1,1);M$i$j:=x;"
        echo "InsertBlock(~Mt,M$i$j,1+snr$i,1+snc$j);"
    done
done
echo "mlist:=["
for i in `seq 0 $((Nh-1))` ; do
    if [ "$i" != 0 ] ; then echo "," ; fi
    echo -n "["
    for j in `seq 0 $((Nv-1))` ; do
        if [ "$j" != 0 ] ; then echo -n ", " ; fi
        echo -n "M$i$j";
    done
    echo -n "]"
done
echo "];"
) > $mdir/placemats.m


$cmd spvector64 < $wdir/Y.0 > $mdir/Y0.m
$cmd spvector64 < $wdir/V0-1.0 > $mdir/V0.m
$cmd spvector64 < $wdir/V0-1.1 > $mdir/V1.m
$cmd spvector64 < $wdir/V0-1.2 > $mdir/V2.m
$cmd spvector64 < $wdir/V0-1.3 > $mdir/V3.m
$cmd spvector64 < $wdir/V0-1.4 > $mdir/V4.m
$cmd spvector64 < $wdir/V0-1.5 > $mdir/V5.m
$cmd spvector64 < $wdir/V0-1.6 > $mdir/V6.m
$cmd spvector64 < $wdir/V0-1.7 > $mdir/V7.m
$cmd spvector64 < $wdir/V0-1.8 > $mdir/V8.m
$cmd spvector64 < $wdir/V0-1.9 > $mdir/V9.m
$cmd spvector64 < $wdir/V0-1.10 > $mdir/V10.m
$cmd spvector64 < $wdir/C.0 > $mdir/C0.m
$cmd spvector64 < $wdir/C.1 > $mdir/C1.m
$cmd spvector64 < $wdir/C.$interval > $mdir/C$interval.m
$cmd x $wdir/X > $mdir/x.m
$cmd spvector64 < $wdir/$afile > $mdir/A.m
$cmd spvector64 < $wdir/F0-$n > $mdir/F0-$n.m
$cmd spvector64 < $wdir/W > $mdir/W.m
