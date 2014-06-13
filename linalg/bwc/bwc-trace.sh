#!/usr/bin/env bash

set -e
set -x

# various configuration variables. environment can be used to override them

: ${scriptpath=$0}
: ${m=128}
: ${n=128}
: ${wdir=/tmp/bwc}
: ${mm_impl=bucket}

export m n scriptpath wdir mm_impl prime=2
exec "`dirname $0`/bwc-ptrace.sh" "$@"

#
: ${Mh=1}
: ${Mv=1}
: ${Th=2}
: ${Tv=2}
# Set the "matrix" variable in order to work on a real matrix.
: ${matrix=}
# Set the "bindir" variable to use pre-built binaries (must point to the
# directories with the bwc binaries)
: ${bindir=}

# This is related to the random matrix which gets automatically created
# if no matrix was given on the cmdline.
: ${random_matrix_size=1000}
# there is no random_matrix_coeffs_per_row ; see random_matrix.c
: ${random_matrix_maxcoeff=10}
: ${random_matrix_minkernel=10}
: ${mats=$HOME/Local/mats}
: ${matsfallback=/local/rsa768/mats}


# XXX note that $wdir is wiped out by this script !
: ${buildopts="MPI=1 DEBUG=1"}
: ${shuffle=1}
: ${nomagma=}
#


Nh=$((Mh*Th))
Nv=$((Mv*Tv))
mpi=${Mh}x${Mv}
thr=${Th}x${Tv}



top=`dirname $0`/../..

# build the software. This is optional, really. The only thing is that
# it's so easy to mess with the location of the binaries that building it
# by ourselves looks like a reasonable thing to do, for this illustrative
# script.

if ! [ "$bindir" ] ; then
    eval "export $buildopts"
    make -s -C $top -j 4
    make -s -C $top -j 4 plingen_p_$nwords
    # make -s -C $top -j 4 plingen_pz
    eval `make -s -C $top show`
    bindir=$top/$build_tree/linalg/bwc
fi

if ! [ -d $mats ] ; then mats=$matsfallback; fi
if [ -d $wdir ] ; then rm -rf $wdir 2>/dev/null ; fi
mkdir $wdir
 
# Create the test matrix if needed.

if ! [ "$matrix" ] ; then
    $bindir/random_matrix  ${random_matrix_size} --kright ${random_matrix_minkernel} > $mats/t${random_matrix_size}p.txt
    matrix=$mats/t${random_matrix_size}p.txt
fi


# $bindir/random_matrix 1000 1000 15 --kleft 4 > $mats/t1000.txt
# $bindir/mf_scan $mats/t1000.txt ofile=$mats/t1000.bin --binary-out --freq
#
# For the purpose of this complete test, a very tiny matrix is a good
# choice. Note though that the code may encounter random failures due to
# singular inputs. In particular, a matrix rank below 64 won't work.
matrix=$mats/t1009.bin
# matrix=$HOME/Local/mats/c59.small.bin
nullspace=left
shuffle=1


# Note that it's better to look for a kernel which is not trivial. Thus
# specifying --kleft for random generation is a good move prior to
# running this script for nullspace=left

if [ "$shuffle" = 1 ] ; then
    shuffle_option=--shuffled-product
fi

$bindir/mf_bal $shuffle_option mfile=$matrix $Nh $Nv out=$wdir/

if [ $(ls $wdir/`basename $matrix .bin`.${Nh}x${Nv}.*.bin | wc -l) != 1 ] ; then
    echo "Weird -- should have only one balancing file as output." >&2
    exit 1
fi
bfile=$(ls $wdir/`basename $matrix .bin`.${Nh}x${Nv}.????????.bin)
echo "Using balancing file $bfile"
checksum=${bfile#$wdir/`basename $matrix .bin`.${Nh}x${Nv}.}
checksum=`basename $checksum .bin`

# $bindir/mf_dobal --matrix /local/rsa768/mats/c72.bin mpi=2x3 /local/rsa768/mats/c72.2x3.6e90c700.bin

common="matrix=$matrix mpi=$mpi thr=$thr balancing=$bfile m=$m n=$n wdir=$wdir"
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
while [ $j0 -lt $n ] ; do
    let j0=$j0+64
    all_splits=$all_splits,$j0
done


# ys=0..64 here is really a hack. It merely has to match the version
# which is used in production.
$bindir/bwc.pl dispatch sanity_check_vector=H1   $common save_submatrices=1 ys=0..64
[ "$?" = 0 ] && $bindir/bwc.pl prep   $common
[ "$?" = 0 ] && $bindir/bwc.pl secure  $common interval=10
[ "$?" = 0 ] && $bindir/bwc.pl :ysplit $common splits=$all_splits
[ "$?" = 0 ] && $bindir/bwc.pl krylov  $common interval=10 end=10 skip_online_checks=1 ys=0..64
[ "$?" = 0 ] && rm -f $wdir/A*
j0=0
while [ $j0 -lt $n ] ; do
    let j1=$j0+64
    [ "$?" = 0 ] && $bindir/bwc.pl krylov  $common interval=10 ys=$j0..$j1
    j0=$j1
done
[ "$?" = 0 ] && $bindir/bwc.pl acollect    $common -- --remove-old
[ "$?" = 0 ] && $bindir/bwc.pl lingen      $common lingen_threshold=64
[ "$?" = 0 ] && $bindir/bwc.pl :fsplit     $common splits=$all_splits
j0=0
while [ $j0 -lt $n ] ; do
    let j1=$j0+64
    [ "$?" = 0 ] && $bindir/bwc.pl mksol  $common interval=10 ys=$j0..$j1
    j0=$j1
done
[ "$?" = 0 ] && $bindir/bwc.pl gather $common interval=10
[ "$?" = 0 ] && $bindir/cleanup --ncols 64 --out $wdir/W $wdir/K.0

failure="$?"

split=${Nh}x${Nv}
b=`basename $matrix .bin`
c=$b.$split.$checksum

mdir=$wdir

cmd=`dirname $0`/convert_magma.pl

rwfile=${matrix%%bin}rw.bin
cwfile=${matrix%%bin}cw.bin
$cmd weights < $rwfile > $mdir/rw.m
$cmd weights < $cwfile > $mdir/cw.m

$cmd balancing < $wdir/$c.bin > $mdir/b.m

$cmd bmatrix < $matrix > $mdir/t.m

# For nullspace=left, the matrices are transposed.
(
echo "nullspace:=\"$nullspace\";"
echo "xtr:=func<x|$transpose_if_left(x)>;"
echo "nh:=$Nh;"
echo "nv:=$Nv;"
echo "nrp:=nv*Ceiling(nr/(nh*nv));"
echo "ncp:=nh*Ceiling(nc/(nh*nv));"
echo "Mt:=Matrix(GF(2),nh*nrp,nv*ncp,[]);"
for i in `seq 0 $((Nh-1))` ; do
    echo "nr$i:=nrp;" # nr div nh + ($i lt nr mod nh select 1 else 0);"
    echo "snr$i:=$i *nrp;" # $i*(nr div nh) + Min($i, nr mod nh);"
    for j in `seq 0 $((Nv-1))` ; do
        $cmd bmatrix < $wdir/$c.h$i.v$j> $mdir/t$i$j.m
        echo "nc$j:=ncp;" #  div nv + ($j lt nc mod nv select 1 else 0);"
        echo "snc$j:=$j*ncp;" # (nc div nv) + Min($j, nc mod nv);"
        echo "load \"$mdir/t$i$j.m\"; M$i$j:=Matrix(GF(2),Matrix(var));"
        echo "x:=RMatrixSpace(GF(2),nr$i,nc$j)!0;InsertBlock(~x,$transpose_if_left(M$i$j),1,1);M$i$j:=x;"
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

$cmd x $wdir/X > $mdir/x.m
# ok, strictly speaking the A files have nothing to do with vectors. But
# the procedure sort of works as is anyway... Same holds for F.
for f in $wdir/[ZADVCYFSWK]* ; do $cmd vector < $f > $mdir/`basename $f`.m ; done
for f in $wdir/[AF]*[0-9] ; do $cmd vector < $f > $mdir/`basename $f`.m ; done
$cmd vector < $wdir/W > $mdir/Wu.m
$cmd vector < $wdir/H1 > $mdir/H1.m

if [ "$failure" != 0 ] ; then
    echo "ERROR: partial failure encountered"
fi
