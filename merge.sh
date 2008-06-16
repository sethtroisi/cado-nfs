#!/bin/sh -
cado=$HOME/CADO-src
mpidir=/share/stymphale/lix/morain/openmpi-1.2.6
if false
then
    echo "##### Using mpi version"
    np=3; hosts="painvin,parmesan,tomme"
    np=4; hosts="painvin,parmesan,tomme,emmental"
    mpirun -v --prefix $mpidir --np $np --host $hosts \
	$cado/linalg/mpimerge $*
    root=`echo $* | awk '{print $NF}'` #### humffffff
    echo ""
    echo "Rebuilding $root from the partial files..."
    $cado/linalg/hismerge $root `expr $np '-' 1` > $root
    echo "... done"
    echo ""
else
    $cado/linalg/merge $*
fi


