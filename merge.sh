#!/bin/sh -
cado=$HOME/CADO-src
mpidir=/share/stymphale/lix/morain/openmpi-1.2.6
hosts="painvin,parmesan,tomme,emmental,abondance,ambiorix,appenzeller"
hosts="emmental,parmesan,tomme,abondance,ambiorix,appenzeller"
np=6; mpiargs="--host $hosts --np $np"
np=11; mpiargs="-hostfile mpi.hosts --np $np"
if true
then
    echo "##### Using mpi version"
    mpirun -v --prefix $mpidir $mpiargs $cado/linalg/mpimerge $*
    root=`echo $* | awk '{print $NF}'` #### humffffff
    echo ""
    echo "Rebuilding $root from the partial files..."
    time $cado/linalg/hismerge $root `expr $np '-' 1` > $root
else
    time $cado/linalg/merge $*
fi


