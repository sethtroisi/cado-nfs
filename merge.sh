#!/bin/sh -
#
# tricky call: $0 -out mergehis ....
#

usempi=false
resume=false

cado=$HOME/CADO-src
mpidir=/share/stymphale/lix/morain/openmpi-1.2.6

hosts="painvin,parmesan,tomme,emmental,abondance,ambiorix,appenzeller"
hosts="parmesan,emmental,tomme,abondance,ambiorix,appenzeller"
np=4; mpiargs="--host $hosts --np $np"

np=4; mpiargs="-hostfile mpi.hosts --np $np"

mpicmd="mpirun -v --prefix $mpidir $mpiargs"

root=$2
if $usempi
then
  if $resume
  then
    echo "##### Merge with 2"
    # last arguments override their first values...
    $cado/linalg/merge $* -out $root.tmp -maxlevel 2 -cwmax 10
    echo "##### Using mpi version on `date`"
    $mpicmd $cado/linalg/mpimerge $* -resume $root.tmp
  else
    $mpicmd $cado/linalg/mpimerge $*
  fi
  echo ""
  echo "Rebuilding $root from the partial files on `date`"
  time $cado/linalg/hismerge $root `expr $np '-' 1` > $root
else
  if $resume 
  then
    echo "##### Trying to stop and resume"
    echo "##### Merge with 2"
    # last arguments override their first values...
    $cado/linalg/merge $* -out $root.tmp -maxlevel 2
    echo "##### Resuming and going to original maxlevel"
    $cado/linalg/merge $* -resume $root.tmp
  else
    time $cado/linalg/merge $*
  fi
fi
