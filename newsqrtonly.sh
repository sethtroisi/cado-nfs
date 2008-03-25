#!/bin/sh -
#
# Typical use: ./newsqrtonly.sh Examples/c90/c90 Examples/c90/c90.20x40 0 30
#
root=$1; name=$2; ndepmin=$3; ndepmax=$4
poly=$root.poly

cado=${cado-.}
sqrt=$cado/sqrt/naive
ndep=$ndepmin
while [ $ndep -lt $ndepmax ]
do
  suf=`echo $ndep | awk '{printf("%03d\n", $1)}'`
  echo "# Dependency $suf"
  rat=$name/dep.rat.$suf
  alg=$name/dep.alg.$suf
  fact=$name/fact.$suf
  $sqrt/algsqrt $alg $rat $poly > $fact
  cat $fact
  res=`grep "Failed" $fact`
  if [ "X$res" == "X" ]
  then
      echo "Stopping, since we have probably finished"
      exit
  fi
  ndep=`expr $ndep '+' 1`
done
