#!/bin/sh -
#
# Typical use: ./sqrtonly.sh Examples/c90/c90 Examples/c90/c90.20x40 0 30
#
root=$1; name=$2; ndepmin=$3; ndepmax=$4
poly=$root.poly

sqrt=sqrt/naive
mag=$name.mag
ndep=$ndepmin
while [ $ndep -lt $ndepmax ]
do
  suf=`echo $ndep | awk '{printf("%03d\n", $1)}'`
  echo "# Dependency $suf"
  rat=$name.dep.rat.$suf
  alg=$name.dep.alg.$suf
  algprod=$name.algside.$suf
  fact=$name.fact.$suf
  if [ -s $algprod ]
  then
    echo "File $algprod already exists..."
  else
    args2="$poly $rat $alg $algprod"
    time $sqrt/fmalgsqrt $args2 $mag
  fi
magma > $fact << EOF
load "$sqrt/jmc.mag";
load "$sqrt/finish.mag";
load "$mag";
finish(f, m, N, "$rat", "$alg", "$algprod");
quit;
EOF
  cat $fact
  res=`grep "PRIMEPRIME" $fact`
  if [ "X$res" != "X" ]
  then
      echo "Stopping, since we have to probable prime cofactors"
      exit
  fi
  ndep=`expr $ndep '+' 1`
done
