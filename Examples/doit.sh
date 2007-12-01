#!/bin/sh -
# We assume that a directory $name exists
#
# Typical use: cd c20; ../doit.sh c20
#
linalg=../../linalg
sqrt=../../sqrt/naive

name=$1; nker=30; nchar=50

rels=$name.rels; nrels=`wc -l $rels | awk '{print $1}'`

$linalg/purge $rels $nrels > $name.purged
nb_merge_max=1000000
argsa="-merge $nb_merge_max -mat $name.purged -maxlevel 3 -wmax 10"
$linalg/merge $argsa > $name.merge.his
$linalg/replay $name.purged $name.merge.his $name.sparse $name.index

echo "Performing the linear algebra phase"

$linalg/linalg $name.sparse 1 > $name.ker_raw
args0="$name.purged $name.ker_raw $name.poly $name.index $name.rels"
args0="$args0 $nker $nchar"
$linalg/characters $args0 > $name.ker
args1="$name.rels $name.purged $name.index $name.ker $name.poly"
$linalg/allsqrt $args1 0 ar $name.dep

echo "Entering the last phase"

mag=$name.mag
ndep=0; ndepmax=`wc -l $name.ker | awk '{print $1}'`
while [ $ndep -lt $ndepmax ]
do
  suf=`echo $ndep | awk '{printf("%03d\n", $1)}'`
  echo "# Dependency $suf"
  rat=$name.dep.rat.$suf
  alg=$name.dep.alg.$suf
  args2="$name.poly $rat $alg $name.algside"
  $sqrt/fmalgsqrt $args2 $mag
magma << EOF
load "$sqrt/finish.mag";
load "$mag";
finish(f, m, N, "$rat", "$name.algside");
quit;
EOF
  ndep=`expr $ndep '+' 1`
done
