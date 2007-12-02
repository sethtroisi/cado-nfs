#!/bin/sh -
# We assume that a directory $name exists
#
# Typical use: cd c20; ../doit.sh c20 5
#
linalg=../../linalg
sqrt=../../sqrt/naive

nker=30; nchar=50; maxlevel=5; cwmax=10

name=$1
if [ $# -ge 2 ]; then maxlevel=$2; fi
if [ $# -ge 3 ]; then cwmax=$3; fi

rels=$name.rels; nrels=`wc -l $rels | awk '{print $1}'`

$linalg/purge $rels $nrels > $name.purged

echo "Performing merges"

nb_merge_max=1000000
argsa="-merge $nb_merge_max -mat $name.purged"
argsa="$argsa -maxlevel $maxlevel -cwmax $cwmax"
time $linalg/merge $argsa > $name.merge.his # 2> $name.merge.err

echo "Replaying merges"

time $linalg/replay $name.purged $name.merge.his $name.small $name.index

echo "Performing the linear algebra phase"

$linalg/linalg $name.small 1 > $name.ker_raw

echo "Adding characters"

args0="$name.purged $name.ker_raw $name.poly $name.index $name.rels"
args0="$args0 $nker $nchar"
time $linalg/characters $args0 > $name.ker

echo "Preparing squareroots"

args1="$name.rels $name.purged $name.index $name.ker $name.poly"
time $linalg/allsqrt $args1 0 ar $name.dep

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
  time $sqrt/fmalgsqrt $args2 $mag
magma << EOF
load "$sqrt/finish.mag";
load "$mag";
finish(f, m, N, "$rat", "$name.algside");
quit;
EOF
  ndep=`expr $ndep '+' 1`
done
