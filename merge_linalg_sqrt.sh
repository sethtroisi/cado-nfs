#!/bin/sh -
# We assume that a directory $name exists
#
# Typical use: ./merge_linalg_sqrt.sh Examples/c20/c20 5 10
#
linalg=linalg
sqrt=sqrt/naive

nkermax=30; nchar=50; maxlevel=5; cwmax=10; rwmax=1000000

root=$1
if [ $# -ge 2 ]; then maxlevel=$2; fi
if [ $# -ge 3 ]; then cwmax=$3; fi
if [ $# -ge 4 ]; then rwmax=$4; fi
name=$root.$maxlevel"x"$cwmax

rels=$root.rels; nrels=`wc -l $rels | awk '{print $1}'`
poly=$root.poly
purged=$root.purged

if [ -s $purged ]
then
  echo "File $purged already exists"
else
  time $linalg/purge $nrels $rels > $purged
  if [ ! -s $purged ]; then echo "zero file $purged"; exit; fi
  excess=`head -1 $purged | awk '{nrows=$1; ncols=$2; print (nrows-ncols)}'`
  echo "excess = $excess"
  if [ $excess -lt 0 ]; then echo "exess < 0, sorry"; exit; fi
fi

echo "Performing merges"

nb_merge_max=1000000
argsa="-merge $nb_merge_max -mat $purged"
argsa="$argsa -maxlevel $maxlevel -cwmax $cwmax -rwmax $rwmax"
time $linalg/merge $argsa > $name.merge.his # 2> $name.merge.err
echo "SIZE(merge.his): `ls -s $name.merge.his`"

echo "Replaying merges"

argsr="$purged $name.merge.his $name.small $name.index"
time $linalg/replay $argsr # 2> $name.replay.err
echo "SIZE(index): `ls -s $name.index`"

echo "Performing the linear algebra phase"

time $linalg/linalg $name.small 1 > $name.ker_raw

if [ ! -s $name.ker_raw ]; then echo "Zerodim kernel, stopping"; exit; fi

nker=`wc -l $name.ker_raw | awk '{print $1}'`
if [ $nker -lt $nkermax ]; then nkermax=$nker; fi

echo "Adding characters"

args0="$purged $name.ker_raw $poly $name.index $rels"
args0="$args0 $nker $nchar"
time $linalg/characters $args0 > $name.ker

ndepmax=`wc -l $name.ker | awk '{print $1}'`
if [ $ndepmax -ge 30 ]; then ndepmax=30; fi

echo "Preparing $ndepmax squareroots"

args1="$rels $purged $name.index $name.ker $poly"
time $linalg/allsqrt $args1 0 $ndepmax ar $name.dep

echo "Entering the last phase"

mag=$name.mag
ndep=0 
while [ $ndep -lt $ndepmax ]
do
  suf=`echo $ndep | awk '{printf("%03d\n", $1)}'`
  echo "# Dependency $suf"
  rat=$name.dep.rat.$suf
  alg=$name.dep.alg.$suf
  args2="$poly $rat $alg $name.algside"
  time $sqrt/fmalgsqrt $args2 $mag
magma << EOF
load "$sqrt/finish.mag";
load "$mag";
finish(f, m, N, "$rat", "$name.algside");
quit;
EOF
  ndep=`expr $ndep '+' 1`
done
