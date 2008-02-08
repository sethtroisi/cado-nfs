#!/bin/sh -
# We assume that a directory $name exists
#
# Typical use: ./merge_linalg_sqrt.sh Examples/c20/c20 5 10 100 [v]
#
linalg=linalg
sqrt=sqrt/naive

# default parameters (see README.params)
nkermax=30; nchar=50; prune=1.0; maxlevel=6; cwmax=10; rwmax=100

root=$1
if [ $# -ge 2 ]; then prune=$2; fi
if [ $# -ge 3 ]; then maxlevel=$3; fi
if [ $# -ge 4 ]; then cwmax=$4; fi
if [ $# -ge 5 ]; then rwmax=$5; fi
name=$root.$prune"x"$maxlevel"x"$cwmax"x"$rwmax
if [ $# -ge 6 ]; then verbose="-v"; fi

echo "Args: $*"

rels=$root.rels
poly=$root.poly
nodup=$root.nodup
purged=$root.purged

if [ -s $nodup -a $nodup -nt $rels ]
then
  echo "File $nodup already exists and is newer than $rels"
else
  nrels=`wc -l $rels | awk '{print $1}'`
  time $linalg/duplicates -nrels $nrels $rels > $nodup
  if [ ! -s $nodup ]; then echo "zero file $nodup"; exit; fi
fi

# update
nrels=`wc -l $nodup | awk '{print $1}'`

if [ -s $purged -a $purged -nt $nodup ]
then
  echo "File $purged already exists and is newer than $nodup"
else
  time $linalg/purge -poly $poly -nrels $nrels $nodup > $purged
  if [ ! -s $purged ]; then echo "zero file $purged"; exit; fi
  excess=`head -1 $purged | awk '{nrows=$1; ncols=$2; print (nrows-ncols)}'`
  echo "excess = $excess"
  if [ $excess -le 0 ]
  then 
      echo "excess <= 0, sorry"
      /bin/rm -f $purged
      exit
  fi
fi

echo "Performing merges"

if [ -d $name ]
then
    echo "Directory $name already exists"
else
    echo "Creating directory $name"
    mkdir $name
fi

nb_merge_max=1000000
argsa="-forbw -prune $prune -merge $nb_merge_max -mat $purged"
argsa="$argsa -maxlevel $maxlevel -cwmax $cwmax -rwmax $rwmax $verbose"
time $linalg/merge $argsa > $name/merge.his # 2> $name.merge.err
echo "SIZE(merge.his): `ls -s $name/merge.his`"

echo "Replaying merges"

bwcostmin=`tail $name/merge.his | grep "BWCOSTMIN:" | awk '{print $NF}'`
argsr="$purged $name/merge.his $name/small $name/index"
time $linalg/replay $argsr $bwcostmin # 2> $name.replay.err

if [ ! -s $name/index ]
then
    echo "Index file $name/index does not exist, stopping"
    exit
fi
echo "SIZE(index): `ls -s $name/index`"

if [ ! -s $name/small ]
then
    echo "Small matrix $name/small does not exist, stopping"
    exit
fi

echo "Performing the linear algebra phase"

./linalg.sh $name

if [ ! -s $name/ker_raw ]; then echo "Zerodim kernel, stopping"; exit; fi

nker=`wc -l < $name/ker_raw`
if [ $nker -lt $nkermax ]; then nkermax=$nker; fi

echo "Adding characters"

args0="$purged $name/ker_raw $poly $name/index $nodup"
args0="$args0 $nker $nchar"
time $linalg/characters $args0 > $name/ker

ndepmax=`wc -l $name/ker | awk '{print $1}'`
if [ $ndepmax -ge 30 ]; then ndepmax=30; fi

echo "Preparing $ndepmax squareroots"

args1="$nodup $purged $name/index $name/ker $poly"
time $linalg/allsqrt $args1 0 $ndepmax ar $name/dep

echo "Entering the last phase"

./newsqrtonly.sh $root $name 0 $ndepmax
## If this fails, use the backup version, based on magma:
##   ./sqrtonly.sh $root $name 0 $ndepmax
