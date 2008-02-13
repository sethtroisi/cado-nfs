#!/bin/sh -
# We assume that a directory $name exists
#
# Typical use: ./merge_linalg_sqrt.sh params rels1 rels2 ... relsk
#
# If no relation files are given, used $root/rels.*
#
linalg=linalg
sqrt=sqrt/naive

echo "Args: $*"

params=$1; shift

# default parameters (see README.params)

#####
root=`grep "root: " $params | awk '{print $NF}'`
if [ "X$root" = "X" ]; then echo "I need at least a name..."; exit; fi

nkermax=`grep "nkermax: " $params | awk '{print $NF}'`
if [ "X$nkermax" = "X" ]; then nkermax=30; fi
nchar=`grep "nchar: " $params | awk '{print $NF}'`
if [ "X$nchar" = "X" ]; then nchar=50; fi
prune=`grep "prune: " $params | awk '{print $NF}'`
if [ "X$prune" = "X" ]; then prune=1.0; fi
maxlevel=`grep "maxlevel: " $params | awk '{print $NF}'`
if [ "X$maxlevel" = "X" ]; then maxlevel=6; fi
cwmax=`grep "cwmax: " $params | awk '{print $NF}'`
if [ "X$cwmax" = "X" ]; then cwmax=10; fi
rwmax=`grep "rwmax: " $params | awk '{print $NF}'`
if [ "X$rwmax" = "X" ]; then rwmax=100; fi
skip=`grep "skip: " $params | awk '{print $NF}'`
if [ "X$skip" = "X" ]; then skip=32; fi
##### multithreading in bw
mt=`grep "mt: " $params | awk '{print $NF}'`
if [ "X$mt" = "X" ]; then mt=0; fi
verbose=`grep "verbose: " $params | awk '{print $NF}'`
if [ "X$verbose" = "X" ]; then verbose="-v"; fi

name=`grep "name: " $params | awk '{print $NF}'`
if [ "X$name" = "X" ]; then name=$root.$prune"x"$maxlevel"x"$cwmax"x"$rwmax; fi

dir=`dirname $root`
if [ $# -eq 0 ]; then allrels=`ls $dir/rels.*`; else allrels="$*"; fi

echo $allrels

poly=$root.poly
nodup=$root.nodup
purged=$root.purged

if [ -s $nodup ]
then
  echo "File $nodup already exists"
else
  nrels=`wc -l $allrels | tail -1 | awk '{print $1}'`
  time $linalg/duplicates -nrels $nrels $allrels > $nodup
  if [ ! -s $nodup ]; then echo "zero file $nodup"; exit; fi
fi

if [ -s $purged -a $purged -nt $nodup ]
then
  echo "File $purged already exists and is newer than $nodup"
else
  nrels=`wc -l $nodup | awk '{print $1}'`
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
keep=`expr 128 '+' $skip`
argsa="-forbw -prune $prune -merge $nb_merge_max -mat $purged -keep $keep"
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

./linalg.sh $name $skip $mt

if [ ! -s $name/ker_raw ]; then echo "Zerodim kernel, stopping"; exit; fi

nker=`wc -l < $name/ker_raw`
if [ $nker -lt $nkermax ]; then nkermax=$nker; fi

echo "Adding characters"

args0="-purged $purged -ker $name/ker_raw -poly $poly -index $name/index"
args0="$args0 -rel $nodup -small $name/small"
args0="$args0 -nker $nker -nchar $nchar -skip $skip"
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
