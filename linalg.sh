#!/bin/sh -
linalg=linalg

if true
then
    time $linalg/linalg $1 1 > $2
else
    tmp=`head -1 $1`
    nrows=`echo $tmp | awk '{print $1}'`
    ncols=`echo $tmp | awk '{print $2}'`
    dep=$1.dep
    args="-d $dep -p $1.pm2 -m $1 -r $nrows -c $ncols -t 128 -v1 -g -n 128"
    time $linalg/wiedemann $args
    $linalg/4cado $nrows $dep $2
fi
