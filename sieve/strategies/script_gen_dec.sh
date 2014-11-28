#!/bin/bash

#to precompute all decompostions of cofactors!

lim=$1
x=$2
mfb=$3
while [ $x -le $mfbr ]
do
    ./gst -gdc -lim0 $limr -mfb0 $x -out decomp_tmp/decomp_${limr}_${x}
    x=$[$x+1]
done
