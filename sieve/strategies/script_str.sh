#!/bin/bash

#merge the strategies for each r0 and each r1 to compute the best
#strategies for each pair of cofactors (r0, r1).

#At the end, we will optain in the directory 'res_matrix' some files
#where each file is a set of optimal strategies for a pair (r0,r1).

fbb0=$1 # say 27
mfb0=$2 # say 110
fbb1=$3 # say 29
mfb1=$4 # say 110

#the directory where the strategies related on each value of r were
#stored.
in=$5 # say 'res_precompt_st'
#the directory where the files results will be stored. One file will
#be generated for each pair of bit size (r0,r1).
out=$6 # say 'res_matrix'

let "r0=0"
while [ $r0 -le $mfb0 ]; do
let "r1=0"
while [ $r1 -le $mfb1 ]; do
./gst -gst -mfb0 $mfb0 -fbb0 $fbb0 -mfb1 $mfb1 -fbb1 $fbb1 -r0 $r0 -r1 $r1 -in $in -out $out &
let "r1 = r1 + 1"
done
let "r0 = r0 + 1"
done



