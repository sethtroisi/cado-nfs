#!/bin/bash

#merge the strategies for each r0 and each r1 to compute the best strategies
#for each pair of cofactors (r0, r1).

#At the end, we will optain in the directory 'res_matrix' some files where
#each file is a set of optimal strategies for a pair (r0,r1).

in='res_precompt_st'
out='res_matrix'
fbb0=27
mfb0=110
fbb1=29
mfb1=110


let "r0=0"
while [ $r0 -le $mfb0 ]; do
let "r1=0"
while [ $r1 -le $mfb1 ]; do
./gst -gst -mfb0 $mfb0 -fbb0 $fbb0 -mfb1 $mfb1 -fbb1 $fbb1 -r0 $r0 -r1 $r1 -in $in -out $out &
let "r1 = r1 + 1"
done
let "r0 = r0 + 1"
done



