#!/bin/bash

#precompute for each size of cofactors and each side (RAT or ALG), the convex hull of
#strategies!


#the sieving region
#remark: when I will embeded the code to compute the decomposition files, 
#I will change fbb0 and fbb1 by lim0 and lim1. 
fbb0=27
lpb0=37
mfb0=110

fbb1=29
lpb1=37
mfb1=110

#directory which store our decomposition files!
decomp='/localdisk/trichard/results/decomp_cofactor/decomp_tmp/'
#file which contain our factoring methods!
methods='Echantillon_21_22'
#directory where the output file will be stored!
out='res_precompt_st/'


let "r0=0"
while [ $r0 -le $mfb0 ]; do
pathname_decomp=$decomp"/decomp_"$fbb0"_"$r0
./gst -gst_r -lpb0 $lpb0 -fbb0 $fbb0 -r $r0 -in $methods -decomp $pathname_decomp -out $out &
let "r0 = r0 + 1"
done

let "r1=0"
while [ $r1 -le $mfb1 ]; do
pathname_decomp=$decomp"/decomp_"$fbb1"_"$r1
./gst -gst_r -lpb0 $lpb1 -fbb0 $fbb1 -r $r1 --in $methods -decomp $pathname_decomp -out $out &
let "r1 = r1 + 1"
done
