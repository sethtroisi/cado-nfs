#!/bin/bash

#precompute for each size of cofactors and each side (RAT or ALG), the
#convex hull of strategies!


#the sieving region
lim0=$1 # say 2^27
lpb0=$2 # say 37
mfb0=$3 # say 110

#directory which store our decomposition files!
decomp=$5 # say '/localdisk/trichard/results/decomp_cofactor/decomp_tmp/'
#file which contains our factoring methods (PM1, PP1, EC)!
methods=$6 # say 'All_methods'
#directory where the output file will be stored, because this binary
#generates a file for each r.
out=$7 # say 'res_precompt_st/'

let "r0=0"
while [ $r0 -le $mfb0 ]; do
pathname_decomp=$decomp"/decomp_"$lim0"_"$r0
./gst -gst_r -lim0 $fbb0 -lpb0 $lpb0 -r0 $r0 -in $methods -decomp $pathname_decomp -out $out &
let "r0 = r0 + 1"
done

