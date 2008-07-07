#!/bin/sh -
cado=${cado-.}
linalg=$cado/linalg

# parameter from merge_linalg_sqrt.sh: 
# $1 is $mat
# $2 is $skip
# $3 is $outdir
# [optional] $4 is $mt

mat=$1; skip=$2; outdir=$3
if [ $# -ge 4 ]; then mt=$4; else mt=0; fi

ker=$outdir/ker_raw

if false ; then
   echo "Calling Gauss"
   time $linalg/linalg $mat 1 > $ker
else
   # use Block-Wiedemann (recommended for large matrices)
   echo "Transposing the matrix"
   if [ -s $mat.tr -a $mat.tr -nt $mat ]
   then
       echo "File $mat.tr exists and is newer than $mat";
   else
       time $linalg/transpose -T $outdir -in $mat -out $mat.tr -skip $skip
   fi

   echo "Calling Block-Wiedemann"
   wdir=$outdir/bw
   solution=$outdir/W

   if [ ! -e $wdir ] ; then
      mkdir $wdir
   fi

   time $linalg/bw/bw.pl mt=$mt matrix=$mat.tr mn=64 vectoring=64 multisols=1 wdir=$wdir solution=$solution remove_input=1
# to use block Lanczos instead of block Wiedemann, just comment the above line
# and uncomment the following one
#   time $linalg/bl/bl.pl matrix=$mat.tr wdir=$wdir solution=$solution

   echo "Converting dependencies to CADO format"
   # bw.pl puts the dependency file W in the directory where the matrix was
   time $linalg/bw/mkbitstrings $solution > $ker

fi
