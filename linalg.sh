#!/bin/sh -
linalg=linalg

# parameter from merge_linalg_sqrt.sh: 
# $1 is $name
# $2 is $skip

name=$1; skip=$2
mat=$name/small
ker=$name/ker_raw

# enabling multithreading
mt=4 # default to 0???

if false ; then
   echo "Calling Gauss"
   time $linalg/linalg $mat 1 > $ker
else
   # use Block-Wiedemann (recommended for large matrices)
   echo "Transposing the matrix"
   time $linalg/transpose -in $mat -out $mat.tr -skip $skip

   echo "Calling Block-Wiedemann"
   if [ ! -e $name/bw ] ; then
      mkdir $name/bw
   fi

   time $linalg/bw/doit.pl mt=$mt matrix=$mat.tr mn=64 vectoring=64 multisols=1 wdir=$name/bw solution=$name/W

   echo "Converting dependencies to CADO format"
   # doit.pl puts the dependency file W in the directory where the matrix was
   time $linalg/bw/mkbitstrings $name/W > $ker

fi
