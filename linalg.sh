#!/bin/sh -
linalg=linalg

# parameter from merge_linalg_sqrt.sh: $1 is $name

name=$1
mat=$name/small
ker=$name/ker_raw

if false ; then
   echo "Calling Gauss"
   time $linalg/linalg $mat 1 > $ker
else
   # use Block-Wiedemann (recommended for large matrices)
   echo "Transposing the matrix"
   time $linalg/transpose $mat $mat.tr

   echo "Calling Block-Wiedemann"
   if [ ! -e $name/bw ] ; then
      mkdir $name/bw
   fi

   time $linalg/bw/doit.pl matrix=$mat.tr mn=64 vectoring=64 multisols=1 wdir=$name/bw solution=$name/W

##   time $linalg/bw/doit.pl matrix=$mat.tr mn=64 vectoring=64 multisols=1 wdir=$name.bw solution=W

   echo "Converting dependencies to CADO format"
   # doit.pl puts the dependency file W in the directory where the matrix was
   time $linalg/bw/mkbitstrings $name/W > $ker

##   root=Examples/c80
##   time $linalg/bw/mkbitstrings $root/W > $ker

fi
