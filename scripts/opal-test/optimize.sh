#!/bin/sh
# this script automatically optimizes sieving parameters
# Usage: optimize.sh params.cxx cxx.polyselect2.poly
# Puts the optimized file in params.cxx.opt in the current directory.
# Remark: the 'nomad' binary must be in $PATH (see README)
# The CADO_BUILD environment variable must contain the CADO-NFS build
# directory (makefb and las are taken from $CADO_BUILD/sieve)

cwd=`pwd`
params=$1
poly=`basename $2`
d=`mktemp -d`
echo "Working directory:" $d
cp $2 las_optimize.py report.py $d
sed "s/c59.polyselect2.poly/$poly/g" las_run.py > $d/las_run.py
rlim=`grep rlim $params | cut -d= -f2`
alim=`grep alim $params | cut -d= -f2`
lpbr=`grep lpbr $params | cut -d= -f2`
lpba=`grep lpba $params | cut -d= -f2`
mfbr=`grep mfbr $params | cut -d= -f2`
mfba=`grep mfba $params | cut -d= -f2`
grep ncurves0 $params > /dev/null
if [ $? -eq 0 ]; then
   ncurves0=`grep ncurves0 $params | cut -d= -f2`
else
   ncurves0=10
fi
grep ncurves1 $params > /dev/null
if [ $? -eq 0 ]; then
   ncurves1=`grep ncurves1 $params | cut -d= -f2`
else
   ncurves1=10
fi
I=`grep I $params | cut -d= -f2`
rlim_min=`expr $rlim / 2`
rlim_max=`expr $rlim \* 2`
alim_min=`expr $alim / 2`
alim_max=`expr $alim \* 2`
lpbr_min=`expr $lpbr - 1`
lpbr_max=`expr $lpbr + 1`
lpba_min=`expr $lpba - 1`
lpba_max=`expr $lpba + 1`
mfbr_min=$lpbr_min
mfbr_max=`expr $lpbr_max \* 3`
mfba_min=$lpba_min
mfba_max=`expr $lpba_max \* 3`
if [ $ncurves0 -gt 3 ]; then
ncurves0_min=`expr $ncurves0 - 3`
else
ncurves0_min=0
fi
ncurves0_max=`expr $ncurves0 + 3`
if [ $ncurves1 -gt 3 ]; then
ncurves1_min=`expr $ncurves1 - 3`
else
ncurves1_min=0
fi
ncurves1_max=`expr $ncurves1 + 3`
I_min=`expr $I - 1`
I_max=`expr $I + 1`
sed "s/rlim_def/$rlim/g" las_decl_template.py | \
sed "s/rlim_min/$rlim_min/g" | sed "s/rlim_max/$rlim_max/g" | \
sed "s/alim_def/$alim/g" | sed "s/alim_min/$alim_min/g" | \
sed "s/alim_max/$alim_max/g" | \
sed "s/lpbr_def/$lpbr/g" | sed "s/lpbr_min/$lpbr_min/g" | \
sed "s/lpbr_max/$lpbr_max/g" | \
sed "s/lpba_def/$lpba/g" | sed "s/lpba_min/$lpba_min/g" | \
sed "s/lpba_max/$lpba_max/g" | \
sed "s/mfbr_def/$mfbr/g" | sed "s/mfbr_min/$mfbr_min/g" | \
sed "s/mfbr_max/$mfbr_max/g" | \
sed "s/mfba_def/$mfba/g" | sed "s/mfba_min/$mfba_min/g" | \
sed "s/mfba_max/$mfba_max/g" | \
sed "s/ncurves0_def/$ncurves0/g" | sed "s/ncurves0_min/$ncurves0_min/g" | \
sed "s/ncurves0_max/$ncurves0_max/g" | \
sed "s/ncurves1_def/$ncurves1/g" | sed "s/ncurves1_min/$ncurves1_min/g" | \
sed "s/ncurves1_max/$ncurves1_max/g" | \
sed "s/I_def/$I/g" | sed "s/I_min/$I_min/g" | sed "s/I_max/$I_max/g" \
> $d/las_declaration.py
cd $d
python las_optimize.py
# optimized parameters are in nomad-solution.nnn.txt
f=`ls -t nomad-solution.*.txt | head -1`
rlim_opt=`head -1 $f`
alim_opt=`head -2 $f | tail -1`
lpbr_opt=`head -3 $f | tail -1`
lpba_opt=`head -4 $f | tail -1`
mfbr_opt=`head -5 $f | tail -1`
mfba_opt=`head -6 $f | tail -1`
ncurves0_opt=`head -7 $f | tail -1`
ncurves1_opt=`head -8 $f | tail -1`
I_opt=`head -9 $f | tail -1`
echo "Optimal parameters:"
echo "rlim=" $rlim_opt
echo "alim=" $alim_opt
echo "lpbr=" $lpbr_opt
echo "lpba=" $lpba_opt
echo "mfbr=" $mfbr_opt
echo "mfba=" $mfba_opt
echo "ncurves0=" $ncurves0_opt
echo "ncurves1=" $ncurves1_opt
echo "I=" $I_opt
cd $cwd
sed "s/rlim.*=.*$/rlim = $rlim_opt/g" $params | \
sed "s/alim.*=.*$/alim = $alim_opt/g" | \
sed "s/lpbr.*=.*$/lpbr = $lpbr_opt/g" | \
sed "s/lpba.*=.*$/lpba = $lpba_opt/g" | \
sed "s/mfbr.*=.*$/mfbr = $mfbr_opt/g" | \
sed "s/mfba.*=.*$/mfba = $mfba_opt/g" | \
sed "s/ncurves0.*=.*$/ncurves0 = $ncurves0_opt/g" | \
sed "s/ncurves1.*=.*$/ncurves1 = $ncurves1_opt/g" | \
sed "s/I.*=.*$/I = $I_opt/g" > $params.opt
/bin/rm -fr $d
