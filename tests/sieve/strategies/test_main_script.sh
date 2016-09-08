#!/bin/bash

CADO_NFS_SOURCE_DIR=$1
CADO_NFS_BINARY_DIR=$2

t=`mktemp -d ${TMPDIR-/tmp}/cado-check.XXXXXXX`
cd $t

cat > c120.poly <<EOF
n: 227010481295437363334259960947493668895875336466084780038173258247009162675779735389791151574049166747880487470296548479
skew: 33200.0
c0: -773897878252443925372185840
c1: -594857401088853151809341
c2: -29938251437178967092
c3: 555583138207303
c4: -3025583902
c5: 99960
Y0: -85694025999552843022178
Y1: 54675511261073897
# MurphyE = 3.06e-09
# lognorm 36.62
# f(x) = 99960*x^5-3025583902*x^4+555583138207303*x^3-29938251437178967092*x^2-594857401088853151809341*x-773897878252443925372185840
# g(x) = 54675511261073897*x-85694025999552843022178
EOF

lim0=2371302
lim1=3837507
lpb0=26
lpb1=26
mfb0=54
mfb1=53
I=12

export CREATE_STRAT_FILE_TEST="yes, please"
$CADO_NFS_SOURCE_DIR/sieve/strategies/create_strat_file.sh $lim0 $lpb0 $mfb0 $lim1 $lpb1 $mfb1 c120.poly $I $CADO_NFS_BINARY_DIR

nl=`wc -l final_st | cut -d " " -f 1`
if [ $nl == "382" ]; then
    rm -rf $t
    exit 0
else
    exit 1
fi
