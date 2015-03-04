#!/bin/bash
# usage: ./script_fm.sh <fbb> <lpb> <cado_bindir>

if [ $# != 3 ]; then
    echo "usage: ./script_fm.sh <fbb> <lpb> <cado_bindir>";
    exit;
fi

fbb=$1
lpb=$2
cadodir=$3

GFM=$cadodir/sieve/strategies/gfm

for m in PM1 ECM-M12 ECM-M16 ECM-B12 PP1-27 PP1-65; do 
    rm -f data_${m}_*
done

x=0
while read m b1min b1max b1step cmin cmax cstep; do
    outfile=`mktemp -p . data_${m}_XXXXXXX`
    sleep 1
    cmd="$GFM -lb $fbb -ub $lpb  -m $m -b1min $b1min -b1max $b1max -b1step $b1step -cmin $cmin -cmax $cmax -cstep $cstep -out $outfile"
    echo $cmd
    let x=$x+1
    if [ $x == 10 ]; then
        x=0
        $cmd
    else
        $cmd &
    fi
done << EOF
PM1      50   300   20  10  50   5 
PM1      300  600   20  10  50   5 
PM1      50   300   20  50  200  10 
PM1      300  600   20  50  200  10 
PM1      600  1000  20  10  50   5 
PM1      600  1000  20  50  200  10
ECM-M12  50   300   20  10  50   5 
ECM-M12  300  600   20  10  50   5 
ECM-M12  50   300   20  50  300  10 
ECM-M12  300  600   20  50  300  10 
ECM-M12  600  1500  20  10  50   5 
ECM-M12  600  1500  20  50  300  10
ECM-M16  50   300   20  10  50   5 
ECM-M16  300  600   20  10  50   5 
ECM-M16  50   300   20  50  300  10 
ECM-M16  300  600   20  50  300  10 
ECM-M16  600  1500  20  10  50   5 
ECM-M16  600  1500  20  50  300  10 
ECM-B12  50   300   20  10  50   5 
ECM-B12  300  600   20  10  50   5 
ECM-B12  50   300   20  50  300  10 
ECM-B12  300  600   20  50  300  10 
ECM-B12  600  1500  20  10  50   5 
ECM-B12  600  1500  20  50  300  10 
PP1-27   50   300   20  10  50   5 
PP1-27   300  600   20  10  50   5 
PP1-27   50   300   20  50  200  10 
PP1-27   300  600   20  50  200  10 
PP1-27   600  1500  20  10  50   5 
PP1-27   600  1500  20  50  200  10 
PP1-65   50   300   20  10  50   5 
PP1-65   300  600   20  10  50   5 
PP1-65   50   300   20  50  200  10 
PP1-65   300  600   20  50  200  10 
PP1-65   600  1500  20  10  50   5 
PP1-65   600  1500  20  50  200  10
EOF

wait

cat data_PM1_* > data_PM1
cat data_PP1-27_* > data_PP1_27
cat data_PP1-65_* > data_PP1_65
cat data_ECM-M12_* > data_ECM_M12
cat data_ECM-M16_* > data_ECM_M16
cat data_ECM-B12_* > data_ECM_B12

$GFM -fch -fch_in data_PP1_27 -fch_out data_PP1_27_ch
$GFM -fch -fch_in data_PM1 -fch_out data_PM1_ch
$GFM -fch -fch_in data_PP1_65 -fch_out data_PP1_65_ch
$GFM -fch -fch_in data_ECM_M12 -fch_out data_ECM_M12_ch
$GFM -fch -fch_in data_ECM_M16 -fch_out data_ECM_M16_ch
$GFM -fch -fch_in data_ECM_B12 -fch_out data_ECM_B12_ch

# for m in PM1 ECM-M12 ECM-M16 ECM-B12 PP1-27 PP1-65; do 
#     rm -f data_${m}_*
# done
