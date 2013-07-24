#!/bin/bash

# usage: argv0 <method> <B1> <B2> <sigma>

method=$1
B1=$2
B2=$3
sigma=$4

## Path where Cado-nfs has been compiled:
CADO="/tmp/cado-nfs-build-master"
BENCHER="$CADO/sieve/ecm/testbench"

## Path where list of primes is to be found
PRIME_DB="/localdisk/primes"

## A prime of size 64, 126, 192 for measuring runtime.
p1w=18446744073709551557
p2w=85070591730234615865843651857942052727
p3w=1569275433846670190958947355801916604025588861116008628213

if [ "$method" = "ecmm12" ]; then
   extra="-ep"
else
   extra=""
fi

## Print header
METHOD=`echo $method | tr '[:lower:]' '[:upper:]'`
echo "METHOD=$METHOD B1=$B1 B2=$B2 sigma=$sigma $extra"

## Bench
echo "  Bench:"
for i in 1 2 3; do
    if [ $i -eq "1" ]; then
        p=$p1w
    elif [ $i -eq "2" ]; then
        p=$p2w
    elif [ $i -eq "3" ]; then
        p=$p3w
    else
        echo "error..."
        exit 1
    fi
    tm=`yes $p | head -1000 | $BENCHER -inp /dev/stdin -$method $B1 $B2 $sigma $extra | awk '/Total/ {print $7/1000.0}'`
    echo "    $i words: $tm"
done


# Evaluate probabilities of finding a prime of each size, for each class
# modulo 12.
#
# This is done by trying at most 2000 primes of the given size and
# congruence, taken from the PRIME_DB directory.
#   Sizes belongs to [15..50] bits
#   Congruences belongs to {1, 5, 7, 11} mod 12

echo "  Success Probability for p of given bit size for 1,5,7,11 mod 12"
for i in `seq 15 50`; do
    echo -n "    $i:"
    for res in 1 5 7 11; do 
        ratio=`head -2000 $PRIME_DB/p${i}.$res | $BENCHER -inp /dev/stdin -cof 1152921504606846883 -$method $B1 $B2 $sigma | awk '/Ratio/ {print $2}'`
        printf " %.2f" $ratio
    done
    echo
done
