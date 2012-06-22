#! /bin/bash

build/f2/makefb $1
build/f2/sieve $1 > yield.txt
grep 'Yield' yield.txt | cut -f3 -d' ' > yield.txt.dig
sed -i -e 's/inf/1/' yield.txt.dig 
n=`wc -l yield.txt.dig | cut -f1 -d' '`
n=`expr $n - 1`
sed -e "s/\(.*\)$/\1\/$n+/" yield.txt.dig > sum.magma
sed -i -e  "s/^\/[0-9]*+/0;/" sum.magma 
magma sum.magma
quit;
rm yield.txt
rm yield.txt.dig
rm sum.magma
exit 0
