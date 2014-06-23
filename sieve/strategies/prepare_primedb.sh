#!/usr/bin/env bash

# usage: argv0 <path_to_db>

PRIME_DB=$1

tmp=`mktemp`

for i in `seq 15 50`; do 
    for res in 1 5 7 11; do
        ./primedb $i $res > $tmp
        sort -R $tmp > $PRIME_DB/p$i.$res
    done
done
