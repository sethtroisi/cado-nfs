#!/bin/sh

d="$1"
shift

p="$1"
shift

for k in "$@" ; do
    echo $k " : "`"$d/$p$k" -nokara 1` | tee /dev/stderr
done | \
awk '/^[0-9]+  : mul took [0-9]+ ns$/ { print $5, $1; }' | \
sort -n | head -1 | \
(read t k ; 	\
echo "Choosing k=$k ($t ns)" ;	\
cp "$d/mul-inlines-k$k.c" "$d/mul-inlines.c" ; \
touch "$d/mul-inlines.c")

