#!/bin/sh

w=$1
file=$2

awk '/^[0-9]+  : mul took [0-9]+ ns$/ { print $5, $1; }' < $file | \
	sort -n | head -1 | \
	(read t k ; 	\
	echo "Choosing k=$k ($t ns)" ;	\
	make -s mul-inlines-k$k.c w=$w ; mv mul-inlines-k$k.c mul-inlines.c)
