#!/usr/bin/env bash

FB="$1"
LAS="$2"
SRCDIR="$3"
TODO="$4"

cmd="${LAS} -poly ${SRCDIR}/params/F9.poly -lim0 2300000 -lim1 1200000 -lpb0 26 -lpb1 26 -mfb0 52 -mfb1 52 -I 12 -fb1 ${FB} -todo ${TODO}"
echo $cmd
## There is a bug on 32-bit archs so that we find only 418 relations
## instead of 420. This is referenced in the tracker as Bug #18751
$cmd | grep -e "Total 420 reports" -e "Total 418 reports" || exit 1

exit 0
