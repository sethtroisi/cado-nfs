#!/usr/bin/env bash                                                             
BINDIR="$1"
DEGREE="$2"
P="$3"
ADMAX="$4"
NQ="$5"
POLYSELECT="${BINDIR}/polyselect/polyselect2l"
cd ${BINDIR}
make polyselect2l
"$POLYSELECT" -N 90377629292003121684002147101760858109247336549001090677693 -degree "$DEGREE" -P "$P" -admax "$ADMAX" -nq "$NQ"
