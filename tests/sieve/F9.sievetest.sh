#!/usr/bin/env bash

# This file defines the sieving parameter and the reference SHA1 value, then calls sievetest.sh

SRCDIR="$1"
BINDIR="$2"
SOURCE_TEST_DIR="`dirname "$0"`"

REFERENCE_SHA1="5ecca4127483d68d1b3cab3c62c8f6ff1ab371af"

rlim=2300000
alim=1200000
lpbr=26
lpba=26
maxbits=10
mfbr=52
mfba=52
rlambda=2.1
alambda=2.2
I=12
q0=1200000
q1=1200200

export rlim alim lpbr lpba maxbits mfbr mfba rlambda alambda I q0 q1
"${SOURCE_TEST_DIR}"/sievetest.sh "${SRCDIR}/params/F9.poly" "${BINDIR}" "${REFERENCE_SHA1}" || exit 1

exit 0
