#!/usr/bin/env bash

# This file defines the sieving parameter and the reference SHA1 value, then calls sievetest.sh

SRCDIR="$1"
BINDIR="$2"
SOURCE_TEST_DIR="`dirname "$0"`"

REFERENCE_SHA1="a5bb1736493c387e415da14d4470761c2c29f463"
# The Git revision that created the REFERENCE_SHA1 hash
REFERENCE_REVISION="9253f54766c5438fe3198053f9edd811ee254a17"

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
"${SOURCE_TEST_DIR}"/sievetest.sh "${SRCDIR}/params/F9.poly" "${BINDIR}" "${REFERENCE_SHA1}" "${REFERENCE_REVISION}" || exit 1

exit 0
