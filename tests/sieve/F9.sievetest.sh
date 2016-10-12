#!/usr/bin/env bash

# This file defines the sieving parameter and the reference SHA1 value, then calls sievetest.sh

FB="$1"
LAS="$2"
SRCDIR="$3"
CHECKSUM_FILE="$4"
SOURCE_TEST_DIR="`dirname "$0"`"
shift 4

REFERENCE_SHA1="c12c10808914fb445dd32e6801ffb8f320db36e6"
# The Git revision that created the REFERENCE_SHA1 hash
REFERENCE_REVISION="0bd0c87d93c709279cc5c3dd471292b28ae04198" # 430 relations
# Previous revisions and number of relations found (most recent first):
# 9253f54766c5438fe3198053f9edd811ee254a17 (421 relations)

lim0=2300000
lim1=1200000
lpb0=26
lpb1=26
maxbits=10
mfb0=52
mfb1=52
lambda0=2.1
lambda1=2.2
I=12
q0=1200000
q1=1200200

export lim0 lim1 lpb0 lpb1 maxbits mfb0 mfb1 lambda0 lambda1 I q0 q1
"${SOURCE_TEST_DIR}"/sievetest.sh "${FB}" "${LAS}" "${SRCDIR}/parameters/polynomials/F9.poly" "${REFERENCE_SHA1}" "${REFERENCE_REVISION}" "${CHECKSUM_FILE}" -v -v "$@" || exit 1

exit 0
