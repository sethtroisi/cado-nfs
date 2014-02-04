#!/usr/bin/env bash

# Input files. First parameter is the polynomial file to use, second is root of the build directory, third is the expected SHA1 value
POLY="$1"
BINDIR="$2"
REFERENCE_SHA1="$3"

if [[ -z "${BINDIR}" || ! -d "${BINDIR}" ]]
then
  echo "Build directory ${BINDIR} is not a directory" >&2
  exit 1
fi

if [ -z "${REFERENCE_SHA1}" ]
then
  echo "No reference SHA1 value specified" >&2
  exit 1
fi

if [[ -z "$rlim" || -z "$alim" || -z "$lpbr" || -z "$lpba" || -z "$maxbits" || -z "$mfbr" || -z "$mfba" || -z "$rlambda" || -z "$alambda" || -z "$I" || -z "$q0" || -z "$q1" ]]
then
  echo "Required shell environment variable not set" >&2
  exit 1
fi

MAKEFB="${BINDIR}/sieve/makefb"
cd ${BINDIR}
make makefb
LAS="${BINDIR}/sieve/las"
BASENAME="`basename "${POLY}"`"
BASENAME="${BASENAME%%.*}"

# Output files
TMPDIR=`mktemp -d /tmp/cadotest.XXXXXXXXXX`
FB="${TMPDIR}/${BASENAME}.roots"
RELS="${TMPDIR}/${BASENAME}.rels"

"$MAKEFB" -poly "$POLY" -alim $alim -maxbits $maxbits -out "${FB}" || exit 1
"$LAS" -poly "$POLY" -fb "${FB}" -I "$I" -rlim "$rlim" -lpbr "$lpbr" -mfbr "$mfbr" -rlambda "$rlambda" -alim "$alim" -lpba "$lpba" -mfba "$mfba" -alambda "$alambda" -q0 "$q0" -q1 "$q1" -out "${RELS}" || exit 1

SHA1=`grep -v "^#" "${RELS}" | sort | sha1sum`
SHA1="${SHA1%% *}"
if [ "${SHA1}" != "${REFERENCE_SHA1}" ]
then
  echo "$0: Got SHA1 of ${SHA1} but expected ${REFERENCE_SHA1}. Files remain in ${TMPDIR}"
  exit 1
fi


rm -f "${FB}" "${RELS}"
rmdir "${TMPDIR}"
