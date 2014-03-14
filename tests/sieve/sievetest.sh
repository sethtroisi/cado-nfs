#!/usr/bin/env bash

# Input files. First parameter is the polynomial file to use, second is root of the build directory, third is the expected SHA1 value
MAKEFB="$1"
LAS="$2"
POLY="$3"
REFERENCE_SHA1="$4"
REFERENCE_REVISION="$5"
shift 5

if [[ -z "${MAKEFB}" || ! -x "${MAKEFB}" ]]
then
  echo "Makefb binary ${MAKEFB} is not a executable file" >&2
  exit 1
fi

if [[ -z "${LAS}" || ! -x "${LAS}" ]]
then
  echo "Las binary ${LAS} is not a executable file" >&2
  exit 1
fi

if [ -z "${REFERENCE_SHA1}" ]
then
  echo "No reference SHA1 value specified" >&2
  exit 1
fi


if [[ -z "$rlim" || -z "$alim" || -z "$lpbr" || -z "$lpba" || -z "$maxbits" \
      || -z "$mfbr" || -z "$mfba" || -z "$rlambda" || -z "$alambda" || -z "$I" \
      || -z "$q0" || ( -z "$q1" && -z "$rho" ) ]]
then
  echo "Required shell environment variable not set" >&2
  exit 1
fi

BASENAME="`basename "${POLY}"`"
BASENAME="${BASENAME%%.*}"

# Output files
TMPDIR=`mktemp -d /tmp/cadotest.XXXXXXXXXX`
# Make temp direcotry world-readable for easier debugging
chmod a+rx "${TMPDIR}"
FB="${TMPDIR}/${BASENAME}.roots"
RELS="${TMPDIR}/${BASENAME}.rels"
FBC="${TMPDIR}/${BASENAME}.fbc"

if [ -n "$q1" ]
then
  end=("-q1" "$q1")
fi
if [ -n "$rho" ]
then
  end=("-rho" "$rho")
fi
echo "$MAKEFB" -poly "$POLY" -alim $alim -maxbits $maxbits -out "${FB}"
"$MAKEFB" -poly "$POLY" -alim $alim -maxbits $maxbits -out "${FB}" || exit 1
# first exercise the -fbc command-line option to create a cache file
echo "$LAS" -poly "$POLY" -fb "${FB}" -I "$I" -rlim "$rlim" -lpbr "$lpbr" -mfbr "$mfbr" -rlambda "$rlambda" -alim "$alim" -lpba "$lpba" -mfba "$mfba" -alambda "$alambda" -q0 "$q0" -q1 "$q0" -out "${RELS}" -fbc "${FBC}" "$@"
"$LAS" -poly "$POLY" -fb "${FB}" -I "$I" -rlim "$rlim" -lpbr "$lpbr" -mfbr "$mfbr" -rlambda "$rlambda" -alim "$alim" -lpba "$lpba" -mfba "$mfba" -alambda "$alambda" -q0 "$q0" -q1 "$q0" -out "${RELS}" -fbc "${FBC}" "$@" || exit 1
# then use the cache file created above
echo "$LAS" -poly "$POLY" -fb "${FB}" -I "$I" -rlim "$rlim" -lpbr "$lpbr" -mfbr "$mfbr" -rlambda "$rlambda" -alim "$alim" -lpba "$lpba" -mfba "$mfba" -alambda "$alambda" -q0 "$q0" "${end[@]}" -out "${RELS}" -fbc "${FBC}" "$@"
"$LAS" -poly "$POLY" -fb "${FB}" -I "$I" -rlim "$rlim" -lpbr "$lpbr" -mfbr "$mfbr" -rlambda "$rlambda" -alim "$alim" -lpba "$lpba" -mfba "$mfba" -alambda "$alambda" -q0 "$q0" "${end[@]}" -out "${RELS}" -fbc "${FBC}" "$@" || exit 1


SHA1BIN=sha1sum
if ! echo | "$SHA1BIN" > /dev/null 2>&1
then
  SHA1BIN=sha1
fi

# Try to make sort produce some well-defined ordering on the integers
export LC_ALL=C
export LANG=C
export LANGUAGE=C

SHA1=`grep -v "^#" "${RELS}" | sort -n | ${SHA1BIN}` || exit 1
SHA1="${SHA1%% *}"
if [ "${SHA1}" != "${REFERENCE_SHA1}" ]
then
  if [ -n "${REFERENCE_REVISION}" ]
  then
    REFMSG=", as created by Git revision ${REFERENCE_REVISION}"
  fi
  echo "$0: Got SHA1 of ${SHA1} but expected ${REFERENCE_SHA1}${REFMSG}. Files remain in ${TMPDIR}"
  exit 1
fi

if [ -z "$KEEP_SIEVETEST" ]
then
  rm -f "${FB}" "${RELS}" "${FBC}"
  rmdir "${TMPDIR}"
else
  echo "Keeping files in ${TMPDIR}"
fi
