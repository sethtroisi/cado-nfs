#!/usr/bin/env bash

# Print and run a command
function run {
  echo "Running: $*"
  "$@"
}

# Input files. First parameter is the polynomial file to use, second is root of the build directory, third is the expected SHA1 value
FB="$1"
LAS="$2"
POLY="$3"
REFERENCE_SHA1="$4"
REFERENCE_REVISION="$5"
CHECKSUM_FILE="$6"
shift 6

if [[ -z "${FB}" || ! -f "${FB}" ]]
then
  echo "Factor base ${FB} is not a file" >&2
  exit 1
fi

if [[ -z "${LAS}" || ! -x "${LAS}" ]]
then
  echo "Las binary ${LAS} is not an executable file" >&2
  exit 1
fi

if [ -z "${REFERENCE_SHA1}" ]
then
  echo "No reference SHA1 value specified" >&2
  exit 1
fi


if [[ -z "$lim0" || -z "$lim1" || -z "$lpb0" || -z "$lpb1" || -z "$maxbits" \
      || -z "$mfb0" || -z "$mfb1" || -z "$lambda0" || -z "$lambda1" || -z "$I" \
      || -z "$q0" || ( -z "$q1" && -z "$rho" ) ]]
then
  echo "Required shell environment variable not set" >&2
  exit 1
fi

if [ "$1" = "-regex" ]
then
  REGEX="$2"
  shift 2
fi

BASENAME="`basename "${POLY}"`"
BASENAME="${BASENAME%%.*}"

# Output files
WORKDIR=`mktemp -d ${TMPDIR-/tmp}/cadotest.XXXXXXXXXX`
# Make temp direcotry world-readable for easier debugging
chmod a+rx "${WORKDIR}"
RELS="${WORKDIR}/${BASENAME}.rels"
FBC="${WORKDIR}/${BASENAME}.fbc"

if [ -n "$q1" ]
then
  end=("-q1" "$q1")
fi
if [ -n "$rho" ]
then
  end=("-rho" "$rho")
fi

bailout() {
    if [ "$EXPECTED_FAIL" ] && ! [ "$KEEP_SIEVETEST" ] ; then
        rm -rf "$WORKDIR"
    fi
    exit 1
}

# first exercise the -fbc command-line option to create a cache file
run "$LAS" -poly "$POLY" -fb "${FB}" -I "$I" -lim0 "$lim0" -lpb0 "$lpb0" -mfb0 "$mfb0" -lambda0 "$lambda0" -lim1 "$lim1" -lpb1 "$lpb1" -mfb1 "$mfb1" -lambda1 "$lambda1" -q0 "$q0" -q1 "$q0" -out "${RELS}" -fbc "${FBC}" "$@" || bailout
# then use the cache file created above
run "$LAS" -poly "$POLY" -fb "${FB}" -I "$I" -lim0 "$lim0" -lpb0 "$lpb0" -mfb0 "$mfb0" -lambda0 "$lambda0" -lim1 "$lim1" -lpb1 "$lpb1" -mfb1 "$mfb1" -lambda1 "$lambda1" -q0 "$q0" "${end[@]}" -out "${RELS}" -fbc "${FBC}" "$@" || bailout


SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
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
  echo "$0: Got SHA1 of ${SHA1} but expected ${REFERENCE_SHA1}${REFMSG}. Files remain in ${WORKDIR}"
  exit 1
fi

if [ -n "${CHECKSUM_FILE}" ]
then
  MYCHECKSUM_FILE="${WORKDIR}/${BASENAME}.checksums"
  grep "# Checksums over sieve region:" "${RELS}" > "${MYCHECKSUM_FILE}"
  if [ -f "${CHECKSUM_FILE}" ]
  then
    # File with checksums already exists, compare
    if diff -b "${CHECKSUM_FILE}" "${MYCHECKSUM_FILE}" > /dev/null
    then
      echo "Checksums agree"
    else
      echo "Error, reference checksums in ${CHECKSUM_FILE} differ from mine in ${MYCHECKSUM_FILE}" >&2
      exit 1
    fi
  else
    # File with checksums does not exists, create it
    cp "${MYCHECKSUM_FILE}" "${CHECKSUM_FILE}"
    echo "Created checksum file"
  fi
fi

if [ -n "${REGEX}" ]
then
  echo "Searching for regex \"${REGEX}\"" >&2
  if ! grep "${REGEX}" "${RELS}" >&2
  then
    echo "Error, regular expression \"${REGEX}\" does not match output file"
    exit 1
  fi
fi

if [ -z "$KEEP_SIEVETEST" ]
then
  rm -f "${RELS}" "${FBC}"
  if [ -n "${MYCHECKSUM_FILE}" ]
  then
    rm -f "${MYCHECKSUM_FILE}"
  fi
  rmdir "${WORKDIR}"
else
  echo "Keeping files in ${WORKDIR}"
fi
