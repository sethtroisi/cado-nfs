#!/usr/bin/env bash

BINDIR="$1"
TESTBENCH="$1/sieve/ecm/testbench"
INPUTFILE="$2"
shift 2

REQUIRED_OUTPUT="`mktemp`"
INPUTNUMBERS="`mktemp`"
# First word on each line is the input number
sed 's/ *#.*$//' < "${INPUTFILE}" | grep . | cut -d " " -f 1 > "${INPUTNUMBERS}" || exit 1
# Remaining words are the required output
sed 's/ *#.*$//' < "${INPUTFILE}" | grep . | cut -d " " -f 2- > "${REQUIRED_OUTPUT}" || exit 1

cd $BINDIR
make testbench

ACTUAL_OUTPUT="`mktemp`"
"${TESTBENCH}" -inp "${INPUTNUMBERS}" "$@" > "${ACTUAL_OUTPUT}"

if ! cmp "${REQUIRED_OUTPUT}" "${ACTUAL_OUTPUT}"
then
  echo "testbench produced output in file \"${ACTUAL_OUTPUT}\", but expected result as in \"${REQUIRED_OUTPUT}\", input numbers are in \"${INPUTNUMBERS}\""
  diff "${REQUIRED_OUTPUT}" "${ACTUAL_OUTPUT}"
  exit 1
fi

rm -f "${REQUIRED_OUTPUT}" "${INPUTNUMBERS}" "${ACTUAL_OUTPUT}"
exit 0
