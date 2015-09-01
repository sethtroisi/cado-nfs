#!/bin/bash

CHECK_RELS="$1"
SOURCE_TEST_DIR="`dirname "$0"`"
TMPDIR=`mktemp -d /tmp/cadotest.XXXXXXXX`
chmod a+rx "${TMPDIR}"

poly="${SOURCE_TEST_DIR}/c60.poly"
rels1="${SOURCE_TEST_DIR}/c60.rels1"
rels2="${SOURCE_TEST_DIR}/c60.rels2"

out="$TMPDIR/out"
fixed="$TMPDIR/fixed"


############## Check with correct relations
if ! ("$CHECK_RELS" -poly $poly -check_primality $rels1 > $out 2>&1) ; then
    echo "$0: check_rels binary failed to recognize correct relations"
    exit 1
fi

nrels=`awk '/read relations/ {print $5}' $out`
correctrels=`awk '/correct relations/ {print $5}' $out`
if [ $nrels != $correctrels ]; then
    echo "$0: check_rels binary failed to recognize correct relations"
    exit 1
fi


############## Check with wrong relations
if ! ("$CHECK_RELS" -poly $poly -check_primality -complete $fixed $rels2 > $out 2>&1) ; then
    echo "$0: check_rels binary failed when testing wrong relations"
    exit 1
fi

nrels=`awk '/read relations/ {print $5}' $out`
correctrels=`awk '/keeped relations/ {print $5}' $out`
composite=`awk '/one non-primes ideal/ {print $2}' $out`
if [ $nrels != $correctrels ]; then
    echo "$0: check_rels binary failed to fix wrong relations"
    exit 1
fi
if [ $composite != "1" ]; then
    echo "$0: check_rels binary failed to detect composite factor"
    exit 1
fi

rm -rf $TMPDIR
