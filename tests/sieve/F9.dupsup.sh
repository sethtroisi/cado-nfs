#!/usr/bin/env bash

if [ "$CADO_DEBUG" ] ; then
    set -x
fi
set -e

RELSD="${WORKDIR:?missing}/rels_dupsup"
RELSR="${WORKDIR:?missing}/rels_raw"
RELSF="${WORKDIR:?missing}/rels_filtered"

RELS=$RELSD "`dirname $0`/sievetest.sh" -dup "$@"
ndup_ref=`grep -wc DUPE $RELSD`
# remove leading spaces (for openbsd 5.3)
let ndup_ref=ndup_ref
if [ "$ndup_ref" != "$DUPE_COUNT" ]; then
    echo "Wrong number of duplicates in las -dup"
    echo "Expected $DUPE_COUNT, got $ndup_ref"
    exit 1
fi


RELS=$RELSR "`dirname $0`/sievetest.sh" "$@"

LAS_BINARY=$DUPSUP_BINARY RELS=$RELSF "`dirname $0`/sievetest.sh" -skew 1.0 "${RELSR}"
ndup=`grep -wc DUPE $RELSF`
# remove leading spaces (for openbsd 5.3)
let ndup=ndup
if [ "$ndup" != "$DUPE_COUNT" ]; then
    echo "Wrong number of duplicates in dupsup"
    echo "Expected $DUPE_COUNT, got $ndup"
    exit 1
fi
