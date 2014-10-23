#!/bin/sh

set -e
set -x

script="$1" ; shift

: ${TMPDIR:=/tmp}

wdir=`mktemp -d $TMPDIR/bwc-test.XXXXXXXXXX`
matdir=`mktemp -d $TMPDIR/bwc-test-mats.XXXXXXXX`

$script mats=$matdir wdir=$wdir "$@"

if ! [ "$CADO_DEBUG" ] ; then
    rm -rf $wdir
    rm -rf $matdir
fi

