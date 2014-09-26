#!/bin/sh

set -e
set -x

script="$1" ; shift

wdir=`mktemp -d /tmp/bwc-test.XXXXXXXXXX`
matdir=`mktemp -d /tmp/bwc-test-mats.XXXXXXXX`

$script mats=$matdir wdir=$wdir "$@"

if ! [ "$CADO_DEBUG" ] ; then
    rm -rf $wdir
    rm -rf $matdir
fi

