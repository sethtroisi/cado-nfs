#!/bin/sh

# This script takes as first input a tar file containing one bwc test. It
# then runs the contained test.sh script, with the environment variable
# $tmp set to the temporary working area, and the variable $bwc pointing
# to the place where bwc binaries are supposed to reside.

tmp=`mktemp -d "/tmp/bwc-test-XXXX"`
export tmp
bwc="`dirname $0`"
export bwc

(cd "$tmp" ; tar xzf - ) < "$1"

"$tmp/test.sh" "$@"

rc=$?

rm -rf $tmp

exit $?
