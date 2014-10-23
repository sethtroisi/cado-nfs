#!/usr/bin/env bash

set -e
export DIR=`mktemp -d /tmp/cado-cov.XXXXXXXXXXXX`
export TMPDIR="$DIR"
export CADO_DEBUG=1
export COV=1
F=$TMPDIR/local.sh
cat > $F <<EOF
CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
CXXFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
EOF
export LOCALFILE=$F
$HOME/NFS/cado/scripts/nightly-test all check
cd $DIR
geninfo --no-checksum --ignore-errors gcov,source -q --output-filename $DIR/cado-nfs.info  ./ --no-external
rm -rf ~/.webdir/cado/ || :
genhtml   -o ~/.webdir/cado/ $DIR/cado-nfs.info
rm -rf $DIR
