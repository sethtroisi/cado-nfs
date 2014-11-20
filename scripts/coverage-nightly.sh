#!/usr/bin/env bash

export DIR=`mktemp -d /tmp/cado-cov.XXXXXXXXXXXX`
F=`mktemp /tmp/cado-cov.local.XXXXXXXXX.sh`

doit() {
    set -e
    export TMPDIR="$DIR"
    export CADO_DEBUG=1
    export COV=1
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
}

rc=0
if ! doit ; then
    echo "FAILED !!!" >&1
    rc=1
fi

rm -rf $DIR $F
exit $rc

