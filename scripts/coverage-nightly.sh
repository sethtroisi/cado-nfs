#!/usr/bin/env bash

: ${WEBDIR=$HOME/.webdir}
: ${TARGET_DIRECTORY=$WEBDIR/cado}
export DIR=`mktemp -d /tmp/cado-cov.XXXXXXXXXXXX`
F=`mktemp /tmp/cado-cov.local.XXXXXXXXX.sh`

while [ $# -gt 0 ] ; do
    if [ "$1" = "--thorough" ] ; then
        CHECKS_EXPENSIVE=1
    elif [ "$1" = "--target-directory" ] ; then
        shift
        TARGET_DIRECTORY="$1"
    else
        echo "Unexpected argument: $1" >&2
        exit 1
    fi
    shift
done


doit() {
    set -e
    export TMPDIR="$DIR"
    export CADO_DEBUG=1
    export COV=1
    cat > $F <<EOF
    CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
    CXXFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
EOF
    if [ "$CHECKS_EXPENSIVE" ] ; then
        cat >> $F <<EOF
CHECKS_EXPENSIVE=1
EOF
    fi
    export LOCALFILE=$F
    $HOME/NFS/cado/scripts/nightly-test all check
    cd $DIR
    geninfo --no-checksum --ignore-errors gcov,source -q --output-filename $DIR/cado-nfs.info  ./ --no-external
    rm -rf "$TARGET_DIRECTORY"/ || :
    genhtml   -o "$TARGET_DIRECTORY"/ $DIR/cado-nfs.info
}

rc=0
if ! doit ; then
    echo "FAILED !!!" >&1
    rc=1
fi

rm -rf $DIR $F

rsync -a  "$TARGET_DIRECTORY"/*png "$TARGET_DIRECTORY"/gcov.css "$TARGET_DIRECTORY"/index.html  "$TARGET_DIRECTORY-$(date +%Y%m%d)"/  || :

exit $rc

