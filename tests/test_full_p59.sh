#!/bin/sh

CADO_NFS_SOURCE_DIR=$1

NCPUS=$("`dirname $0`/ncpus.sh)
export t=`mktemp -d /tmp/cado-check.XXXXXXX`
${CADO_NFS_SOURCE_DIR}/factor.sh 43341748620473677010074177283795146221310971425909898235183 -dlp -t $NCPUS && rm -rf $t
