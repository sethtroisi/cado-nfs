#!/bin/sh

CADO_NFS_SOURCE_DIR=$1

NCPUS=$("`dirname $0`/ncpus.sh)
export t=`mktemp -d /tmp/cado-check.XXXXXXX`
${CADO_NFS_SOURCE_DIR}/factor.sh 100000000000000000039 -dlp -gfpext 2 -t $NCPUS -ell 164354743277891 && rm -rf $t
