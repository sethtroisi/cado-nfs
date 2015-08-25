#!/bin/sh

CADO_NFS_SOURCE_DIR=$1

NCPUS=$("`dirname $0`/ncpus.sh)
export t=`mktemp -d /tmp/cado-check.XXXXXXX`
${CADO_NFS_SOURCE_DIR}/factor.sh 427545955933445497303117859123 -dlp -ell 24982234190338056404295773 -t $NCPUS target=62720951885220719832970790782  && rm -rf $t
