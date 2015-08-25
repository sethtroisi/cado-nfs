#!/bin/sh

CADO_NFS_SOURCE_DIR=$1

NCPUS=$("`dirname $0`"/ncpus.sh)
# we cap tasks.threads to 2, because las and polyselect behave quite
# badly for larger values, especially for smallish tests.
export t=`mktemp -d /tmp/cado-check.XXXXXXX`
${CADO_NFS_SOURCE_DIR}/factor.sh 43341748620473677010074177283795146221310971425909898235183 -dlp -t 2 slaves.nrclients=$((NCPUS/2)) tasks.linalg.bwc.threads=$NCPUS && rm -rf $t
