#!/bin/sh

CADO_NFS_BUILD_DIR=$1

cd $CADO_NFS_BUILD_DIR/gf2x
make check
