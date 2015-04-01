#!/usr/bin/env bash

set -e
# set -x

# various configuration variables. environment can be used to override them

: ${scriptpath=$0}
: ${m=128}
: ${n=128}
: ${wdir=/tmp/bwc}
: ${mm_impl=bucket}
: ${nullspace=left}

export m n scriptpath wdir mm_impl prime=2 nullspace
exec "`dirname $0`/bwc-ptrace.sh" "$@"
