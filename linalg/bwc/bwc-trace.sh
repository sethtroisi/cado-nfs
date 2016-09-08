#!/usr/bin/env bash

# NOTE: magma debug was in bf225ce125bc21c5ae5f4c3b2f8acc652e248f46, and
# it seems it has not been ported completely.

set -e
# set -x

# various configuration variables. environment can be used to override them

: ${scriptpath=$0}
: ${m=128}
: ${n=128}
: ${wdir=/tmp/bwc}
# : ${mm_impl=bucket}
: ${nullspace=left}

export m n scriptpath wdir mm_impl prime=2 nullspace
exec "`dirname $0`/bwc-ptrace.sh" "$@"
