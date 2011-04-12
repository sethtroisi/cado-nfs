#!/bin/sh

# This script must be present on the remote site. It achieves the
# resetting of the source tree, and only later calls the client script.

set -e
set -x

reference=origin

if [ "$1" = "--git-reference" ] ; then
    shift
    reference="$1"
    shift
fi

cd "$HOME/cado"
while ! git reset --hard HEAD  ; do echo "waiting for git" ; sleep 1 ; done
while ! git pull --force $reference ; do echo "waiting for git" ; sleep 1 ; done
while ! git reset --hard HEAD ; do echo "waiting for git" ; sleep 1 ; done

exec "$@"
