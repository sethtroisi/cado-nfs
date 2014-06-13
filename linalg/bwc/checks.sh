#!/usr/bin/env bash

# Usage: ./checks.sh < checks.list ; or ./checks.sh < checks.list.light

# Note that the typical failure situation yields a job which never
# finishes.

if [ -f Makefile.local ] ; then
    source Makefile.local
    if [ -n "$MPI" ] ; then
        export PATH="$MPI:$PATH"
    fi
fi

# This is for mpich2 only.
if [ -z "$HYDRA_HOST_FILE" ] ; then
    export HYDRA_USE_LOCALHOST=1
fi

ntests=0
while read cmd ; do
    if [ -z "$cmd" ] ; then continue; fi
    let ntests+=1
    echo -n "$cmd --> "
    if sh -c "$cmd" < /dev/null | grep OK ; then
        :
    else
        echo "Test $ntests FAILED"
        exit 1
    fi
done

echo "All $ntests tests passed"
