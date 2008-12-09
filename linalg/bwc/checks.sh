#!/bin/bash

# Usage: ./checks.sh < checks.list ; or ./checks.sh < checks.list.light

# Note that the typical failure situation yields a job which never
# finishes.
ntests=0
while read cmd ; do
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
