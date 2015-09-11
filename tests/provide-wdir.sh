#!/usr/bin/env bash

# This shell script will do the following.
#
# - create a temp directory.
# - run the given command line, by providing the infomration about the
#   location of the temp directory in any desired way.
# - clean up the temp directory

# If CADO_DEBUG is set to something, the temp directory is not cleaned.


t=`mktemp -d /tmp/XXXXXXXXXXXXXX`

extra=()

while [ $# -gt 0 ] ; do
    if [ "$1" = "--" ] ; then
        shift
        break;
    elif [ "$1" = "--env" ] ; then
        shift
        var="$1"
        eval "$var=$t"
        eval "export $var"
        shift
    elif [ "$1" = "--arg" ] ; then
        shift
        extra=("${extra[@]}" "$1=$t")
        shift
    else
        # Then we finish processing, presumably the user omitted the --
        # separator.
        break
    fi
done

"$@" "${extra[@]}"

rc=$?

if [ "$CADO_DEBUG" ] ; then
    echo "debug mode, data left in $t"
    exit $rc
fi

if [ $rc != 0 ] ; then
    # what do we do when the script failed ?
    rm -rf "$t"
else
    rm -rf "$t"
fi

exit $rc
