#!/bin/sh

# This displays as much information as might be relevant. A trailing +
# indicates that local modifications exist.

cd "`dirname $0`"

commit="`git show --pretty=format:%h 2>/dev/null | head -1`"
if [ "$commit" != "" ] ; then
    # git diff --name-status is svn status -q
    # Note that checked out copies which have GIT_DIR set for
    # some reason could end up triggering errors on git diff
    if [ "`git diff --name-status 2>/dev/null`" != "" ] ; then
        echo "${commit}+"
    else
        echo $commit
    fi
else
    head -1 ../.git_version
fi
