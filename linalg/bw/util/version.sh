#!/bin/sh

# This displays as much information as might be relevant, either within
# an svn checkout or from my local git-svn managed branch.

cd `dirname $0`
if [ -d .svn ] ; then
    svnversion
    svn status -q
else
    commit=`git-show --pretty=format:%h 2>/dev/null | head -1`
    if [ $commit != "" ] ; then
        ver=`git svn log | head -n2 | tail -n1 | cut -d\  -f1`
        echo "svn $ver -- git $commit"
        git-diff --name-status
    else
        echo "exported"
    fi
fi
