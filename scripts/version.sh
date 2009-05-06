#!/bin/sh

# This displays as much information as might be relevant, either within
# an svn checkout or from my local git-svn managed branch.

cd "`dirname $0`"
if [ -d .svn ] ; then
    svnversion
    # svn status -q # removed, since it pollutes the revision number
else
    commit="`git show --pretty=format:%h 2>/dev/null | head -1`"
    if [ "$commit" != "" ] ; then
        # ver="`git svn log | head -n2 | tail -n1 | cut -d\  -f1`"
        # much faster
        perl='while (<>) { /trunk@(\d+)/ && print "$1\n"; }'
        ver=`git show refs/remotes/trunk 2>/dev/null | perl -e "$perl"`
        if [ "$ver" ] ; then
            echo -n "svn $ver -- git $commit"
        else
            echo -n "git $commit"
        fi
        # git diff --name-status is svn status -q
        if [ "`git diff --name-status`" != "" ] ; then
            echo " +mods"
        else
            echo
        fi
    else
        echo "exported"
    fi
fi
