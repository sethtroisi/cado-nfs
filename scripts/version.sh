#!/bin/sh

# This displays as much information as might be relevant, either within
# an svn checkout or from my local git-svn managed branch.

cd "`dirname $0`"

if [ "$1" = "--svn" ] ; then
    if [ -d .svn ] ; then
        svnversion .| cut -d\: -f2
    else
        commit="`git show --pretty=format:%h 2>/dev/null | head -1`"
        if [ "$commit" != "" ] ; then
            perl='while (<>) { /trunk@(\d+)/ && print "$1\n"; }'
            ver=`git show refs/remotes/trunk 2>/dev/null | perl -e "$perl"`
            if [ "$ver" ] ; then
                echo "$ver"
            else
                echo "$commit"
            fi
        fi
    fi
else
    if [ -d .svn ] ; then
        svnversion .
        # svn status -q # removed, since it pollutes the revision number
    else
        commit="`git show --pretty=format:%h 2>/dev/null | head -1`"
        if [ "$commit" != "" ] ; then
            # ver="`git svn log | head -n2 | tail -n1 | cut -d\  -f1`"
            # much faster
            perl='while (<>) { /trunk@(\d+)/ && print "$1\n"; }'
            ver=`git show refs/remotes/trunk 2>/dev/null | perl -e "$perl"`
            if [ "$ver" ] ; then
                commit="svn $ver -- git $commit"
            else
                commit="git $commit"
            fi
            # git diff --name-status is svn status -q
            # Note that checked out copies which have GIT_DIR set for
            # some reason could end up triggering errors on git diff
            if [ "`git diff --name-status 2>/dev/null`" != "" ] ; then
                echo "$commit +mods"
            else
                echo "$commit"
            fi
        else
            echo "exported"
        fi
    fi
fi
