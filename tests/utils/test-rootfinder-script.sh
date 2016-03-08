#!/usr/bin/env bash

# This script eats a testbase on stdin.  Normally it gives nothing on
# stdout, and gives progress report on stderr. If it prints anything,
# that means bug.

# testbases are built for examble by the magma script, which may be
# modified.

n=0
nok=0
while read in x ; do
    y=$(./test-rootfinder `echo $x`)
    read out z
    if [ "$y" != "$z" ] ; then
        (echo in $x ; echo exp $z ; echo got $y)
    else
        let nok+=1
    fi
    let n+=1
    if [ `expr $n % 100` = 0 ] ; then
        echo -e "# OK $nok/$n\r" >&2
    fi
done

# That's our return code.
[ $nok = $n ]
