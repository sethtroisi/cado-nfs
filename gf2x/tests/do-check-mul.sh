#!/bin/sh

cat "$srcdir/check-mul.res" | while read n1 n2 v s ; do
    echo -n "${n1}x${n2} "
    got=`./check-mul $n1 $n2`
    expected="$n1 $n2 $v $s"
    if [ "$got" != "$expected" ] ; then
        echo "failed check for ${n1}x${n2} : '$got' != '$expected'" >&2
        echo "failed : '$got' != '$expected'"
        exit 1
    fi
done
rc=$?
echo
exit $rc
