#!/bin/sh

if [ -z "$CANTOR" ] ; then CANTOR=./cantor ; fi

for i in "$@" ; do
  "$CANTOR" --test $i > /tmp/toto
  if `magma < debug/check_cantor.m | grep true`; then
    echo "ok $i"
  else
    echo "failed $i"
    exit 1
  fi
done

