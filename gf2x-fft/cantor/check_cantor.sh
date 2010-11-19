#!/bin/sh

if [ -z "$CANTOR" ] ; then CANTOR=./cantor ; fi

for i in "$@" ; do
  "$CANTOR" $i > /tmp/toto
  if `magma < check_cantor.m | grep true`; then
    echo "ok $i"
  else
    echo "failed $i"
    exit 1
  fi
done

