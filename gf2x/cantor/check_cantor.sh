#!/bin/sh

for i in "$@" ; do
  ./cantor $i > /tmp/toto
  if `magma < cantor128.m | grep true`; then
    echo "ok $i"
  else
    echo "failed $i"
    exit 1
  fi
done

