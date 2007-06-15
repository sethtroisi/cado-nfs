#!/bin/sh

echo "*** Running rational side sqrt..."
./ratsqrt $@ > /tmp/ratside

echo "*** Computing algebraic side square..."
./algsqrt $@ > /tmp/algside

echo "*** Computing algebraic side sqrt..."
magma < finish.mag

rm /tmp/ratside /tmp/algside /tmp/formagma
