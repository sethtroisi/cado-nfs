#!/bin/sh

# usage:
# ./tryonedep.sh -depnum 0 matfile kerfile polyfile
# This creates temp files in /tmp with fixed names (sorry...) If they
# already exist and you cannot overwrite them, this won't work.


echo "*** Running rational side sqrt..."
./ratsqrt $@ > /tmp/ratside

echo "*** Computing algebraic side square..."
./algsqrt $@ > /tmp/algside

echo "*** Computing algebraic side sqrt..."
magma < finish.mag

rm /tmp/ratside /tmp/algside /tmp/formagma
