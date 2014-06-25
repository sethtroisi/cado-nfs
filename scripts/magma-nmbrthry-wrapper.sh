#!/usr/bin/env bash

unset P
unset ELL
unset POLY
unset BAD
unset BADINFO

while [ -n "$1" ]
do
  if [ "$1" = "-ell" ]
  then
    ELL="$2"
    shift 2
  elif [ "$1" = "-poly" ]
  then
    POLY="$2"
    shift 2
  elif [ "$1" = "-p" ]
  then
    P="$2"
    shift 2
  elif [ "$1" = "-badfile" ]
  then
    BAD="$2"
    shift 2
  elif [ "$1" = "-badinfofile" ]
  then
    BADINFO="$2"
    shift 2
  else
    echo "Unknown parameter!"
    exit 1
  fi
done

# ell=0 means undefined, please provide one valid ell
: ${ELL:="0"}

rm -f $BAD
rm -f $BADINFO

# where am I ?
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CMD="magma -b ell:=$ELL mainp:=$P polyfile:=$POLY badfile:=$BAD badinfofile:=$BADINFO $DIR/badideals.mag" 

echo $CMD >&2

$CMD
