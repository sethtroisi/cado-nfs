#!/bin/bash

unset ELL
unset NMAPS
unset SPARSE
unset SM
unset KER

while [ -n "$1" ]
do
  if [ "$1" = "-ell" ]
  then
    ELL="$2"
    shift 2
  elif [ "$1" = "-nmaps" ]
  then
    NMAPS="$2"
    shift 2
  elif [ "$1" = "-sparsemat" ]
  then
    SPARSE="$2"
    shift 2
  elif [ "$1" = "-sm" ]
  then
    SM="$2"
    shift 2
  elif [ "$1" = "-ker" ]
  then
    KER="$2"
    shift 2
  else
    echo "Unknown parameter!"
    exit 1
  fi
done

# where am I ?
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if true; then
CMD="magma ell:=$ELL nmaps:=$NMAPS sparsefile:=$SPARSE smfile:=$SM kerfile:=$KER $DIR/linalg.mag"

echo $CMD

$CMD

else
./linalg4catrel -ell $ELL -sparse $SPARSE -sm $SM -mode sparse -nthreads 2 > $KER
fi
