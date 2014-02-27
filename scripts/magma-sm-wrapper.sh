#!/bin/bash

unset POLY
unset PURGED
unset INDEX
unset OUT
unset ELL
unset SMEXP
unset NMAPS
unset MT

while [ -n "$1" ]
do
  if [ "$1" = "-poly" ]
  then
    POLY="$2"
    shift 2
  elif [ "$1" = "-purged" ]
  then
    PURGED="$2"
    shift 2
  elif [ "$1" = "-index" ]
  then
    INDEX="$2"
    shift 2
  elif [ "$1" = "-out" ]
  then
    OUT="$2"
    shift 2
  elif [ "$1" = "-gorder" ]
  then
    ELL="$2"
    shift 2
  elif [ "$1" = "-smexp" ]
  then
    SMEXP="$2"
    shift 2
  elif [ "$1" = "-nsm" ]
  then
    NMAPS="$2"
    shift 2
  elif [ "$1" = "-mt" ]
  then
    MT="$2"
    shift 2
  else
    echo "Unknown parameter: $1"
    exit 1
  fi
done

# where am I ?
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CMD="$DIR/../filter/sm -poly $POLY -purged $PURGED -index $INDEX -out $OUT -gorder $ELL -smexp $SMEXP -nsm $NMAPS -mt $MT"

echo $CMD

$CMD
