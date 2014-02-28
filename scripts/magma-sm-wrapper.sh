#!/bin/bash

unset POLY
unset PURGED
unset INDEX
unset OUT
unset ELL
unset SMEXP
unset NMAPS
unset MT
## for units
unset RELSDIR
EXPLICIT="no"

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
  elif [ "$1" = "-rels" ]
  then
    RELSDIR="$2"
    shift 2
  elif [ "$1" = "-explicit_units" ]
  then
    EXPLICIT="yes"
    shift 1
  else
    echo "Unknown parameter: $1"
    exit 1
  fi
done

# where am I ?
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $EXPLICIT == "no" ]; then
    CMD="$DIR/../filter/sm -poly $POLY -purged $PURGED -index $INDEX -out $OUT -gorder $ELL -smexp $SMEXP -nsm $NMAPS -mt $MT"
else
## RELSDIR=....../p59.sieving.204000-205000.1_1jgx.gz
## make it ....../p59
    RELSDIR=`echo $RELSDIR | sed 's/\.sieving.*//'`
    CMD="magma -b polyfile:=$POLY purged:=$PURGED index:=$INDEX ficunits:=$OUT relsdir:=$RELSDIR $DIR/nfsunits.mag" 
fi

echo $CMD

$CMD
