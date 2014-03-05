#!/bin/bash

unset POLY
unset RENUMBER
unset BADIDEALINFO
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
  elif [ "$1" = "-renumber" ]
  then
    RENUMBER="$2"
    shift 2
  elif [ "$1" = "-badidealinfo" ]
  then
    BADIDEALINFO="$2"
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
    # operates in 3 steps:
    # 1. extract (p, r) ideals from $RENUMBER and put them in ...algpr.gz
    # 2. from ...algpr.gz, compute a generator for each ideal and put it in
    #    ...generators.gz
    # 3. compute the unit contribution for each relation set from $INDEX
    #
    # algpr.gz and generators.gz should be computed only once and can be
    # done in the meantime for instance.
    #
    # TODO: it might be possible to skip the algpr file (together with the
    # debug_renumber prgm) by rewriting / using this more cleverly in C.

    prgm=$DIR/../misc/debug_renumber

    if [ ! -s $prgm ]; then echo "Please compile $prgm"; exit; fi

    # guess names
    algpr=`echo $RENUMBER | sed 's/freerel.renumber/algpr/'`
    generators=`echo $RENUMBER | sed 's/freerel.renumber/generators/'`

    if [ ! -s $algpr ]; then
	echo "Building file $algpr using debug_renumber"
	$prgm -poly $POLY -renumber $RENUMBER |\
        grep " alg " | sed 's/alg side//' | sed 's/=/ /g' |\
        awk '{print $2, $6, $8}' |\
        gzip -c > $algpr
    fi
    if [ ! -s $generators ]; then
	echo "Building file $generators using Magma"
	CMD="magma -b polyfile:=$POLY renumber:=$RENUMBER algpr:=$algpr generators:=$generators badidealinfo:=$BADIDEALINFO purged:=$PURGED index:=$INDEX ficunits:=$OUT ww:=true $DIR/nfsunits.mag"
	echo $CMD; $CMD
    fi

    echo "## Using magma to compute units (version 2)"
    # results are put in $OUT
    CMD="magma -b polyfile:=$POLY renumber:=$RENUMBER algpr:=$algpr generators:=$generators badidealinfo:=$BADIDEALINFO purged:=$PURGED index:=$INDEX ficunits:=$OUT ww:=false $DIR/nfsunits.mag"
fi

echo $CMD

$CMD
