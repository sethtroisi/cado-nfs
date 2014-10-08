#!/usr/bin/env bash

unset POLY
unset RENUMBER
unset BADIDEALINFO
unset PURGED
unset INDEX
unset OUT
unset ELL
unset SMEXP0
unset SMEXP1
unset NMAPS0
unset NMAPS1
unset ABUNITS
unset MT

EXPLICIT0="no"
EXPLICIT1="no"

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
  elif [ "$1" = "-smexp0" ]
  then
    SMEXP0="$2"
    shift 2
  elif [ "$1" = "-smexp1" ]
  then
    SMEXP1="$2"
    shift 2
  elif [ "$1" = "-nsm0" ]
  then
    NMAPS0="$2"
    shift 2
  elif [ "$1" = "-nsm1" ]
  then
    NMAPS1="$2"
    shift 2
  elif [ "$1" = "-abunits" ]
  then
    ABUNITS="$2"
    shift 2
  elif [ "$1" = "-t" ]
  then
    MT="$2"
    shift 2
  elif [ "$1" = "-explicit_units0" ]
  then
    EXPLICIT0="yes"
    shift 1
  elif [ "$1" = "-explicit_units1" ]
  then
    EXPLICIT1="yes"
    shift 1
  else
    echo "Unknown parameter: $1"
    exit 1
  fi
done

# where am I ?
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## set parameters
side=-1
if [ $EXPLICIT0 = "yes" ]; then
    smopts="-smexp0 $SMEXP0 -nsm0 0"
    side=0
    usef=false
else
    smopts="-smexp0 $SMEXP0 -nsm0 $NMAPS0"
fi
if [ $EXPLICIT1 = "yes" ]; then
    smopts="$smopts -smexp1 $SMEXP1 -nsm1 0"
    side=1
    usef=true
else
    smopts="$smopts-smexp1 $SMEXP1 -nsm1 $NMAPS1"
fi

## build required SM's
if [ $NMAPS0 -gt 0 -o $NMAPS1 -gt 0 ]; then
    if [ -s $OUT ]; then
	echo "File $OUT already exists"
    else
	CMD="$DIR/../filter/sm -poly $POLY -purged $PURGED -index $INDEX -out $OUT -gorder $ELL $smopts -t $MT"
	echo $CMD; $CMD
    fi
else
    echo "No SM's to be computed"
fi

## now, use units if needed: we suppose that we use either 0 or 1, not both
if [ $side -ne -1 ]; then
    # operates in 3 steps:
    # 1. extract (p, r) ideals from $RENUMBER and put them in ...algpr.gz
    # 2. from ...algpr.gz, compute a generator for each ideal and put it in
    #    ...generators.gz
    # 3. compute the unit contribution for each relation set from $INDEX
    #
    # algpr.$side.gz and generators.$side.gz should be computed only once and
    # can be done in the meantime for instance.
    #
    # TODO: it might be possible to skip the algpr file (together with the
    # debug_renumber prgm) by rewriting / using this more cleverly in C.

    prgm=$DIR/../misc/debug_renumber

    if [ ! -s $prgm ]; then echo "Please compile $prgm"; exit; fi

    # guess names
    algpr=`echo $RENUMBER | sed "s/freerel.renumber/algpr.$side/"`
    generators=`echo $RENUMBER | sed "s/freerel.renumber/generators.$side/"`
    relsdel=`echo $PURGED | sed 's/purged/relsdel/'`
    abunits=`echo $RENUMBER | sed "s/freerel.renumber.gz/units.abunits.$side/"`

    # extract (p, r)
    if [ ! -s $algpr ]; then
	echo "Building file $algpr for side $side using debug_renumber"
	# we have to make it work for rat/alg or side 0/1
	$prgm -poly $POLY -renumber $RENUMBER |\
        egrep "(alg side|side $side)" | sed 's/=/ /g' |\
        awk '{print $2, $6, $8}' |\
        gzip -c > $algpr
    fi
    # find generators for (p, r)
    if [ ! -s $generators ]; then
	echo "Building file $generators using Magma for side $side"
	CMD="magma -b polyfile:=$POLY renumber:=$RENUMBER algpr:=$algpr generators:=$generators badidealinfo:=$BADIDEALINFO purged:=$PURGED relsdel:=$relsdel index:=$INDEX ficunits:=$OUT abunits:=$abunits usef:=$usef ww:=true $DIR/nfsunits.mag"
	echo $CMD; $CMD
    fi

    if [ -s $OUT.$side ]; then
	echo "File $OUT.$side already exists"
    else
	echo "## Using magma to compute units for side $side"
        # results are put in $OUT.$side
	CMD="magma -b polyfile:=$POLY renumber:=$RENUMBER algpr:=$algpr generators:=$generators badidealinfo:=$BADIDEALINFO purged:=$PURGED relsdel:=$relsdel index:=$INDEX ficunits:=$OUT.$side abunits:=$abunits usef:=$usef ww:=false $DIR/nfsunits.mag"
	echo $CMD; $CMD

	# split abunits
	# typicaly ABUNITS="...../p3dd15-f4g3-GJL-1.sm.abunits.dir"
	abunitsdir=$ABUNITS".$side"
	if [ -d $abunitsdir ]; then
	    echo "Directory $abunitsdir already exists"
	else
	    echo "I create directory $abunitsdir"
	    sort -k1n -k2n $abunits > /tmp/foo$$
	    mkdir $abunitsdir
	    echo "I split /tmp/foo$$ into files + index"
	    $DIR/split.py /tmp/foo$$ $abunitsdir/ 100
	    /bin/rm -f /tmp/foo$$
	fi

        # paste files if needed
	if [ $side -eq 1 ]; then
	    # paste $OUT.$side and $OUT to get a new $OUT
	    f1=$OUT.$side; f2=$OUT
	else
	    # paste $OUT and $OUT.$side to get a new $OUT
	    f1=$OUT; f2=$OUT.$side
	fi
	paste -d" " $f1 $f2 | awk '{n++; if(n==1){print $1}else{print}}' > $OUT.$$
	mv $OUT.$$ $OUT
    fi
fi
