#!/usr/bin/env bash

if [ -z "$1" ]
then
  echo "Specify which polynomial to process [a-h]"
  exit 1
fi

WHICH_POLY="$1"
KEEP=100
ADMAX=400000
MINP=1000
STEPP=100
MAXP=20000

if [ -z "${POLYSELECT2L}" ]
then
  POLYSELECT2L="$HOME/bin/polyselect2l.33f2069"
  echo "Using default polyselect2l at ${POLYSELECT2L}"
fi

# Check if an output file is needed, i.e., does not exist, and print a message
# Returns true (=0) if it does not exists, and false (=1) if it does
function test_if_needed {
  if [ -f "$1" ]
  then
    echo "Not generating $1, already exists"
    return 1
  else
    echo "Generating $1"
    return 0
  fi
}

grep "^${WHICH_POLY} " ../c120.list | while read IDX N
do
  # Array with the file names of the .keep files, which are used several times
  unset KEEPFILENAMES
  declare -a KEEPFILENAMES

  # For each P, generate one file with size-optimized polynomials
  for P in `seq "${MINP}" "${STEPP}" "${MAXP}"`
  do
    OUTPUTFILE="c120${IDX}.polyselect1.${P}.0-${ADMAX}"
    if test_if_needed "${OUTPUTFILE}"
    then
      "${POLYSELECT2L}" -P "${P}" -N "${N}" -degree 5 -r -t 2 -admin 0 -admax "${ADMAX}" -incr 60 -nq 1000 -keep "${KEEP}" > "${OUTPUTFILE}"
    fi

    # Generate file of the best KEEP polynomials, sorted by increasing log norm
    KEEPFILE="c120${IDX}.polyselect1.${P}.0-${ADMAX}.keep_${KEEP}"
    KEEPFILENAMES+=("${KEEPFILE}") # bash 3.1+ syntax
    if test_if_needed "${KEEPFILE}"
    then
      ../parse_poly.py -keep "${KEEP}" "${OUTPUTFILE}" > "${KEEPFILE}"
    fi
  done

  # Generate a file containing the smallest lognorm from each output file, i.e.,
  # the smallest lognorm produced by each P value, one per line
  BESTLOGNORMFILE="c120${IDX}.best_lognorm"
  if test_if_needed "${BESTLOGNORMFILE}"
  then
    for KEEPFILE in "${KEEPFILENAMES[@]}"
    do
      grep "# lognorm" "${KEEPFILE}" | head -n 1 || exit 1
    done | sed "s/# lognorm //" >| "${BESTLOGNORMFILE}"
  fi

  # Like before, but takes the 10-th best lognorm from each P value, which is
  # less jittery
  BEST10LOGNORMFILE="c120${IDX}.10th_best_lognorm"
  if test_if_needed "${BEST10LOGNORMFILE}"
  then
    for KEEPFILE in "${KEEPFILENAMES[@]}"
    do
      grep "# lognorm" "$KEEPFILE" | head -n 10 | tail -n 1
    done | sed "s/# lognorm //" >| "${BEST10LOGNORMFILE}"
  fi

  # Like before, but takes the average of the 20-th best lognorm from each P value,
  # which is even less jittery
  AVGLOGNORMFILE="c120${IDX}.avg_lognorm"
  if test_if_needed "${AVGLOGNORMFILE}"
  then
    for KEEPFILE in "${KEEPFILENAMES[@]}"
    do
      echo -n "( "
      grep "# lognorm" "$KEEPFILE" | head -n 20 | sed "s/# lognorm //" | tr "\n" "+"
      echo "0 ) / 10"
    done | bc -l > "${AVGLOGNORMFILE}"
  fi

  # Generate a file containing all unique polynomials from all the KEEP files
  UNIQFILE="c120${IDX}.polyselect1.all.0-${ADMAX}.keep_${KEEP}.uniq"
  if test_if_needed "${UNIQFILE}"
  then
    ../parse_poly.py -uniq "${KEEPFILENAMES[@]}" > "${UNIQFILE}"
  fi

  ROPTFILE="${UNIQFILE}.ropt"
  if [ -f "${ROPTFILE}" ]
  then
    NEWFILE="${UNIQFILE}.new"
    echo "Generating ${NEWFILE} with polynomials in ${UNIQFILE} but not in ${ROPTFILE}"
    if [ -f "${NEWFILE}" ]
    then
      echo "Deleting existing ${NEWFILE}"
      rm -f "${NEWFILE}"
    fi
    ../parse_poly.py -nomatch "${ROPTFILE}" "${UNIQFILE}" > "${NEWFILE}"
  else
    # If no ROPTFILE exists, take all unique polynomials as input
    NEWFILE="${UNIQFILE}"
  fi

  if [ -s "${NEWFILE}" ] # if the NEWFILE is non-empty
  then
    echo "Root-optimizing new polynomials from ${NEWFILE} and appending to ${ROPTFILE}"
    "${POLYSELECT2L}" -rootsieve "${NEWFILE}" -t 2 -keep 10000 -area 67108864000000.0 -Bf 8000000.0 -Bg 4000000.0 >> "${ROPTFILE}"
  fi
done
