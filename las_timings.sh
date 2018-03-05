#!/bin/sh
USAGE="Usage: ./las_timings.sh POLYFILE PARAMFILE"

if [ -n "$1" ] ; then
  POLY=$1
else
  echo "Error, missing poly file"
  echo ${USAGE}
  exit 1
fi
if [ -n "$2" ] ; then
  PARAM=$2
else
  echo "Error, missing param file"
  echo ${USAGE}
  exit 1
fi

filename()
{
  echo $(git show -s --format='%ci' $1 \
          | awk '//{ print $1, $2;}' \
          | tr -d -c 0-9)-$(git show -s --format='%h' $1)
}

OUTDIR="las_timings"
mkdir -p ${OUTDIR} 
echo "# outdir = ${OUTDIR}"

NQ=25 # TODO make it a command line parameter

I=$(grep tasks.I ${PARAM} | tr -d " " | cut -d= -f 2)
LIM0=$(grep lim0 ${PARAM} | tr -d " " | cut -d= -f 2)
LIM1=$(grep lim1 ${PARAM} | tr -d " " | cut -d= -f 2)
LPB0=$(grep lpb0 ${PARAM} | tr -d " " | cut -d= -f 2)
LPB1=$(grep lpb1 ${PARAM} | tr -d " " | cut -d= -f 2)
MFB0=$(grep mfb0 ${PARAM} | tr -d " " | cut -d= -f 2)
MFB1=$(grep mfb1 ${PARAM} | tr -d " " | cut -d= -f 2)

#-powlim1 32767
#-powlim0 32767
echo "# param: I = ${I}"
echo "# param: lim0 = ${LIM0}"
echo "# param: lim1 = ${LIM1}"
echo "# param: lpb0 = ${LPB0}"
echo "# param: lpb1 = ${LPB1}"
echo "# param: mfb0 = ${MFB0}"
echo "# param: mfb1 = ${MFB1}"

SQSIDE=1 # Assume sqside = 1
Q0=${LIM1}
T=1
H="${I} ${LIM0} ${LIM1} ${LPB0} ${LPB1} ${MFB0} ${MFB1} ${SQSIDE} ${Q0} ${NQ} ${T}"

POLYSTR=$(cat ${POLY})
prefix_fb=$(echo ${POLYSTR} ${I} ${LIM1} | md5sum | cut -c1-8)
prefix=$(filename HEAD)-$(echo ${POLYSTR} ${H} | md5sum | cut -c1-8)
echo "# prefix_fb = ${prefix_fb}"
echo "# prefix = ${prefix}"

ROOTSFILE="${OUTDIR}/${prefix_fb}.roots"
LOG="${OUTDIR}/${prefix}.log"

##### With new default strategy
cp --backup=numbered local.sh local.sh.bak 2> /dev/null
echo 'CFLAGS="-O3 -funroll-loops -DNDEBUG"' > local.sh
echo 'CXXFLAGS="-O3 -funroll-loops -DNDEBUG"' >> local.sh
echo 'build_tree="${up_path}/build/cofact.twed12"' >> local.sh
BUILDDIR="./build/cofact.twed12"
RELSFILE=${OUTDIR}/${prefix}.twed12.rels

echo "# twed12: make cmake"
make cmake > ${LOG} 2>&1
echo "# twed12: make makefb las"
make -j2 makefb las >> ${LOG} 2>&1

if [ -e "${ROOTSFILE}" ] ; then
  echo "# twed12: skipping makefb"
else
  echo "# twed12: doing makefb...."
  ${BUILDDIR}/sieve/makefb -lim ${LIM1} -maxbits ${I} -poly ${POLY} -t 2 > ${ROOTSFILE}
  echo "# twed12: done"
fi

echo "# twed12: doing las..."
${BUILDDIR}/sieve/las -poly ${POLY} -sqside ${SQSIDE} -lim0 ${LIM0} \
                      -lim1 ${LIM1} -lpb0 ${LPB0} -lpb1 ${LPB1} -mfb0 ${MFB0} \
                      -mfb1 ${MFB1} -fb1 ${ROOTSFILE} -I ${I} -q0 ${Q0} \
                      -nq ${NQ} -t ${T} 2>&1 > ${RELSFILE}
echo "# twed12: done"

TWED_COFACT_TIME=$(grep "^# Total cpu time" ${RELSFILE} \
                                    | sed 's/^.*+\s\([0-9]*\.[0-9]*\))]$/\1/')
echo "# twed12: cofact time: ${TWED_COFACT_TIME}"

###### With old strategy
echo 'CFLAGS="-O3 -funroll-loops -DNDEBUG -DUSE_LEGACY_DEFAULT_STRATEGY=1"' > local.sh
echo 'CXXFLAGS="-O3 -funroll-loops -DNDEBUG -DUSE_LEGACY_DEFAULT_STRATEGY=1"' >> local.sh
echo 'build_tree="${up_path}/build/cofact.monty"' >> local.sh
BUILDDIR="./build/cofact.monty"
RELSFILE=${OUTDIR}/${prefix}.monty.rels

echo "# monty: make cmake"
make cmake >> ${LOG} 2>&1
echo "# monty: make makefb las"
make -j2 makefb las >> ${LOG} 2>&1

if [ -e "${ROOTSFILE}" ] ; then
  echo "# monty: skipping makefb"
else
  echo "# monty: doing makefb...."
  ${BUILDDIR}/sieve/makefb -lim ${LIM1} -maxbits ${I} -poly ${POLY} -t 2 > ${ROOTSFILE}
  echo "# monty: done"
fi

echo "# monty: doing las..."
${BUILDDIR}/sieve/las -poly ${POLY} -sqside ${SQSIDE} -lim0 ${LIM0} \
                      -lim1 ${LIM1} -lpb0 ${LPB0} -lpb1 ${LPB1} -mfb0 ${MFB0} \
                      -mfb1 ${MFB1} -fb1 ${ROOTSFILE} -I ${I} -q0 ${Q0} \
                      -nq ${NQ} -t ${T} 2>&1 > ${RELSFILE}
echo "# monty: done"

MONTY_COFACT_TIME=$(grep "^# Total cpu time" ${RELSFILE} \
                                    | sed 's/^.*+\s\([0-9]*\.[0-9]*\))]$/\1/')
echo "# monty: cofact time: ${MONTY_COFACT_TIME}"

echo "### twed12: cofact time: ${TWED_COFACT_TIME}" >> ${LOG}
echo "### monty: cofact time: ${MONTY_COFACT_TIME}" >> ${LOG}

##### restore local.sh
rm local.sh
mv local.sh.bak local.sh 2> /dev/null
