#!/usr/bin/env bash

# The main functions to compute SM are tested in utils/test_sm_utils.
# This test is here to check that the multi-threaded and mono-threaded versions
# give the same results and that the -nsm option works correctly.

SRCDIR="$1"
BINDIR="$2"
SOURCE_TEST_DIR="`dirname "$0"`"
TMPDIR=`mktemp -d /tmp/cadotest.XXXXXXXX`
# Make temp direcotry world-readable for easier debugging
chmod a+rx "${TMPDIR}"

poly="${SOURCE_TEST_DIR}/testsm.p59.poly"
purged="${SOURCE_TEST_DIR}/testsm.p59.purged.gz"
id="${SOURCE_TEST_DIR}/testsm.p59.index.gz"
go="2926718140519"
smexp="8565679074042993029589360"

SM="${BINDIR}/filter/sm"
cd ${SRCDIR}
make sm

args="-poly ${poly} -purged ${purged} -index ${id} -gorder ${go} -smexp ${smexp}"

#without -nsm (ie nsm = deg F) and -mt 1
${SM} ${args} -out ${TMPDIR}/sm.5.1 -mt 1
if [ "$?" -ne "0" ] ; then
  echo "$0: sm binary failed with -mt 1 and without -nsm. Files remain in ${TMPDIR}"
  exit 1
fi
#without -nsm (ie nsm = deg F) and -mt 2
${SM} ${args} -out ${TMPDIR}/sm.5.2 -mt 2
if [ "$?" -ne "0" ] ; then
  echo "$0: sm binary failed with -mt 2 and without -nsm. Files remain in ${TMPDIR}"
  exit 1
fi

#with -nsm 2 and -mt 1
${SM} ${args} -out ${TMPDIR}/sm.2.1 -mt 1 -nsm 2
if [ "$?" -ne "0" ] ; then
  echo "$0: sm binary failed with -mt 1 and -nsm 2. Files remain in ${TMPDIR}"
  exit 1
fi
#with -nsm 2 and -mt 2
${SM} ${args} -out ${TMPDIR}/sm.2.2 -mt 2 -nsm 2
if [ "$?" -ne "0" ] ; then
  echo "$0: sm binary failed with -mt 2 and -nsm 2. Files remain in ${TMPDIR}"
  exit 1
fi


diff -q ${TMPDIR}/sm.5.1 ${TMPDIR}/sm.5.2
if [ "$?" -ne "0" ] ; then
  echo "$0: Mono-threaded and multi-threaded versions do not match (without -nsm). Files remain in ${TMPDIR}"
  exit 1
fi

diff -q ${TMPDIR}/sm.2.1 ${TMPDIR}/sm.2.2
if [ "$?" -ne "0" ] ; then
  echo "$0: Mono-threaded and multi-threaded versions do not match (with -nsm 2). Files remain in ${TMPDIR}"
  exit 1
fi

cut -d " " -f 1-2 ${TMPDIR}/sm.5.1 > ${TMPDIR}/sm.5.1.short
diff -q ${TMPDIR}/sm.5.1.short ${TMPDIR}/sm.2.1
if [ "$?" -ne "0" ] ; then
  echo "$0: First two SMs computed without -nsm do not match SMs computed with -nsm 2). Files remain in ${TMPDIR}"
  exit 1
fi

rm -f ${TMPDIR}/sm.5.1 ${TMPDIR}/sm.5.2 ${TMPDIR}/sm.2.1 ${TMPDIR}/sm.2.2
rm -f ${TMPDIR}/sm.5.1.short
rmdir ${TMPDIR}
