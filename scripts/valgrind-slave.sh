# This is a snapshot of the shell script which is run on the continuous
# integration platform to do valgrind checking of the cado-nfs code base.

# This is only informational. The authoritative script is the one on the
# jenkins master (and *not* this one!)

# Move old test results to Testing-old/ to prevent the xUnit plug-in from finding old result files and complaining
#BUILD=`ls -d build/*/`
#mkdir -p "$BUILD"/Testing-old
#mv --backup=numbered "$BUILD"/Testing/20* "$BUILD"/Testing-old/ || true
#make check ARGS="-T Test"
CHECKS_EXPENSIVE=1 make cmake
CHECKS_EXPENSIVE=1 make -j4 full_c59_dependencies
export FOO=$(uuidgen)
mkdir valgrind.$BUILD_NUMBER.$FOO
# Don't use --error-exitcode, so that we get a chance to be notified of all potential errors at once.
VALGRIND="valgrind --suppressions=cado-nfs.supp --trace-children=yes --trace-children-skip=/bin/sh,gzip,*perl*,*python* --log-file=valgrind.%q{BUILD_NUMBER}.%q{FOO}/pid-%p --leak-check=full"
set +e
./factor.sh 90377629292003121684002147101760858109247336549001090677693 -t 2 tasks.runprefix="$VALGRIND" tasks.linalg.bwc.precmd="$VALGRIND"
rc=$?
set -e
(cd valgrind.$BUILD_NUMBER.$FOO ; mkdir ok nok ; grep -l 'ERROR SUMMARY: [^0]' pid-* | xargs -r mv -t nok ; ls | grep pid | xargs -r grep -l 'ERROR SUMMARY: 0' | xargs -r mv -t ok)
set +e
ls valgrind.${BUILD_NUMBER}.${FOO}/nok | grep -q .
found_nok_files=$?
set -e
ls valgrind.${BUILD_NUMBER}.${FOO}/nok | while read f ; do
  echo "Errors in file valgrind.${BUILD_NUMBER}.${FOO}/nok/$f"
  cat valgrind.${BUILD_NUMBER}.${FOO}/nok/$f
done
tar cvzf valgrind.${BUILD_NUMBER}.${FOO}.tar.gz valgrind.${BUILD_NUMBER}.${FOO}/
rm -rf valgrind.${BUILD_NUMBER}.${FOO}
if [ $rc != 0 ] ; then
 echo "Cadofactor exit code was $rc"
 exit $rc
fi
if [ $found_nok_files != 0 ] ; then
  echo "Found valgrind errors"
  echo "See archive of log files in `hostname`:$PWD/valgrind.${BUILD_NUMBER}.${FOO}.tar.gz"
  exit 1
fi
