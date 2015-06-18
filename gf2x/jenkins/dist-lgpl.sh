. "`dirname $0`"/common.sh
autoreconf -i
src="$PWD"
TMP=`mktemp -d /tmp/${BUILD_TAG}-XXXXXXX`
if ! (cd "$TMP" ; $src/configure $configure_extra && make untaint && make dist) ; then
   echo "FAILED"
   rm -rf "$TMP"
   exit 1
fi
version=$(grep AC_INIT configure.ac | perl -ne '/(\d+(?:\.\d+)+)/ && print "$1\n";')
cd "$TMP"
tar xzf gf2x-$version.tar.gz
cd gf2x-$version
if ! (./configure $configure_extra && make && make check) ; then
   echo "FAILED"
   cd "$src"
   rm -rf "$TMP"
   exit 1
fi
cd "$src"
rm -rf "$TMP"
