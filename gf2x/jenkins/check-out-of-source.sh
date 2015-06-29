. "`dirname $0`"/common.sh
autoreconf -i
src="$PWD"
TMP=`mktemp -d /tmp/${BUILD_TAG}-XXXXXXX`
cd "$TMP"
if ! ($src/configure $configure_extra && make && make check) ; then
   echo "FAILED"
   cd "$src"
   rm -rf "$TMP"
   exit 1
fi
cd "$src"
rm -rf "$TMP"
