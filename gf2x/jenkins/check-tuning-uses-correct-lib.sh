#!/usr/bin/env bash
. "`dirname $0`"/common.sh
wdir="`mktemp -d /tmp/XXXXXXXXXX`"
cleanup() { rm -rf "$wdir" ; }
autoreconf -i
if [ "$DEBUG_SCRIPTS" ] ; then
    echo "data will be left in $wdir"
    set -x
else
    trap cleanup EXIT
fi

set -e

# we need a tarball. build it out of source while we're at it...
src="$PWD"
mkdir $wdir/prepare
cd $wdir/prepare
$src/configure $configure_extra
make dist
eval `egrep '^PACKAGE_(TARNAME|VERSION)=' $src/configure`
cd $wdir
tar xzf "$wdir/prepare/$PACKAGE_TARNAME-$PACKAGE_VERSION".tar.gz
cd "$PACKAGE_TARNAME-$PACKAGE_VERSION"

# first do a bastardized build. install it.
# some ex binaries insist on having a TERM variable.
TERM=xterm ex toom.c <<EOF
/void gf2x_mul_kara
/{
a
fprintf(stderr, "Hello, world\n");
abort();
.
wq
EOF
$src/configure --prefix=$wdir/inst $configure_extra
make
make install

# now rebuild the good source, and try to see whether the included
# binaries pass. Make sure we rebuild the whole thing !
cd $wdir
tar xzf "$wdir/prepare/$PACKAGE_TARNAME-$PACKAGE_VERSION".tar.gz
cd "$PACKAGE_TARNAME-$PACKAGE_VERSION"
make distclean
$src/configure --prefix=$wdir/inst $configure_extra
make
# it's a bit unfortunate as far as check times are concerned, but we need
# to do testing at least that far in order to have a valid threshold at
# the end of the day.
make tune-toom TOOM_TUNING_LIMIT=64

if [ "$DEBUG_SCRIPTS" ] ; then
    echo "data left in $wdir"
fi
