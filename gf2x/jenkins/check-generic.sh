. "`dirname $0`"/common.sh
autoreconf -i
./configure $configure_extra
make
make check
if [ -f "`dirname $0`"/extra-"`basename $0`" ] ; then
    . "`dirname $0`"/extra-"`basename $0`"
fi

# do that on all architectures.
make distclean
./configure $configure_extra --disable-hardware-specific-code
make
make check
make tune-lowlevel
make tune-toom TOOM_TUNING_LIMIT=64
if [ -f "`dirname $0`"/extra-"`basename $0`" ] ; then
    . "`dirname $0`"/extra-"`basename $0`"
fi

