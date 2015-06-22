. "`dirname $0`"/common.sh
autoreconf -i
./configure $configure_extra
make
make check
