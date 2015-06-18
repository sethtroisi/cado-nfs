. "`dirname $0`"/common.sh
autoreconf -i
./configure $configure_extra
make
make tune-lowlevel
make tune-toom
