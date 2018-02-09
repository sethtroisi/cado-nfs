#!/bin/sh

# 32-56 bits
/tmp/cado-nfs-build/sieve/ecm/testbench -p -cof 4294967311 -ecm 315 5355 11 3 10000000

/tmp/cado-nfs-build/sieve/ecm/testbench -p -cof 4294967311 -ecmem12 315 5355 1 3 10000000

# 64-88 bits
/tmp/cado-nfs-build/sieve/ecm/testbench -p -cof 18446744073709551629 -ecm 315 5355 11 3 10000000

/tmp/cado-nfs-build/sieve/ecm/testbench -p -cof 18446744073709551629 -ecmem12 315 5355 1 3 10000000

# 96-120 bits
/tmp/cado-nfs-build/sieve/ecm/testbench -p -cof 79228162514264337593543950397 -ecm 315 5355 11 3 10000000

/tmp/cado-nfs-build/sieve/ecm/testbench -p -cof 79228162514264337593543950397 -ecmem12 315 5355 1 3 10000000

# 128-152 bits
/tmp/cado-nfs-build/sieve/ecm/testbench -p -cof 340282366920938463463374607431768211507 -ecm 315 5355 11 3 10000000

/tmp/cado-nfs-build/sieve/ecm/testbench -p -cof 340282366920938463463374607431768211507 -ecmem12 315 5355 1 3 10000000

# 196-220 bits
/tmp/cado-nfs-build/sieve/ecm/testbench -p -cof 100433627766186892221372630771322662657637687111424552206357 -ecm 315 5355 11 3 10000000

/tmp/cado-nfs-build/sieve/ecm/testbench -p -cof 100433627766186892221372630771322662657637687111424552206357 -ecmem12 315 5355 1 3 10000000


