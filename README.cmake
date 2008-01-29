In order to compile CADO-NFS with cmake:

1) get cmake installed on your computer.
2) create a build directory somewhere:
     mkdir /tmp/build
3) cd to this build dir and type:
     cd /tmp/build
     cmake /path/to/cadonfs
     make

You will get in the /tmp/build directory a directory tree with binaries
in appropriate positions. At some point, we should probably move them
automatically to a bin subdir, I guess.

WARNING: Don't build with cmake inside the source tree, since there is no
equivalent to "make distclean" that destroys files created by cmake. Your
source tree will be a mess if you do "cmake ." .

For passing flags to the compiler:
     cmake /path/to/cadonfs -DCMAKE_C_FLAGS="-O4 -pedantic -sse7"

