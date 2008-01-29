PROJECT(CADO_NFS)

## You can force a path to gmp.h using -DGMP_INCLUDE_PATH=/path/to/gmph
FIND_PATH(GMP_INCLUDE_PATH gmp.h
    /usr/include
    /usr/local/include)

## You can force the gmp lib using -DGMP_LIB=/opt/lib/libgmp.a
FIND_LIBRARY(GMP_LIB gmp
    /usr/lib
    /usr/lib64
    /usr/local/lib
    /usr/local/lib64)

## You can override these defaults with -DCMAKE_C_FLAGS="-g -fomit-frame-pointer"
IF (CMAKE_C_FLAGS STREQUAL " ")
    SET (CMAKE_C_FLAGS "-Wall -W -O2")
ENDIF (CMAKE_C_FLAGS STREQUAL " ")

INCLUDE_DIRECTORIES (${CADO_NFS_SOURCE_DIR} ${GMP_INCLUDE_PATH}
    ${CADO_NFS_SOURCE_DIR}/utils)
LINK_DIRECTORIES (${CADO_NFS_BINARY_DIR}/utils)

ADD_SUBDIRECTORY (utils)
ADD_SUBDIRECTORY (polyselect)
ADD_SUBDIRECTORY (sieve)
ADD_SUBDIRECTORY (postsieve/checknorms)
ADD_SUBDIRECTORY (postsieve/tifa)
ADD_SUBDIRECTORY (linalg)
ADD_SUBDIRECTORY (sqrt/naive)
