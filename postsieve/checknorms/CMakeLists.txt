SET (sources 
    checknorms.c
    checknorms.2lp.c)

## Add -std=c99 and Revision stuff to compile flags for this target
SET (CFLAGS "-std=c99 -fkeep-inline-functions -DREV=\\\"`svnversion ${CMAKE_CURRENT_SOURCE_DIR}`\\\"")

SET_SOURCE_FILES_PROPERTIES (${sources}
    COMPILE_FLAGS ${CFLAGS})

LINK_DIRECTORIES (${CADO_NFS_BINARY_DIR}/postsieve/tifa)
SET (TIFADIR ${CADO_NFS_SOURCE_DIR}/postsieve/tifa)
INCLUDE_DIRECTORIES (${TIFADIR} ${TIFADIR}/lib
    ${TIFADIR}/lib/algo/include ${TIFADIR}/lib/data/include
    ${TIFADIR}/lib/utils/include)

ADD_EXECUTABLE (checknorms checknorms.c)
TARGET_LINK_LIBRARIES (checknorms tifa utils ${GMP_LIB} m)

ADD_EXECUTABLE (checknorms.2lp checknorms.2lp.c)
TARGET_LINK_LIBRARIES (checknorms.2lp tifa utils ${GMP_LIB} m)
