SET (sources 
    makefb.c
    sieve.c
    fb.c)

## Add specific compile flags for these targets
SET (CFLAGS "-fkeep-inline-functions -DREV=\\\"`svnversion ${CMAKE_CURRENT_SOURCE_DIR}`\\\"")

SET_SOURCE_FILES_PROPERTIES (${sources}
    COMPILE_FLAGS ${CFLAGS})

ADD_EXECUTABLE (sieve sieve.c fb.c)
TARGET_LINK_LIBRARIES (sieve utils ${GMP_LIB} m)

ADD_EXECUTABLE (makefb makefb.c)
TARGET_LINK_LIBRARIES (makefb utils ${GMP_LIB} m)

