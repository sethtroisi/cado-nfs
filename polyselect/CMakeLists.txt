SET (sources 
    polyselect.c)

## Add -std=c99 and Revision stuff to compile flags for this target
SET (CFLAGS "-std=c99 -DREV=\\\"`svnversion ${CMAKE_CURRENT_SOURCE_DIR}`\\\"")

SET_SOURCE_FILES_PROPERTIES (${sources}
    COMPILE_FLAGS ${CFLAGS})

ADD_EXECUTABLE (polyselect ${sources})
TARGET_LINK_LIBRARIES (polyselect utils ${GMP_LIB} m)

