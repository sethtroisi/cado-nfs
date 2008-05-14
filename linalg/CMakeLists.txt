SET (sources 
    hashpair.c
    freerel.c
    duplicates.c
    prune.c
    purge.c
    merge.c
    replay.c
    files.c
    linalg.c
    gauss.c
    readmat.c
    transpose.c
    characters.c
    sparse.c
    allsqrt.c)

## Add Revision stuff to compile flags for this target
SET (CFLAGS "-DREV=\\\"`svnversion ${CMAKE_CURRENT_SOURCE_DIR}`\\\"")

SET_SOURCE_FILES_PROPERTIES (${sources}
    COMPILE_FLAGS ${CFLAGS})

SET_SOURCE_FILES_PROPERTIES (gauss.c
    COMPILE_FLAGS "-DMULTI_ROW=3 -DNO_MAIN")

SET (LIBS utils ${GMP_LIB} m)

ADD_EXECUTABLE(allsqrt allsqrt.c files.c hashpair.c)
TARGET_LINK_LIBRARIES(allsqrt ${LIBS})

ADD_EXECUTABLE(replay replay.c sparse.c)
TARGET_LINK_LIBRARIES(replay ${LIBS})

ADD_EXECUTABLE(merge merge.c sparse.c prune.c)
TARGET_LINK_LIBRARIES(merge ${LIBS})

ADD_EXECUTABLE(purge purge.c files.c hashpair.c)
TARGET_LINK_LIBRARIES(purge ${LIBS})

ADD_EXECUTABLE(freerel freerel.c ../sieve/fb.c hashpair.c)
TARGET_LINK_LIBRARIES(freerel ${LIBS})

ADD_EXECUTABLE(linalg linalg.c gauss.c readmat.c)
TARGET_LINK_LIBRARIES(linalg ${LIBS})

ADD_EXECUTABLE(characters characters.c files.c gauss.c)
TARGET_LINK_LIBRARIES(characters ${LIBS})

ADD_EXECUTABLE(transpose transpose.c readmat.c)
TARGET_LINK_LIBRARIES(transpose ${LIBS})
