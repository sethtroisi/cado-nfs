
find_program(MAGMA magma HINTS ENV MAGMA)

if(MAGMA)
    cado_define_test(
        SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test-numbertheory.sh
        --binary ${CADO_NFS_BINARY_DIR}/numbertheory/ant
        --magma ${MAGMA}
        --fields 4
        TARGET_DEPENDENCIES
        ant
        )
endif()
