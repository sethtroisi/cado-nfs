
include_directories(${CADO_NFS_SOURCE_DIR}/numbertheory)

add_executable(test-ant test-ant.cpp)
target_link_libraries(test-ant ant)

find_program(MAGMA magma HINTS ENV MAGMA)

if(MAGMA)
    cado_define_test(
        SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test-p-maximal-order.sh
        --binary
        ${CMAKE_CURRENT_BINARY_DIR}/numbertheory/test-ant
        --magma ${MAGMA}
        --fields 4
        TARGET_DEPENDENCIES
        test-ant
        )
    cado_define_test(
        SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test-factorization-of-prime.sh
        --binary
        ${CMAKE_CURRENT_BINARY_DIR}/numbertheory/test-ant
        --magma ${MAGMA}
        --fields 4
        TARGET_DEPENDENCIES
        test-ant
        )
endif()
