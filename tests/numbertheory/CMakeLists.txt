
include_directories(${CADO_NFS_SOURCE_DIR}/numbertheory)

add_executable(test-ant test-ant.cpp)
target_link_libraries(test-ant ant)

cado_define_test(
    numbertheory-p-maximal-order
    PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test-ant
    -test p-maximal-order-batch -batch
    ${CMAKE_CURRENT_SOURCE_DIR}/p_maximal_order.data
    TARGET_DEPENDENCIES
    test-ant
    )
cado_define_test(
    numbertheory-factorization-of-prime
    PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test-ant
    -test p-maximal-order-batch -batch
    ${CMAKE_CURRENT_SOURCE_DIR}/factorization_of_prime.data
    TARGET_DEPENDENCIES
    test-ant
    )

find_program(MAGMA magma HINTS ENV MAGMA)
if(MAGMA AND CHECKS_EXPENSIVE)
cado_define_test(
    numbertheory-p-maximal-order-expensive
    PROGRAM
    bash -c "${CMAKE_CURRENT_BINARY_DIR}/test-ant -test p-maximal-order-batch -batch <(${MAGMA} -b -S 1 ${CMAKE_CURRENT_SOURCE_DIR}/gentest_p_maximal_order.m < /dev/null)"
    TARGET_DEPENDENCIES
    test-ant
    )
cado_define_test(
    numbertheory-factorization-of-prime-expensive
    PROGRAM
    bash -c "${CMAKE_CURRENT_BINARY_DIR}/test-ant -test p-maximal-order-batch -batch <(${MAGMA} -b -S 1 ${CMAKE_CURRENT_SOURCE_DIR}/gentest_factorization_of_prime.m < /dev/null)"
    TARGET_DEPENDENCIES
    test-ant
    )
endif()
