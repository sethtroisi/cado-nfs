
include_directories(${CADO_NFS_SOURCE_DIR}/numbertheory)

add_executable(test-ant test-ant.cpp)
target_link_libraries(test-ant ant)

cado_define_test(
    p-maximal-order
    PROGRAM
    ${CMAKE_CURRENT_BINARY_DIR}/test-ant
    -test p-maximal-order-batch -batch
    ${CMAKE_CURRENT_SOURCE_DIR}/p_maximal_order_quick.data
    TARGET_DEPENDENCIES
    test-ant
    )
