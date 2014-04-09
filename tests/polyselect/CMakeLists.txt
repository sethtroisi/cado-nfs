# Build polyselect2l before running check
add_dependencies(check polyselect2l)

add_test(test_polsel2 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/polysel.sh\" \"${CADO_NFS_BINARY_DIR}\" 2 1000 100 8")
add_test(test_polsel3 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/polysel.sh\" \"${CADO_NFS_BINARY_DIR}\" 3 1000 500 10")
add_test(test_polsel4 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/polysel.sh\" \"${CADO_NFS_BINARY_DIR}\" 4 1000 1000 10")
add_test(test_polsel5 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/polysel.sh\" \"${CADO_NFS_BINARY_DIR}\" 5 1000 1000 10")
add_test(test_polsel6 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/polysel.sh\" \"${CADO_NFS_BINARY_DIR}\" 6 1000 1000 10")
add_test(test_polsel7 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/polysel.sh\" \"${CADO_NFS_BINARY_DIR}\" 7 1000 1000 10")
add_executable(test_auxiliary test_auxiliary.c)
target_link_libraries(test_auxiliary utils tests polyselect_common)
run_test(test_auxiliary test_auxiliary)
build_and_run_test(test_rho.c "polyselect_common;tests")
