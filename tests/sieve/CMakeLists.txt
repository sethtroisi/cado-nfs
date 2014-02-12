# CMake 2.8.4 and newer has a WORKING_DIRECTORY property, but we may want
# to stay compatible with older versions
add_test(F9_sievetest "${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh" "${CADO_NFS_SOURCE_DIR}" "${CADO_NFS_BINARY_DIR}")
add_executable(test_trialdiv test_trialdiv.c "${CADO_NFS_SOURCE_DIR}/sieve/trialdiv.c")
target_link_libraries(test_trialdiv utils tests)
run_test(test_trialdiv test_trialdiv)

add_subdirectory (ecm)
