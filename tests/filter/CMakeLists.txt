# Build sm before running check

# Get path to the sm binary and store in variable sm_path
# see comment in tests/sieve/CMakeLists.txt
set(sm_path ${CADO_NFS_BINARY_DIR}/filter/sm)

add_test(testsm bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/testsm.sh\" \"${sm_path}\"")
set_test_build_depends(testsm sm)
