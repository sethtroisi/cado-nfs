# CMake 2.8.4 and newer has a WORKING_DIRECTORY property, but we may want
# to stay compatible with older versions
add_test(F9_sievetest "${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh" "${CADO_NFS_SOURCE_DIR}" "${CADO_NFS_BINARY_DIR}")
