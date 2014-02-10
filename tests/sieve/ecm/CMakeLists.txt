# Tests for the P-1 algorithm.
add_test(test_pm1 "${CMAKE_CURRENT_SOURCE_DIR}/test_testbench.sh" "${CADO_NFS_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/pm1test.inp" -q -vf -vcf -pm1 100 1000)
