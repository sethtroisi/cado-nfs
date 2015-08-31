# Build polyselect2l before running check

cado_define_test(test_polsel SCRIPT
    ${CMAKE_CURRENT_SOURCE_DIR}/polysel.sh
    NO_DEFAULT_RUN
    TARGET_DEPENDENCIES polyselect2l)

cado_divert_test(test_polsel 2 ${CADO_NFS_BINARY_DIR} 2 1000 100 8)
cado_divert_test(test_polsel 3 ${CADO_NFS_BINARY_DIR} 3 1000 500 10)
cado_divert_test(test_polsel 4 ${CADO_NFS_BINARY_DIR} 4 1000 1000 10)
cado_divert_test(test_polsel 5 ${CADO_NFS_BINARY_DIR} 5 1000 1000 10)
cado_divert_test(test_polsel 6 ${CADO_NFS_BINARY_DIR} 6 1000 1000 10)
cado_divert_test(test_polsel 7 ${CADO_NFS_BINARY_DIR} 7 1000 1000 10)

cado_define_test(test_auxiliary.c LIBRARIES utils tests polyselect_common)
cado_define_test(test_rho.c LIBRARIES polyselect_common tests)
