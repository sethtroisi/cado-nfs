
include_directories(${CADO_NFS_SOURCE_DIR}/linalg)

cado_define_test(test_gauss.c ${CADO_NFS_SOURCE_DIR}/linalg/gauss.c
    LIBRARIES ${gmp_libname})


add_subdirectory (bwc)
