cado_define_test(SCRIPT
    ${CMAKE_CURRENT_SOURCE_DIR}/test_gf2x.sh ${CADO_NFS_BINARY_DIR}
    TARGET_DEPENDENCIES gf2x-build)

link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft)
include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft/mpfq)

cado_define_test(test_cantor.c
    ARGUMENTS 100000
    LIBRARIES gf2x tests ${gmp_libname})
# This **MUST** match the #define in gf2x/Makefile.am
set_source_files_properties(test_cantor.c COMPILE_FLAGS
        "-DCANTOR_BASE_FIELD_SIZE=128")
