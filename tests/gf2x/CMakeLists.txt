run_test(gf2x ${CMAKE_CURRENT_SOURCE_DIR}/test_gf2x.sh ${CADO_NFS_BINARY_DIR})
set_test_build_depends(gf2x gf2x-build)
link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft)
include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/fft/mpfq)
add_executable(test_cantor test_cantor.c)
# This **MUST** match the #define in gf2x/Makefile.am
set_source_files_properties(test_cantor.c COMPILE_FLAGS
        "-DCANTOR_BASE_FIELD_SIZE=128")
target_link_libraries(test_cantor gf2x tests ${gmp_libname})
run_target_as_a_test(test_cantor 100000)
