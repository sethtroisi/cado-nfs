add_executable(test_array test_array.c ../../nfs-hd/sieving_interval.c
  ../../nfs-hd/int64_vector.c ../../nfs-hd/array.c)
FIND_LIBRARY(ECM_LIB ecm ${CADO_NFS_SOURCE_DIR}/nfs-hd)
target_link_libraries(test_array ${ECM_LIB} utils)
include_directories(${CADO_NFS_SOURCE_DIR}/nfs-hd)
run_test(test_array test_array)

# -----
add_executable(test_int64_poly test_int64_poly.c ../../nfs-hd/int64_poly.c
	../../nfs-hd/sieving_interval.c ../../nfs-hd/utils_int64.c
	../../utils/mpz_poly.c)
target_link_libraries(test_int64_poly ${ECM_LIBS} utils)
include_directories(${CADO_NFS_SOURCE_DIR}/nfs-hd)
run_test(test_int64_poly test_int64_poly)

# -----
add_executable(test_ideal test_ideal.c ../../nfs-hd/ideal.c
	../../utils/mpz_poly.c)
target_link_libraries(test_ideal ${ECM_LIBS} utils)
include_directories(${CADO_NFS_SOURCE_DIR}/nfs-hd)
run_test(test_ideal test_ideal)

# -----
add_executable(test_makefb test_makefb.c ../../nfs-hd/makefb.c
	../../utils/mpz_poly.c ../../nfs-hd/ideal.c ../../nfs-hd/utils_int64.c
  ../../nfs-hd/factor_base.c)
target_link_libraries(test_makefb ${ECM_LIBS} utils)
include_directories(${CADO_NFS_SOURCE_DIR}/nfs-hd)
run_test(test_makefb test_makefb)
