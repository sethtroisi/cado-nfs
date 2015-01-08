add_executable(test_array test_array.c ../../nfs-hd/sieving_interval.c
  ../../nfs-hd/int64_vector.c ../../nfs-hd/array.c)
FIND_LIBRARY(ECM_LIB ecm ${CADO_NFS_SOURCE_DIR}/nfs-hd)
target_link_libraries(test_array ${ECM_LIB} utils)
include_directories(${CADO_NFS_SOURCE_DIR}/nfs-hd)
run_test(test_array test_array)
