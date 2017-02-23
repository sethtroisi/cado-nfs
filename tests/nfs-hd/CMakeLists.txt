include_directories(${CADO_NFS_SOURCE_DIR})
include_directories(${CADO_NFS_SOURCE_DIR}/nfs-hd)

# -----

cado_define_test(test_array.c ${CADO_NFS_SOURCE_DIR}/nfs-hd/sieving_bound.c
	${CADO_NFS_SOURCE_DIR}/nfs-hd/int64_vector.c
	${CADO_NFS_SOURCE_DIR}/nfs-hd/array.c LIBRARIES utils)

# -----
cado_define_test(test_int64_poly.c ${CADO_NFS_SOURCE_DIR}/nfs-hd/int64_poly.c
	${CADO_NFS_SOURCE_DIR}/nfs-hd/sieving_bound.c
	${CADO_NFS_SOURCE_DIR}/nfs-hd/utils_int64.c
	${CADO_NFS_SOURCE_DIR}/utils/mpz_poly.cpp LIBRARIES utils)

# -----
cado_define_test(test_ideal.c ${CADO_NFS_SOURCE_DIR}/nfs-hd/ideal.c
	${CADO_NFS_SOURCE_DIR}/utils/mpz_poly.cpp LIBRARIES utils)

# -----
cado_define_test(test_utils_int64.c ${CADO_NFS_SOURCE_DIR}/nfs-hd/utils_int64.c
	LIBRARIES utils)

# -----
cado_define_test(test_mat_int64.c
	${CADO_NFS_SOURCE_DIR}/nfs-hd/mat_int64.c
	${CADO_NFS_SOURCE_DIR}/nfs-hd/utils_int64.c
	${CADO_NFS_SOURCE_DIR}/nfs-hd/mat_Z.c
	${CADO_NFS_SOURCE_DIR}/nfs-hd/int64_poly.c
	${CADO_NFS_SOURCE_DIR}/nfs-hd/int64_vector.c LIBRARIES utils)

# -----
cado_define_test(test_alpha3d.c
	${CADO_NFS_SOURCE_DIR}/nfs-hd/alpha3d.c
	LIBRARIES utils)

# -----
#cado_define_test(test_resultant.c ${CADO_NFS_SOURCE_DIR}/nfs-hd/sieving_bound.c
	#${CADO_NFS_SOURCE_DIR}/utils/mpz_poly.cpp
	#${CADO_NFS_SOURCE_DIR}/nfs-hd/int64_vector.c
	#${CADO_NFS_SOURCE_DIR}/nfs-hd/mat_int64.c
	#${CADO_NFS_SOURCE_DIR}/nfs-hd/int64_poly.c
	#${CADO_NFS_SOURCE_DIR}/nfs-hd/mat_Z.c
	#${CADO_NFS_SOURCE_DIR}/utils/double_poly.c LIBRARIES utils)
