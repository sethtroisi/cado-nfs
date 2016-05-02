set (SQ_SRC special-q_sieve.c mat_int64.c int64_vector.c uint64_array.c
    mat_Z.c ideal.c makefb.c utils_mpz.c int64_poly.c sieving_bound.c
    array.c utils_int64.c factor_base.c utils_lattice.c double_vector.c
    mat_double.c list_int64_vector.c list_double_vector.c
    list_int64_vector_index.c utils_norm.c utils_cofactorisation.c
    utils_mpz_poly.c)

add_executable(special-q_sieve ${SQ_SRC})

#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS TIME_SIEVES ENUM_LATTICE_INSTEAD_OF_SPACE_SIEVE)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS SKEW_LLL_SPQ SPQ_DEFINED)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS Q_BELOW_FBB)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS ASSERT_FACTO
#    MEAN_NORM_BOUND MEAN_NORM_BOUND_SIEVE TEST_MQR NUMBER_HIT)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS COME_BACK_FK
#	TEST_MQR SPECIAL_Q_IDEAL_U)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS
	#ASSERT_FACTO TEST_MQR SKEW_LLL)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS
	#ASSERT_FACTO TEST_MQR)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS
	#ASSERT_FACTO ASSERT_SIEVE TRACE_POS=1000)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS ASSERT_FACTO
#	NUMBER_HIT TEST_MQR)# TRACE_POS=0	ASSERT_SIEVE)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS ASSERT_FACTO
#	MEAN_NORM_BOUND MEAN_NORM_BOUND_SIEVE NUMBER_HIT NUMBER_SURVIVALS
#	PRINT_PARAMETERS TIMER_SIEVE=100 TEST_MQR TRACE_POS=0
#	MEAN_NORM ASSERT_SIEVE TRACE_POS=0 ENUM_LATTICE)

include_directories (${CMAKE_SOURCE_DIR}/sieve)
target_link_libraries(special-q_sieve facul utils)

# Makefb

set(MAKEFB_SRC makefb.c factor_base.c ideal.c utils_int64.c utils_mpz_poly.c)

add_executable(makefb-hd ${MAKEFB_SRC})

set_property(TARGET makefb-hd PROPERTY COMPILE_DEFINITIONS MAIN)

target_link_libraries(makefb-hd utils)

# Polynomial selection

#set (POL_SRC polyselect-hd.c polynomials.c sieving_bound.c utils_mpz.c)

#add_executable(polyselect-hd ${POL_SRC})

#target_link_libraries(polyselect-hd utils polyselect_common)

#set_property(TARGET polyselect-hd PROPERTY COMPILE_DEFINITIONS EEA_BOUND)

# Parameters selection

#set (PAR_SRC paraselect-hd.c parameters.c polynomials.c sieving_bound.c utils_mpz.c)

#add_executable(paraselect-hd ${PAR_SRC})

#target_link_libraries(paraselect-hd utils polyselect_common)

# Alpha 3d
set (ALP_SRC alpha3d.c utils_mpz_poly.c)
add_executable(alpha3d ${ALP_SRC})
set_property(TARGET alpha3d PROPERTY COMPILE_DEFINITIONS MAIN_ALPHA3D)
target_link_libraries(alpha3d utils)

# Murphy-E
set (MUR_SRC murphyE3d.c alpha3d.c utils_mpz_poly.c)
add_executable(murphyE3d ${MUR_SRC})
set_property(TARGET murphyE3d PROPERTY COMPILE_DEFINITIONS MAIN_MURPHYE3D)
target_link_libraries(murphyE3d utils polyselect_common)
