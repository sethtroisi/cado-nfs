
set (SQ_SRC special-q_sieve.c mat_int64.c int64_vector.c uint64_array.c
    mat_Z.c ideal.c makefb.c utils_mpz.c int64_poly.c sieving_bound.c
    array.c utils_int64.c factor_base.c utils_lattice.c double_vector.c
    mat_double.c list_int64_vector.c list_double_vector.c
    list_int64_vector_index.c)

add_executable(special-q_sieve ${SQ_SRC})

#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS ASSERT_FACTO
#    MEAN_NORM_BOUND MEAN_NORM_BOUND_SIEVE TEST_MQR NUMBER_HIT)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS COME_BACK_FK
#	TEST_MQR SPECIAL_Q_IDEAL_U)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS
	#ASSERT_FACTO TEST_MQR SKEW_LLL)
set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS
	SPACE_SIEVE SPACE_SIEVE_CONTRIBUTION NEW_NORM)
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

set(MAKEFB_SRC makefb.c factor_base.c ideal.c utils_int64.c)

add_executable(makefb-hd ${MAKEFB_SRC})

set_property(TARGET makefb-hd PROPERTY COMPILE_DEFINITIONS MAIN)

target_link_libraries(makefb-hd utils)
