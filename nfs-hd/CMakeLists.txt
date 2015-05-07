# TODO: one day, this should go up to the main CMakelists
FIND_LIBRARY(ECM_LIB ecm ${CMAKE_CURRENT_SOURCE_DIR})

set (SQ_SRC special-q_sieve.c mat_int64.c int64_vector.c uint64_array.c
    mat_Z.c ideal.c makefb.c utils_mpz.c int64_poly.c sieving_bound.c
    array.c utils_int64.c factor_base.c utils_lattice.c double_vector.c
    mat_double.c)

add_executable(special-q_sieve ${SQ_SRC})

#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS ASSERT_FACTO
#    MEAN_NORM_BOUND MEAN_NORM_BOUND_SIEVE TEST_MQR NUMBER_HIT)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS COME_BACK_FK
#	TEST_MQR SPECIAL_Q_IDEAL_U)
set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS
	ASSERT_FACTO TEST_MQR SKEW_LLL)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS ASSERT_FACTO
#	NUMBER_HIT TEST_MQR)# TRACE_POS=0	ASSERT_SIEVE)
#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS ASSERT_FACTO
#	MEAN_NORM_BOUND MEAN_NORM_BOUND_SIEVE NUMBER_HIT NUMBER_SURVIVALS
#	PRINT_PARAMETERS TIMER_SIEVE=100 TEST_MQR TRACE_POS=0
#	MEAN_NORM ASSERT_SIEVE TRACE_POS=0 ENUM_LATTICE)

target_link_libraries(special-q_sieve ${ECM_LIB} utils)

# Makefb

set(MAKEFB_SRC makefb.c factor_base.c ideal.c utils_int64.c)

add_executable(makefb-hd ${MAKEFB_SRC})

set_property(TARGET makefb-hd PROPERTY COMPILE_DEFINITIONS MAIN)

target_link_libraries(makefb-hd ${ECM_LIB} utils)

# Polynomial selection

set (POL_SRC polyselect-hd.c polynomials.c sieving_bound.c utils_mpz.c)

add_executable(polyselect-hd ${POL_SRC})

target_link_libraries(polyselect-hd utils polyselect_common ${ECM_LIB})

#set_property(TARGET polyselect-hd PROPERTY COMPILE_DEFINITIONS EEA_BOUND)

# Parameters selection

set (PAR_SRC paraselect-hd.c parameters.c polynomials.c sieving_bound.c utils_mpz.c)

add_executable(paraselect-hd ${PAR_SRC})

target_link_libraries(paraselect-hd utils polyselect_common ${ECM_LIB})
