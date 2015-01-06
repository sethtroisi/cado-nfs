## TODO: one day, this should go up to the main CMakelists
FIND_LIBRARY(ECM_LIB ecm ${CMAKE_CURRENT_SOURCE_DIR})

set (SQ_SRC special-q_sieve.c mat_int64.c int64_vector.c uint64_array.c
    mat_Z.c ideal.c makefb.c utils_mpz.c int64_poly.c sieving_interval.c
    array.c utils_int64.c factor_base.c)

add_executable(special-q_sieve ${SQ_SRC})

#set_property(TARGET special-q_sieve PROPERTY COMPILE_DEFINITIONS TRACE_POS=3)
target_link_libraries(special-q_sieve ${ECM_LIB} utils)

