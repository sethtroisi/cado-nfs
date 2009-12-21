
set (sources 
        cado_poly.c relation.c timing.c plain_poly.c getprime.c mod_ul.c
        modredc_ul.c modredc_15ul.c modredc_2ul2.c modul_poly.c gmp_aux.c
        mpz_poly.c fpoly.c rootfinder.c params.c gcd_int64.c gcd_uint64.c
        discriminant.c random_generation.c mpz_array.c gzip.c misc.c
        hashpair.c poly.c crc.c gmp_gcd.c gmp_factorization_complete.c
		gmp_prob_prime.c
)

add_library (utils ${sources})
target_link_libraries(utils gmp)


add_executable(check_rels check_rels.c)
install(TARGETS check_rels RUNTIME DESTINATION bin/utils)
target_link_libraries(check_rels utils)

add_executable(test_mod_ul test_mod_ul.c)
add_executable(test_modredc_ul test_modredc_ul.c)
add_executable(test_modredc_15ul test_modredc_15ul.c)
add_executable(test_modredc_2ul2 test_modredc_2ul2.c)
add_executable(gmp_gcd gmp_gcd.c)
add_executable(gmp_factorization_complete gmp_factorization_complete.c)
add_executable(gmp_prob_prime gmp_prob_prime.c)

target_link_libraries(test_mod_ul utils)
target_link_libraries(test_modredc_ul utils)
target_link_libraries(test_modredc_15ul utils)
target_link_libraries(test_modredc_2ul2 utils)
target_link_libraries(gmp_gcd gmp)
target_link_libraries(gmp_factorization_complete gmp)
target_link_libraries(gmp_prob_prime gmp)

include(CMakeLists-nodist.txt OPTIONAL)
