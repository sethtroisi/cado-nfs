
SET (sources 
        cado_poly.c relation.c timing.c plain_poly.c getprime.c mod_ul.c
        modredc_ul.c modredc_15ul.c modredc_2ul2.c modul_poly.c gmp_aux.c
        mpz_poly.c fpoly.c rootfinder.c params.c gcd_int64.c gcd_uint64.c
        discriminant.c random_generation.c mpz_array.c gzip.c misc.c
)

ADD_LIBRARY (utils ${sources})
target_link_libraries(utils gmp)

ADD_EXECUTABLE(cache cache.c timing.c)
ADD_EXECUTABLE(cachesize cachesize.c)
ADD_EXECUTABLE(check_rels check_rels.c)
target_link_libraries(check_rels utils)

macro(add_alex_test mayberedc typewidth)
add_executable(test_mod_${mayberedc}${typewidth} test_mod.c)
target_link_libraries(test_mod_${mayberedc}${typewidth} utils)
set_property(TARGET test_mod_${mayberedc}${typewidth} APPEND PROPERTY
    COMPILE_FLAGS
    -DARITHMETIC=\"\\\"mod${mayberedc}_${typewidth}_default.h\\\"\")
endmacro(add_alex_test)

add_executable(test-rootfinder test-rootfinder.c)
target_link_libraries(test-rootfinder utils)

add_alex_test("" ul)
add_alex_test(redc ul)
add_alex_test(redc 15ul)
add_alex_test(redc 2ul2)

# test_mod_ul
# test_mod_redcul
# test_mod_redc15ul
# test_mod_redc2ul2
# test-rootfinder-script.sh < testbase-rootfinder.quick

