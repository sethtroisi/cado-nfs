
set (sources 
        cado_poly.c relation.c timing.cpp plain_poly.c getprime.c mod_ul.c
        modredc_ul.c modredc_15ul.c modredc_2ul2.c modul_poly.c gmp_aux.c
        mpz_poly.c double_poly.c rootfinder.c params.c gcd.c
        discriminant.c mpz_array.c gzip.c misc.c
        poly.c crc.c purgedfile.c bit_vector.c
        barrier.c usp.c fix-endianness.c memusage.c
        cachesize_cpuid.c cachesize_guess.c roots_mod.c mod_mpz.c renumber.c
        ringbuf.c
)
if (NOT MINGW)
    set (sources ${sources} cado_popen.c)
endif (NOT MINGW)

# the utils library is used by some shared objects (prominently the
# matmul shared libs in bwc), and therefore must be shared.
add_library (utils ${CADO_SHARED} ${sources})
target_link_libraries(utils pthread ${gmp_libname} m)
install(TARGETS utils DESTINATION lib/${package_basename})

add_executable(cut_n_roots cut_n_roots.c)
target_link_libraries(cut_n_roots utils)

add_executable(roots roots.c)
target_link_libraries(roots utils)

add_executable(test_mod_ul ${EXCLUDE_FROM_RELEASE_BUILD} test_mod_ul.c)
add_executable(test_modredc_ul ${EXCLUDE_FROM_RELEASE_BUILD} test_modredc_ul.c)
add_executable(test_modredc_15ul ${EXCLUDE_FROM_RELEASE_BUILD} test_modredc_15ul.c)
add_executable(test_modredc_2ul2 ${EXCLUDE_FROM_RELEASE_BUILD} test_modredc_2ul2.c)
add_executable(test-rootfinder ${EXCLUDE_FROM_RELEASE_BUILD} test-rootfinder.c)
add_executable(test_roots_mod ${EXCLUDE_FROM_RELEASE_BUILD} test_roots_mod.c)
add_executable(test_omega ${EXCLUDE_FROM_RELEASE_BUILD} test_omega.c)
add_executable(gmp_prob_prime gmp_prob_prime.c)
add_executable(gmp_nextprime gmp_nextprime.c)
add_executable(random_integer random_integer.c)

# To avoid the warning: implicit declaration of function nanosleep. It's POSIX, non c99.
set (CMAKE_OLD_C_FLAGS "${CMAKE_C_FLAGS}")
set (CMAKE_C_FLAGS "prout -D _POSIX_C_SOURCE=200809L ${CMAKE_C_FLAGS}")
add_executable(antebuffer antebuffer.c timing.cpp memusage.c)
set (CMAKE_C_FLAGS "${CMAKE_OLD_C_FLAGS}")

target_link_libraries(test_mod_ul utils)
target_link_libraries(test_modredc_ul utils)
target_link_libraries(test_modredc_15ul utils)
target_link_libraries(test_modredc_2ul2 utils)
target_link_libraries(test-rootfinder utils)
target_link_libraries(test_roots_mod utils)
target_link_libraries(gmp_prob_prime ${gmp_libname})
target_link_libraries(gmp_nextprime ${gmp_libname})
target_link_libraries(random_integer ${gmp_libname})
target_link_libraries(test_omega utils)
target_link_libraries(antebuffer pthread)

include(CMakeLists-nodist.txt OPTIONAL)

install(TARGETS gmp_prob_prime random_integer RUNTIME DESTINATION bin/utils)
install(TARGETS cut_n_roots roots gmp_nextprime antebuffer RUNTIME DESTINATION bin/utils)
