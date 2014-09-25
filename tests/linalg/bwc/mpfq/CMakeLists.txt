
include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc/mpfq)

foreach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})

    add_executable(mpfq_test_${gfp_layer} mpfq_test_${gfp_layer}.c
        ${PROJECT_SOURCE_DIR}/linalg/bwc/mpfq/mpfq_${gfp_layer}.c)

    set(extra)
    string(REGEX MATCH "p_([0-9]+)" t "${gfp_layer}")
    if(t)
        set(width ${CMAKE_MATCH_1})
    else(t)
        if(gfp_layer STREQUAL "pz")
            set(t TRUE)
            set(width "z")
            set_target_properties(mpfq_test_${gfp_layer} PROPERTIES COMPILE_FLAGS "-DVARIABLE_SIZE_PRIME")
        endif(gfp_layer STREQUAL "pz")
    endif(t)
    set(p "${prime_${width}}")


    target_link_libraries(mpfq_test_${gfp_layer} ${gmp_libname})
    run_test(mpfq_test_${gfp_layer} mpfq_test_${gfp_layer} -s 1 -p ${p})
endforeach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})

# mpfq_test_fixmp is not seeded as of now, because it uses random2, and
# there is no mpn_rrandomb as of gmp-6.0.0
add_executable(mpfq_test_fixmp mpfq_test_fixmp.c)
target_link_libraries(mpfq_test_fixmp ${gmp_libname})
run_test(mpfq_test_fixmp mpfq_test_fixmp 40)
