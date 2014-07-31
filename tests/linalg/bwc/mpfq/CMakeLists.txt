
include_directories(${PROJECT_SOURCE_DIR}/linalg/bwc/mpfq)

foreach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})

    set(extra)
    string(REGEX MATCH "p_([0-9]+)" t "${gfp_layer}")
    if(t)
        set(width ${CMAKE_MATCH_1})
    else(t)
        if(gfp_layer STREQUAL "pz")
            set(t TRUE)
            set(width "z")
        endif(gfp_layer STREQUAL "pz")
    endif(t)
    set(p "${prime_${width}}")

    add_executable(mpfq_test_${gfp_layer} mpfq_test_${gfp_layer}.c
        ${PROJECT_SOURCE_DIR}/linalg/bwc/mpfq/mpfq_${gfp_layer}.c)

    target_link_libraries(mpfq_test_${gfp_layer} gmp)
    run_test(mpfq_test_${gfp_layer} mpfq_test_${gfp_layer} -s 1 -p ${p})
endforeach(gfp_layer ${BWC_GFP_ARITHMETIC_BACKENDS})

