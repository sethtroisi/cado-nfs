
# SSSE-3
message(STATUS "Testing whether ssse-3 code can be used")
# We used to base this test on sse-3, but apparently there are some
# compiler / binutils versions (gcc-4.7.2 on x86_64-unknown-openbsd5.3,
# binutils 2.15) for which ssse3 and sse4.1 fail with no apparent
# explanation.
if (HAVE_SSE2)
    try_run(ssse3_runs ssse3_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/ssse3.c)
    if(ssse3_compiles)
        if (ssse3_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether ssse-3 code can be used -- No")
            set (HAVE_SSSE3 0)
        else (ssse3_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether ssse-3 code can be used -- Yes")
            set (HAVE_SSSE3 1)
        endif (ssse3_runs MATCHES FAILED_TO_RUN)
    else(ssse3_compiles)
        try_run(ssse3_runs ssse3_compiles
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/ssse3.c
            COMPILE_DEFINITIONS -mssse3)
        if(ssse3_compiles)
            if (ssse3_runs MATCHES FAILED_TO_RUN)
                message(STATUS "Testing whether ssse-3 code can be used -- No")
                set (HAVE_SSSE3 0)
            else (ssse3_runs MATCHES FAILED_TO_RUN)
                message(STATUS "Testing whether ssse-3 code can be used -- Yes, with -mssse3")
                set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mssse3")
                set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mssse3")
                set (HAVE_SSSE3 1)
            endif (ssse3_runs MATCHES FAILED_TO_RUN)
        else(ssse3_compiles)
            message(STATUS "Testing whether ssse-3 code can be used -- No")
            set (HAVE_SSSE3 0)
        endif(ssse3_compiles)
    endif(ssse3_compiles)
else (HAVE_SSE2)
    message(STATUS "Testing whether ssse-3 code can be used -- skipped")
endif (HAVE_SSE2)
