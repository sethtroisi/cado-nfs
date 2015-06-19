
# pclmul
message(STATUS "Testing whether pclmul code can be used")
try_run(pclmul_runs pclmul_compiles
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/pclmul.c)
if(pclmul_compiles)
    if (pclmul_runs MATCHES FAILED_TO_RUN)
        message(STATUS "Testing whether pclmul code can be used -- No")
        set (HAVE_PCLMUL 0)
    else (pclmul_runs MATCHES FAILED_TO_RUN)
        message(STATUS "Testing whether pclmul code can be used -- Yes")
        set (HAVE_PCLMUL 1)
    endif (pclmul_runs MATCHES FAILED_TO_RUN)
else(pclmul_compiles)
    try_run(pclmul_runs pclmul_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/pclmul.c
        COMPILE_DEFINITIONS -mpclmul)
    if(pclmul_compiles)
        if (pclmul_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether pclmul code can be used -- No")
            set (HAVE_PCLMUL 0)
        else (pclmul_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether pclmul code can be used -- Yes, with -mpclmul")
            set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mpclmul")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpclmul")
            set (HAVE_PCLMUL 1)
        endif (pclmul_runs MATCHES FAILED_TO_RUN)
    else(pclmul_compiles)
        message(STATUS "Testing whether pclmul code can be used -- No")
        set (HAVE_PCLMUL 0)
    endif(pclmul_compiles)
endif(pclmul_compiles)
