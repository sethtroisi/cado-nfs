
# SSE-2
message(STATUS "Testing whether sse-2 code can be used")
try_compile(sse2_compiles
            ${CADO_NFS_BINARY_DIR}/config
            ${CADO_NFS_SOURCE_DIR}/config/sse2.c)
if(sse2_compiles)
    message(STATUS "Testing whether sse-2 code can be used -- Yes")
else(sse2_compiles)
    try_compile(sse2_compiles
        ${CADO_NFS_BINARY_DIR}/config
        ${CADO_NFS_SOURCE_DIR}/config/sse2.c
        COMPILE_DEFINITIONS -msse2)
    if(sse2_compiles)
        message(STATUS "Trying whether sse-2 code can be used -- Yes, with -msse2")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    else(sse2_compiles)
        message(STATUS "Trying whether sse-2 code can be used -- No")
    endif(sse2_compiles)
endif(sse2_compiles)

