
message(STATUS "Testing for log2")
try_compile(log2_compiles
            ${CADO_NFS_BINARY_DIR}/config
            ${CADO_NFS_SOURCE_DIR}/config/log2.c)
if(log2_compiles)
message(STATUS "Testing for log2 -- found")
set(HAVE_LOG2 1)
else(log2_compiles)
message(STATUS "Testing for log2 -- not found")
set(HAVE_LOG2 0)
endif(log2_compiles)
