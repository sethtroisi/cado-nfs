
message(STATUS "Testing for cabsl")
try_compile(cabsl_compiles
            ${CADO_NFS_BINARY_DIR}/config
            ${CADO_NFS_SOURCE_DIR}/config/cabsl.c)
if(cabsl_compiles)
message(STATUS "Testing for cabsl -- found")
set(HAVE_CABSL 1)
else(cabsl_compiles)
message(STATUS "Testing for cabsl -- not found")
set(HAVE_CABSL 0)
endif(cabsl_compiles)
