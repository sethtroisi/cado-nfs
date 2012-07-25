# Test the existence of __sync_fetch functions
message(STATUS " Test the existence of __sync_fetch functions")
try_compile(HAVE_SYNC_FETCH
    ${CADO_NFS_BINARY_DIR}/config
    ${CADO_NFS_SOURCE_DIR}/config/sync_fetch.c)
if(HAVE_SYNC_FETCH)
    message(STATUS "Trying to compile gcc-style x86_64 assembly -- Success")
else(HAVE_SYNC_FETCH)
    message(STATUS "Trying to compile gcc-style x86_64 assembly -- Failed")
endif(HAVE_SYNC_FETCH)

