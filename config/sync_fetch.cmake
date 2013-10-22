#Testing the existence of __sync_fetch functions
message(STATUS "Testing the existence of __sync_fetch functions")
try_compile(HAVE_SYNC_FETCH
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/sync_fetch.c)
if(HAVE_SYNC_FETCH)
    message(STATUS "Testing the existence of __sync_fetch functions -- Success")
else(HAVE_SYNC_FETCH)
    message(STATUS "Testing the existence of __sync_fetch functions -- Failed")
endif(HAVE_SYNC_FETCH)

