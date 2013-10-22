#Testing the existence of asprintf/vasprintf functions
message(STATUS "Testing the existence of asprintf/vasprintf functions")
try_compile(HAVE_ASPRINTF
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/asprintf.c)
if(HAVE_ASPRINTF)
    message(STATUS "Testing the existence of asprintf/vasprintf functions -- Success")
else(HAVE_asprintf)
    message(STATUS "Testing the existence of asprintf/vasprintf functions -- Failed")
endif(HAVE_ASPRINTF)

