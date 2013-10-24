# RUSAGE_THREAD
message(STATUS "Testing whether RUSAGE_THREAD can be used")
try_compile(HAVE_RUSAGE_THREAD
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/rusage_thread.c)
if(HAVE_RUSAGE_THREAD)
   message(STATUS "Trying to compile with RUSAGE_THREAD -- Success")
else(HAVE_RUSAGE_THREAD)
   message(STATUS "Trying to compile with RUSAGE_THREAD -- Failed")
endif(HAVE_RUSAGE_THREAD)
