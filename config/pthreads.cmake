include(${PROJECT_SOURCE_DIR}/config/search_for_function.cmake)

set(CMAKE_REQUIRED_LIBRARIES)
search_for_function(pthread_create HAVE_PTHREAD_CREATE pthread)
if(HAVE_PTHREAD_CREATE)
    # OK. Assume that we have the bare minimum for using threads, falling
    # back on workalikes for barrier synchronization waits if needed
    # (like we used to do in the past, anyway). Thus we can already set
    # the proper flags.  We're using WITH_PTHREADS in the top-level
    # substitution, so it needs to escape its scope and go into the cache
    # right now.
    set(WITH_PTHREADS 1 CACHE INTERNAL "pthreads are being used")
    add_definitions(-DWITH_PTHREADS)
    set(pthread_libs ${CMAKE_REQUIRED_LIBRARIES})
    search_for_function(pthread_barrier_wait HAVE_PTHREAD_BARRIER_WAIT rt)
    if(HAVE_PTHREAD_BARRIER_WAIT)
        set(pthread_libs ${CMAKE_REQUIRED_LIBRARIES})
    endif(HAVE_PTHREAD_BARRIER_WAIT)
else(HAVE_PTHREAD_CREATE)
    message(FATAL_ERROR "POSIX threads not found. Basic thread support is required by cado-nfs")
endif(HAVE_PTHREAD_CREATE)

