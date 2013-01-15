# Test for some POSIX headers and functions that may not exist on non-Unix systems

INCLUDE (CheckIncludeFiles)

CHECK_INCLUDE_FILES (sys/resource.h HAVE_RESOURCE_H)

# Unset the CMake variable that search_for_function() interprets
set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_DEFINITIONS)
set(CMAKE_REQUIRED_INCLUDES)
set(CMAKE_REQUIRED_LIBRARIES)

include(${CADO_NFS_SOURCE_DIR}/config/search_for_function.cmake)
search_for_function(getc_unlocked HAVE_GETC_UNLOCKED)
search_for_function(nanosleep HAVE_NANOSLEEP)
search_for_function(usleep HAVE_USLEEP)
search_for_function(popen HAVE_POPEN)
search_for_function(pclose HAVE_PCLOSE)
search_for_function(getrusage HAVE_GETRUSAGE)
search_for_function(lrand48 HAVE_LRAND48)
