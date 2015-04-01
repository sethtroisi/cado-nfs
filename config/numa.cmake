
# You can force a path to numa.h using the environment variables NUMA, or
# NUMA_INCDIR and NUMA_LIBDIR (note that curl/curl.h is also searched)

if (DEFINED ENV{NUMA})
    message(STATUS "Adding $ENV{NUMA} to the search path for numa")
    set(NUMA_INCDIR_HINTS ${NUMA_INCDIR_HINTS} "$ENV{NUMA}/include")
    set(NUMA_INCDIR_HINTS ${NUMA_INCDIR_HINTS} "$ENV{NUMA}"        )
    set(NUMA_LIBDIR_HINTS ${NUMA_LIBDIR_HINTS} "$ENV{NUMA}/lib"    )
    set(NUMA_LIBDIR_HINTS ${NUMA_LIBDIR_HINTS} "$ENV{NUMA}/.libs"  )
endif(DEFINED ENV{NUMA})

if (DEFINED ENV{NUMA_INCDIR})
    message(STATUS "Adding $ENV{NUMA_INCDIR} to the search path for numa")
    # prepend !
    set(NUMA_INCDIR_HINTS "$ENV{NUMA_INCDIR}" ${NUMA_INCDIR_HINTS})
endif (DEFINED ENV{NUMA_INCDIR})

if (DEFINED ENV{NUMA_LIBDIR})
    message(STATUS "Adding $ENV{NUMA_LIBDIR} to the search path for numa")
    # prepend !
    set(NUMA_LIBDIR_HINTS "$ENV{NUMA_LIBDIR}" ${NUMA_LIBDIR_HINTS})
endif(DEFINED ENV{NUMA_LIBDIR})

# Try in three passes, otherwise cmake gets in the way...
find_path   (NUMA_INCDIR numa.h HINTS ${NUMA_INCDIR_HINTS} DOC "numa headers"
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH
    )
if(NOT NUMA_INCDIR)
find_path   (NUMA_INCDIR numa.h HINTS ${NUMA_INCDIR_HINTS} DOC "numa headers"
        NO_DEFAULT_PATH
        )
endif(NOT NUMA_INCDIR)
if(NOT NUMA_INCDIR)
find_path   (NUMA_INCDIR numa.h HINTS ${NUMA_INCDIR_HINTS} DOC "numa headers")
endif(NOT NUMA_INCDIR)

find_library(NUMA_LIB    numa   HINTS ${NUMA_LIBDIR_HINTS} DOC "numa library"
    NO_DEFAULT_PATH
    )
if(NOT NUMA_LIB)
find_library(NUMA_LIB    numa   HINTS ${NUMA_LIBDIR_HINTS} DOC "numa library")
endif(NOT NUMA_LIB)

# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
message(STATUS "NUMA_INCDIR=${NUMA_INCDIR}")
message(STATUS "NUMA_LIB=${NUMA_LIB}")
string(COMPARE NOTEQUAL "${NUMA_INCDIR}" NUMA_INCDIR-NOTFOUND NUMA_INCDIR_OK)
string(COMPARE NOTEQUAL "${NUMA_LIB}" NUMA_LIB-NOTFOUND NUMA_LIBDIR_OK)

get_filename_component(NUMA_LIBDIR ${NUMA_LIB} PATH)

if(NUMA_INCDIR_OK AND NUMA_LIBDIR_OK)
include_directories(${NUMA_INCDIR})
link_directories(${NUMA_LIBDIR})
include(CheckCSourceCompiles)
set(CMAKE_REQUIRED_LIBRARIES "-lnuma")
set(CMAKE_REQUIRED_FLAGS "-L${NUMA_LIBDIR}")
set(CMAKE_REQUIRED_INCLUDES ${NUMA_INCDIR})
CHECK_C_SOURCE_COMPILES("
    #include <numa.h>
    int main(void)
    {
    return (numa_num_configured_cpus() << 16) + numa_num_configured_nodes();
    }
" HAVE_NUMA)
endif(NUMA_INCDIR_OK AND NUMA_LIBDIR_OK)
