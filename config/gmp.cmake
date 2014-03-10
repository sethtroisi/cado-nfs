
# You can force a path to gmp.h using the environment variables GMP, or
# GMP_INCDIR and GMP_LIBDIR
string(COMPARE NOTEQUAL "$ENV{GMP}" "" HAS_GMP_OVERRIDE)
if (HAS_GMP_OVERRIDE)
    message(STATUS "Adding $ENV{GMP} to the search path for Gnu MP")
    set(GMP_INCDIR_HINTS "$ENV{GMP}/include" ${GMP_INCDIR_HINTS})
    set(GMP_INCDIR_HINTS "$ENV{GMP}"         ${GMP_INCDIR_HINTS})
    set(GMP_LIBDIR_HINTS "$ENV{GMP}/lib"     ${GMP_LIBDIR_HINTS})
    set(GMP_LIBDIR_HINTS "$ENV{GMP}/.libs"   ${GMP_LIBDIR_HINTS})
endif(HAS_GMP_OVERRIDE)
string(COMPARE NOTEQUAL "$ENV{GMP_INCDIR}" "" HAS_GMP_INCDIR_OVERRIDE)
if (HAS_GMP_INCDIR_OVERRIDE)
    message(STATUS "Adding $ENV{GMP_INCDIR} to the search path for Gnu MP")
    set(GMP_INCDIR_HINTS "$ENV{GMP_INCDIR}" ${GMP_INCDIR_HINTS})
endif(HAS_GMP_INCDIR_OVERRIDE)
string(COMPARE NOTEQUAL "$ENV{GMP_LIBDIR}" "" HAS_GMP_LIBDIR_OVERRIDE)
if (HAS_GMP_LIBDIR_OVERRIDE)
    message(STATUS "Adding $ENV{GMP_LIBDIR} to the search path for Gnu MP")
    set(GMP_LIBDIR_HINTS "$ENV{GMP_LIBDIR}"     ${GMP_LIBDIR_HINTS})
endif(HAS_GMP_LIBDIR_OVERRIDE)

# First try overrides, really. We want cmake to shut up.
if (NOT GMP_INCDIR)
find_path   (GMP_INCDIR gmp.h PATHS ${GMP_INCDIR_HINTS} DOC "Gnu MP headers"
        NO_DEFAULT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        NO_CMAKE_FIND_ROOT_PATH)
endif(NOT GMP_INCDIR)
if (NOT GMP_INCDIR)
find_path   (GMP_INCDIR gmp.h HINTS ${GMP_INCDIR_HINTS} DOC "Gnu MP headers"
        NO_DEFAULT_PATH
    )
endif(NOT GMP_INCDIR)
if (NOT GMP_INCDIR)
find_path   (GMP_INCDIR gmp.h HINTS ${GMP_INCDIR_HINTS} DOC "Gnu MP headers")
endif(NOT GMP_INCDIR)

find_library(GMP_LIB    gmp   HINTS ${GMP_LIBDIR_HINTS} DOC "Gnu MP library" NO_DEFAULT_PATH)
if(NOT GMP_LIBDIR)
find_library(GMP_LIB    gmp   HINTS ${GMP_LIBDIR_HINTS} DOC "Gnu MP library")
endif(NOT GMP_LIBDIR)

# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
get_filename_component(GMP_LIBDIR ${GMP_LIB} PATH)
message(STATUS "GMP_INCDIR=${GMP_INCDIR}")
message(STATUS "GMP_LIBDIR=${GMP_LIBDIR}")
if(GMP_INCDIR)
include_directories(${GMP_INCDIR})
else(GMP_INCDIR)
message(FATAL_ERROR "gmp.h cannot be found. Please install Gnu MP, and specify its install prefix in local.sh")
endif(GMP_INCDIR)
if(GMP_LIBDIR)
link_directories(${GMP_LIBDIR})
else(GMP_LIBDIR)
message(FATAL_ERROR "GMP library cannot be found. Please install Gnu MP, and specify its install prefix in local.sh")
endif(GMP_LIBDIR)

set(gmp_libname "gmp")
