string(COMPARE NOTEQUAL "$ENV{ECM}" "" HAS_ECM_OVERRIDE)
if (HAS_ECM_OVERRIDE)
    message(STATUS "Adding $ENV{ECM} to the search path")
    set(ECM_INCDIR_HINTS "$ENV{ECM}/include" ${ECM_INCDIR_HINTS})
    set(ECM_LIBDIR_HINTS "$ENV{ECM}/lib"     ${ECM_LIBDIR_HINTS})

    find_path(ECM_INCDIR ecm.h PATHS ${ECM_INCDIR_HINTS} DOC "GMP-ECM headers"
            NO_DEFAULT_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_PATH
            NO_CMAKE_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            NO_CMAKE_FIND_ROOT_PATH)

    find_library(ECM_LIB    ecm   HINTS ${ECM_LIBDIR_HINTS} DOC "GMP-ECM library" NO_DEFAULT_PATH)

    # Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
    # helpful documentation there :-((
    message(STATUS "ECM_INCDIR=${ECM_INCDIR}")
    message(STATUS "ECM_LIBDIR=${ECM_LIBDIR}")
    if(ECM_INCDIR)
    include_directories(${ECM_INCDIR})
    else(ECM_INCDIR)
    message(FATAL_ERROR "ecm.h cannot be found. Please install GMP-ECM, and specify its install prefix in local.sh")
    endif(ECM_INCDIR)
    if(ECM_LIBDIR)
    link_directories(${ECM_LIBDIR})
    else(ECM_LIBDIR)
    message(FATAL_ERROR "libecm cannot be found. Please install GMP-ECM, and specify its install prefix in local.sh")
    endif(ECM_LIBDIR)
endif(HAS_ECM_OVERRIDE)
