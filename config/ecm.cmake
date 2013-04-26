string(COMPARE NOTEQUAL "$ENV{ECM}" "" HAS_ECM_OVERRIDE)
if (HAS_ECM_OVERRIDE)
    message(STATUS "Adding $ENV{ECM} to the search path")
    set(ECM_INCDIR_HINTS "$ENV{ECM}/include" ${ECM_INCDIR_HINTS})
    set(ECM_LIBDIR_HINTS "$ENV{ECM}/lib"     ${ECM_LIBDIR_HINTS})
endif(HAS_ECM_OVERRIDE)

set(ECM_INCDIR "$ENV{ECM}/include")
set(ECM_LIBDIR "$ENV{ECM}/lib")


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
message(FATAL_ERROR "ecm.h cannot be found. Please install GMP-ECM, and specify its install prefix in local.sh")
endif(ECM_LIBDIR)

