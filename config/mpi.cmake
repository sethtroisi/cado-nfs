
#############################################################
# mpi
# Don't use the FindMPI module, it's buggy.

if($ENV{MPI})
    set(WANT_MPI 1)
else($ENV{MPI})
    # cmake claims this:
#         if(variable)
#
#       True if the variable's value is not empty, 0, N, NO, OFF, FALSE,
#       NOTFOUND, or <variable>-NOTFOUND.
    # However it's not true, and a path evaluates to false.
    string(REGEX MATCH "^/.*$" WANT_MPI "$ENV{MPI}")
    if(WANT_MPI STREQUAL "")
        set(WANT_MPI 0)
    else(WANT_MPI STREQUAL "")
        set(WANT_MPI 1)
    endif(WANT_MPI STREQUAL "")
endif($ENV{MPI})

if(${WANT_MPI})
    set(findprog_flags
        NO_DEFAULT_PATH
        NO_CMAKE_ENVIRONMENT_PATH
        NO_CMAKE_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH)
    if("$ENV{MPI}" MATCHES "^(1|YES|yes|ON|on|)$")
        set(findprog_flags)
    endif("$ENV{MPI}" MATCHES "^(1|YES|yes|ON|on|)$")
    find_program(MPI_C_COMPILER mpicc HINTS $ENV{MPI} "$ENV{MPI}/bin"
        ${findprog_flags})
    find_program(MPI_CXX_COMPILER
        NAMES mpic++ mpicxx mpiCC
        HINTS $ENV{MPI} "$ENV{MPI}/bin"
        ${findprog_flags})
    find_program(MPIEXEC
        mpiexec
        HINTS $ENV{MPI} "$ENV{MPI}/bin"
        ${findprog_flags})
    if (MPI_C_COMPILER AND MPI_CXX_COMPILER AND MPIEXEC)
        message(STATUS "Using MPI C compiler ${MPI_C_COMPILER}")
        message(STATUS "Using MPI C++ compiler ${MPI_CXX_COMPILER}")
        message(STATUS "Using MPI driver ${MPIEXEC}")
        get_filename_component(HAVE_MPI ${MPIEXEC} PATH)
        # We're using this variable in the top-level substitution, so it needs
        # to escape its scope and go into the cache right now.
        set(WITH_MPI 1 CACHE INTERNAL "MPI is being used (for relevant code parts)")
    else(MPI_C_COMPILER AND MPI_CXX_COMPILER AND MPIEXEC)
        message(FATAL_ERROR "Cannot find all of mpicc/mpic++/mpiexec with MPI=$ENV{MPI}")
    endif(MPI_C_COMPILER AND MPI_CXX_COMPILER AND MPIEXEC)
else(${WANT_MPI})
    message(STATUS "MPI is not enabled")
    set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
    set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
endif(${WANT_MPI})


