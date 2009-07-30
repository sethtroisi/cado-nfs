
include_directories(..)

include(${CMAKE_SOURCE_DIR}/config/search_for_function.cmake)

#############################################################
# pthreads
if($ENV{PTHREADS})
    set(CMAKE_REQUIRED_LIBRARIES)
    search_for_function(pthread_create HAVE_PTHREAD_CREATE -lpthread)
    if(HAVE_PTHREAD_CREATE)
        # OK. Assume that we have the bare minimum for using threads, falling
        # back on workalikeas for barrier synchronization waits if needed
        # (like we used to do in the past, anyway). Thus we can already set
        # the proper flags.  We're using WITH_PTHREADS in the top-level
        # substitution, so it needs to escape its scope and go into the cache
        # right now.
        set(WITH_PTHREADS 1 CACHE INTERNAL "pthreads are being used")
        add_definitions(-DWITH_PTHREADS)
        set(pthread_libs ${pthread_libs} ${CMAKE_REQUIRED_LIBRARIES})
        set(CMAKE_REQUIRED_LIBRARIES)
        search_for_function(pthread_barrier_wait HAVE_PTHREAD_BARRIER_WAIT -lrt)
        if(HAVE_PTHREAD_BARRIER_WAIT)
            set(pthread_libs ${pthread_libs} ${CMAKE_REQUIRED_LIBRARIES})
        endif(HAVE_PTHREAD_BARRIER_WAIT)
    endif(HAVE_PTHREAD_CREATE)
else($ENV{PTHREADS})
    message(STATUS "PTHREADS are not enabled")
endif($ENV{PTHREADS})

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
        set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
        set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
        get_filename_component(HAVE_MPI ${MPIEXEC} PATH)
        # We're using this variable in the top-level substitution, so it needs
        # to escape its scope and go into the cache right now.
        set(WITH_MPI 1 CACHE INTERNAL "MPI is being used")
        add_definitions(-DWITH_MPI)
    else(MPI_C_COMPILER AND MPI_CXX_COMPILER AND MPIEXEC)
        message(FATAL_ERROR "Cannot find all of mpicc/mpic++/mpiexec with MPI=$ENV{MPI}")
    endif(MPI_C_COMPILER AND MPI_CXX_COMPILER AND MPIEXEC)
else(${WANT_MPI})
    message(STATUS "MPI is not enabled")
endif(${WANT_MPI})


configure_file(bwc.pl ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl ESCAPE_QUOTES @ONLY)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
configure_file(bwc_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/bwc_config.h ESCAPE_QUOTES @ONLY)

add_executable(test_hello EXCLUDE_FROM_ALL
    test_hello.c parallelizing_info.c barrier.c)
add_executable(random random.c)
add_executable(split split.c bw-common.c)
add_executable(shell
    shell.c parallelizing_info.c bw-common.c bw-common-mpi.c barrier.c)
add_executable(acollect acollect.c bw-common.c)

install(TARGETS random split shell acollect RUNTIME DESTINATION bin/linalg/bwc)

target_link_libraries(test_hello utils)
target_link_libraries(random utils m)
target_link_libraries(split utils)
target_link_libraries(shell utils ${pthread_libs})
target_link_libraries(acollect utils)

add_library(bwc_base
    parallelizing_info.c
    intersections.cpp
    xvectors.c
    bw-common.c
    bw-common-mpi.c
    rusage.c
    async.c
    debug.c
    barrier.c
    worker-threads.c)
target_link_libraries(bwc_base ${pthread_libs})

macro(push_abase_cflags v target)
    set_property(TARGET ${v}_${target}
        APPEND PROPERTY COMPILE_FLAGS -DSELECT_ABASE_${v})
endmacro(push_abase_cflags)

set(all_variants)
macro(configure_bwc_variant v)
    set(${v}_matmul_sources 
        matmul-common.c matmul-basic.c matmul-sliced.cpp
        matmul-threaded.c matmul-bucket.cpp matmul.c)

    string(COMPARE EQUAL "${v}" "u64" doing_u64)
    if(doing_u64)
        if(HAVE_GCC_STYLE_AMD64_ASM)
            set(${v}_matmul_sources ${${v}_matmul_sources}
                    matmul-sub-small1.S
                    matmul-sub-small2.S
                    matmul-sub-large-fbi.S
                    matmul-sub-large-fbd.S
                    matmul-sub-vsc-combine.S
                    matmul-sub-vsc-dispatch.S
                    )
            set_property(SOURCE matmul-sub-small1.S PROPERTY LANGUAGE C)
            set_property(SOURCE matmul-sub-small2.S PROPERTY LANGUAGE C)
            set_property(SOURCE matmul-sub-large-fbi.S PROPERTY LANGUAGE C)
            set_property(SOURCE matmul-sub-large-fbd.S PROPERTY LANGUAGE C)
            set_property(SOURCE matmul-sub-vsc-dispatch.S PROPERTY LANGUAGE C)
            set_property(SOURCE matmul-sub-vsc-combine.S PROPERTY LANGUAGE C)
        endif(HAVE_GCC_STYLE_AMD64_ASM)
    endif(doing_u64)

    add_library(${v}_matmul ${${v}_matmul_sources})
    target_link_libraries(${v}_matmul bwc_base)

    push_abase_cflags(${v} matmul)
    target_link_libraries(${v}_matmul mat_io)

    add_library(${v}_bwc abase-${v}.c info_file.c matmul_top.c)
    push_abase_cflags(${v} bwc)

    target_link_libraries(${v}_bwc bwc_base ${v}_matmul utils pthread)

    set(all_variants ${all_variants} ${v})
endmacro(configure_bwc_variant)

configure_bwc_variant(u64)
configure_bwc_variant(u64k)
configure_bwc_variant(u64n)
configure_bwc_variant(u128)

set(build_objs build.c)
set(bench_objs bench.c)
set(prep_objs prep.c xymats.c ../gauss.c)
set_source_files_properties (../gauss.c COMPILE_FLAGS "-DMULTI_ROW=3 -DNO_MAIN")
set(krylov_objs krylov.c xymats.c xdotprod.c)
set(mksol_objs mksol.c xymats.c xdotprod.c)
set(gather_objs gather.c)
set(secure_objs secure.c)

macro(configure_program prg)
    foreach(v ${ARGN})
        add_executable(${v}_${prg} ${${prg}_objs})
        install(TARGETS ${v}_${prg} RUNTIME DESTINATION bin/linalg/bwc)
        target_link_libraries(${v}_${prg} ${v}_bwc)
        push_abase_cflags(${v} ${prg})
        set_target_properties(${v}_${prg} PROPERTIES LINKER_LANGUAGE CXX)
    endforeach(v)
endmacro(configure_program)

configure_program(krylov ${all_variants})
configure_program(mksol  ${all_variants})
configure_program(build  ${all_variants})
configure_program(bench  ${all_variants})
configure_program(prep   u64n)
configure_program(gather u64n)
configure_program(secure u64)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl DESTINATION bin/linalg/bwc)

# Also lingen.

include_directories(${CADO_NFS_SOURCE_DIR}/gf2x-fft)
add_executable(lingen lingen_binary.cpp)
install(TARGETS lingen RUNTIME DESTINATION bin/linalg/bwc)
target_link_libraries(lingen bwc_base utils gf2x-fft)

include(CMakeLists-nodist.txt OPTIONAL)
