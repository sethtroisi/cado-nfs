
include_directories(. ..)
# include_directories(mpfq)

# Use with care. Normally, binaries not needed MPI for real should not
# emit MPI symbols since they don't do anything depending on this config
# variable. However this mechanism is somewhat fragile.
if (${WITH_MPI})
add_definitions(-DWITH_MPI)
endif (${WITH_MPI})

add_library(mf mf.c)
target_link_libraries(mf balancing utils)
add_executable(mf_scan ${EXCLUDE_FROM_RELEASE_BUILD} mf_scan.c)
target_link_libraries(mf_scan mf)

include_directories(${CADO_NFS_SOURCE_DIR}/linalg)
add_executable(mf_bal mf_bal.c ../rowset_heap.cpp)
target_link_libraries(mf_bal mf m)
install(TARGETS mf_bal RUNTIME DESTINATION bin/linalg/bwc)

configure_file(bwc.pl ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl ESCAPE_QUOTES @ONLY)

add_executable(random_matrix ${EXCLUDE_FROM_RELEASE_BUILD} random_matrix.c)

# balancing.c itself is in fact standalone, and lightweight. We prefer to
# avoid having the whole thing entangled with too many dependencies all
# around, since all that split.c actually needs is very, very little.
add_executable(split split.c bw-common.c balancing.c)

add_executable(acollect acollect.c bw-common.c)

install(TARGETS split acollect RUNTIME DESTINATION bin/linalg/bwc)

target_link_libraries(random_matrix utils m)
target_link_libraries(split utils)
target_link_libraries(acollect utils)

set(balancing_sources
    balancing.c
    balancing_workhorse.c
    balancing_file_source.c
    balancing_mpi_source.c)

if (HAVE_CURL)
set(balancing_sources ${balancing_sources} balancing_curl_source.c)
endif (HAVE_CURL)

add_library(balancing ${balancing_sources})

set_target_properties(balancing PROPERTIES
        COMPILE_FLAGS "--mpi"
        LINK_FLAGS "--mpi")

add_executable(bcast-file ${EXCLUDE_FROM_RELEASE_BUILD} bcast-file.c)
set_target_properties(bcast-file PROPERTIES
        COMPILE_FLAGS "--mpi"
        LINK_FLAGS "--mpi")

target_link_libraries(balancing utils)
if (HAVE_CURL)
target_link_libraries(balancing curl)
endif (HAVE_CURL)

set(bwc_abases u64k1 u64k2)

set(bwc_mpfq_sources mpfq/abase_vbase.c)
foreach(v ${bwc_abases})
    set(bwc_mpfq_sources ${bwc_mpfq_sources} mpfq/abase_${v}.c mpfq/abase_${v}_t.c )
endforeach(v)


if (ENABLE_MPFQ_PRIME_FIELDS_FOR_DLOG)
foreach(v p_1 p_2 p_3 p_4 p_8)
    set(bwc_mpfq_sources ${bwc_mpfq_sources} mpfq/abase_${v}.c mpfq/abase_${v}_t.c )
endforeach(v)
endif (ENABLE_MPFQ_PRIME_FIELDS_FOR_DLOG)

set(bwc_base_sources
    parallelizing_info.c
    intersections.c
    xvectors.c
    bw-common.c
    bw-common-mpi.c
    async.c
    rolling.c
    matmul.c
    matmul_top.c
)

add_library(bwc_base ${bwc_base_sources})
add_library(bwc_mpfq ${CADO_SHARED} ${bwc_mpfq_sources})
install(TARGETS bwc_mpfq DESTINATION "lib/${package_basename}")

IF(UNIX)
FIND_PROGRAM(CMAKE_UNAME uname /bin /usr/bin /usr/local/bin )
IF(CMAKE_UNAME)
EXEC_PROGRAM(uname ARGS -s OUTPUT_VARIABLE CMAKE_SYSTEM_NAME)
ENDIF(CMAKE_UNAME)
ENDIF(UNIX)
IF(CMAKE_SYSTEM_NAME MATCHES BSD)
target_link_libraries(bwc_base balancing utils bwc_mpfq ${pthread_libs})
ELSE(CMAKE_SYSTEM_NAME MATCHES BSD)
if(BUILD_DYNAMICALLY_LINKABLE_BWC)
target_link_libraries(bwc_base balancing utils bwc_mpfq ${pthread_libs} dl)
else(BUILD_DYNAMICALLY_LINKABLE_BWC)
target_link_libraries(bwc_base balancing utils bwc_mpfq ${pthread_libs})
endif(BUILD_DYNAMICALLY_LINKABLE_BWC)
ENDIF(CMAKE_SYSTEM_NAME MATCHES BSD)

set_target_properties(bwc_base PROPERTIES
        COMPILE_FLAGS "--mpi"
        LINK_FLAGS "--mpi")

target_link_libraries(bwc_mpfq ${gmp_libname} utils)
if (ENABLE_MPFQ_PRIME_FIELDS_FOR_DLOG)
    set_target_properties(bwc_mpfq PROPERTIES
        COMPILE_FLAGS "--mpi -DCOMPILE_MPFQ_PRIME_FIELDS")
else (ENABLE_MPFQ_PRIME_FIELDS_FOR_DLOG)
    set_target_properties(bwc_mpfq PROPERTIES
        COMPILE_FLAGS "--mpi")
endif (ENABLE_MPFQ_PRIME_FIELDS_FOR_DLOG)
set_target_properties(bwc_mpfq PROPERTIES
        LINK_FLAGS "--mpi")

add_library(matmul_common ${CADO_SHARED} matmul-common.c)
set_target_properties(matmul_common PROPERTIES
    COMPILE_FLAGS "--mpi"
    LINK_FLAGS "--mpi")

set(matmul_basic_sources    matmul-basic.c)
set(matmul_basicp_sources   matmul-basicp.c)
set(matmul_sliced_sources   matmul-sliced.cpp)
set(matmul_bucket_sources   matmul-bucket.cpp)
set(matmul_threaded_sources matmul-threaded.c worker-threads.c)

macro(CONFIGURE_MATMUL_LIB v comma i)
    string(COMPARE EQUAL "${v}" "u64k1" doing_u64)
    string(COMPARE EQUAL "${i}" "bucket" doing_bucket)
    set(sources matmul_facade.c ${matmul_${i}_sources})
    if(doing_u64 AND doing_bucket)
        if(HAVE_GAS_SYNTAX_ASSEMBLY_SOURCES)
            set(sources ${sources}
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
        endif(HAVE_GAS_SYNTAX_ASSEMBLY_SOURCES)
    endif(doing_u64 AND doing_bucket)

    add_library(matmul_${v}_${i} ${CADO_SHARED} ${sources})
    target_link_libraries(matmul_${v}_${i} matmul_common bwc_mpfq)
    set_target_properties(matmul_${v}_${i} PROPERTIES
        COMPILE_FLAGS "--mpi -DSELECT_ABASE_${v} -DMM_IMPL=${i} -DMM_ABASE=${v}"
            LINK_FLAGS "--mpi")
    set(all_matmul_libs ${all_matmul_libs} "matmul_${v}_${i}")
    install(TARGETS matmul_${v}_${i} DESTINATION "lib/${package_basename}")
endmacro(CONFIGURE_MATMUL_LIB)

if(BUILD_DYNAMICALLY_LINKABLE_BWC)
set(MATMUL_LIBS_PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
set(MATMUL_LIBS_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
endif(BUILD_DYNAMICALLY_LINKABLE_BWC)

configure_file(matmul-libnames.h.in matmul-libnames.h ESCAPE_QUOTES @ONLY)
include_directories (${CMAKE_CURRENT_BINARY_DIR})


# This list must match the one in linalg/bwc/matmul.c
CONFIGURE_MATMUL_LIB(u64k1     , bucket)
CONFIGURE_MATMUL_LIB(u64k2     , bucket)
CONFIGURE_MATMUL_LIB(u64k1     , basic)
CONFIGURE_MATMUL_LIB(u64k2     , basic)
CONFIGURE_MATMUL_LIB(u64k1     , sliced)
CONFIGURE_MATMUL_LIB(u64k2     , sliced)
CONFIGURE_MATMUL_LIB(u64k1     , threaded)
CONFIGURE_MATMUL_LIB(u64k2     , threaded)
if (ENABLE_MPFQ_PRIME_FIELDS_FOR_DLOG)
CONFIGURE_MATMUL_LIB(p_1       , basicp)
CONFIGURE_MATMUL_LIB(p_2       , basicp)
CONFIGURE_MATMUL_LIB(p_3       , basicp)
CONFIGURE_MATMUL_LIB(p_4       , basicp)
CONFIGURE_MATMUL_LIB(p_8       , basicp)
endif (ENABLE_MPFQ_PRIME_FIELDS_FOR_DLOG)

if(NOT BUILD_DYNAMICALLY_LINKABLE_BWC)
target_link_libraries(bwc_base ${all_matmul_libs})
endif(NOT BUILD_DYNAMICALLY_LINKABLE_BWC)

macro(bwc_program v)
    add_executable(${v} ${ARGN})
    install(TARGETS ${v} RUNTIME DESTINATION bin/linalg/bwc)
    set_target_properties(${v} PROPERTIES
            LINKER_LANGUAGE CXX
            COMPILE_FLAGS "--mpi"
            LINK_FLAGS "--mpi")
        # -Wl,-rpath,${CMAKE_CURRENT_BINARY_DIR}")
    target_link_libraries(${v} bwc_base)
endmacro(bwc_program)

add_library(matmul_mf matmul-mf.c)
    set_target_properties(matmul_mf PROPERTIES
            COMPILE_FLAGS "--mpi"
            LINK_FLAGS "--mpi")
target_link_libraries(matmul_mf mf)
bwc_program(build build.c)
target_link_libraries(build matmul_mf)
bwc_program(bench bench.c worker-threads.c)
target_link_libraries(bench matmul_mf)
set_source_files_properties (../gauss.c COMPILE_FLAGS "-DMULTI_ROW=3 -DNO_MAIN")
bwc_program(prep prep.c xymats.c ../gauss.c)
bwc_program(krylov krylov.c xymats.c xdotprod.c)
bwc_program(dispatch dispatch.c xymats.c xdotprod.c)
bwc_program(mksol mksol.c xymats.c xdotprod.c)
bwc_program(gather gather.c)
bwc_program(secure secure.c)
bwc_program(spmv_test spmv_test.c)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl DESTINATION bin/linalg/bwc)

# Also lingen.

link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
link_directories (${CADO_NFS_BINARY_DIR}/gf2x/gf2x-fft/.libs)

include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/gf2x-fft)
add_executable(lingen lingen_binary.cpp)
install(TARGETS lingen RUNTIME DESTINATION bin/linalg/bwc)
target_link_libraries(lingen bwc_base utils gf2x-fft gf2x cantor128)
add_dependencies(lingen gf2x-build)

set_target_properties(lingen PROPERTIES
        COMPILE_FLAGS "--mpi"
        LINK_FLAGS "--mpi")

add_executable(short_matmul ${EXCLUDE_FROM_RELEASE_BUILD} short_matmul.c)
target_link_libraries(short_matmul balancing)

add_executable(cleanup cleanup.c ../blockmatrix.c ../bit_matrices.c
    ../gauss.c)
target_link_libraries(cleanup utils m)

install(TARGETS cleanup RUNTIME DESTINATION bin/linalg/bwc)

set(plingen_prime_nwords 8)

if(WITH_MPIR)
    file(GLOB flint_fft_files flint-fft/*.c)
    add_library(flint-fft ${flint_fft_files})
    add_executable(test-flint test-flint.c)
    target_link_libraries(test-flint flint-fft ${gmp_libname})
endif(WITH_MPIR)


if (ENABLE_MPFQ_PRIME_FIELDS_FOR_DLOG)
add_executable(plingen EXCLUDE_FROM_ALL
                bw-common.c
                bw-common-mpi.c
                plingen.c
                plingen-tuning.cpp
                lingen-matpoly.c
                lingen-bigmatpoly.c
                lingen-polymat.c
                lingen-bigpolymat.c
                mpfq/abase_p_${plingen_prime_nwords}.c
        )
set_source_files_properties(plingen-tuning.cpp COMPILE_FLAGS "-std=c++11")

target_link_libraries(plingen ${gmp_libname} utils)
set_target_properties(plingen PROPERTIES
    COMPILE_FLAGS "--mpi -DSELECT_ABASE_p_${plingen_prime_nwords}"
    LINK_FLAGS "--mpi"
    )
# install(TARGETS plingen RUNTIME DESTINATION bin/linalg/bwc)
endif (ENABLE_MPFQ_PRIME_FIELDS_FOR_DLOG)

include(CMakeLists-nodist.txt OPTIONAL)

