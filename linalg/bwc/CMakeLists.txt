
include_directories(. ..)
# include_directories(mpfq)

# Use with care. Normally, binaries not needed MPI for real should not
# emit MPI symbols since they don't do anything depending on this config
# variable. However this mechanism is somewhat fragile.
if (${WITH_MPI})
add_definitions(-DWITH_MPI)
endif (${WITH_MPI})

### bwc-specific configuration.

# unset(BUILD_DYNAMICALLY_LINKABLE_BWC)
# Dynamically linkable bwc needs shared libs.
set(BUILD_DYNAMICALLY_LINKABLE_BWC ${ENABLE_SHARED})
# having both shared and static is possible but messy. One needs two
# add_library calls, and then I'm not sure whether
# target_link_libraries picks the static or shared ones...

# This is a default setting
set(BWC_GF2_ARITHMETIC_BACKENDS u64k1 u64k2
    CACHE STRING
    "which low-level arithmetic backends to use for GF(2) linear algebra")
set(BWC_GF2_MATMUL_BACKENDS bucket basic sliced threaded
    CACHE STRING
    "which low-level matrix multiplication backends to use for GF(2) linear algebra")
set(BWC_GFP_ARITHMETIC_BACKENDS p_1 p_2 p_3 p_4 p_8 pz
    CACHE STRING
    "which GF(p) backends to compile for DLOG mode")
set(BWC_GFP_MATMUL_BACKENDS basicp
    CACHE STRING
    "which low-level matrix multiplication backends to use for GF(2) linear algebra")

# In case there's a user override, just make mention of the situation.
if(NOT BWC_GFP_ARITHMETIC_BACKENDS)
    message(STATUS "DL mode disabled for linear algebra")
else(NOT BWC_GFP_ARITHMETIC_BACKENDS)
    if(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
        message(STATUS "GF(p) implementations enabled for DL linear algebra: ${BWC_GFP_ARITHMETIC_BACKENDS}")
    else(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
        # The code in linalg/bwc/mpfq/ is 64-bit only, and relies on GNU
        # CC idioms. As it is not used in production yet, we don't want
        # to expose it on 32-bit platforms.
        set(BWC_GFP_ARITHMETIC_BACKENDS CACHE STRING
                "which GF(p) backends to compile for DLOG mode" FORCE)
        message(STATUS "Warning: DL mode forcibly disabled for linear algebra (amd64-only)")
    endif(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
endif(NOT BWC_GFP_ARITHMETIC_BACKENDS)

set(COOKED_BWC_BACKENDS)
foreach(v ${BWC_GF2_ARITHMETIC_BACKENDS})
    foreach(w ${BWC_GF2_MATMUL_BACKENDS})
        set(COOKED_BWC_BACKENDS ${COOKED_BWC_BACKENDS} "DO(${v}, ${w})")
    endforeach(w)
endforeach(v)
foreach(v ${BWC_GFP_ARITHMETIC_BACKENDS})
    foreach(w ${BWC_GFP_MATMUL_BACKENDS})
        set(COOKED_BWC_BACKENDS ${COOKED_BWC_BACKENDS} "DO(${v}, ${w})")
    endforeach(w)
endforeach(v)

configure_file(bwc_config_h.in bwc_config.h ESCAPE_QUOTES @ONLY)

### end bwc-specific configuration.


### bwc_base

# Note that because the balancing code contains balancing_get_matrix_u32,
# which calls the parallelizing_info layer, it's difficult to make them
# two distinct libraries. So far it's the only stopper though.
set(balancing_sources
    balancing.c
    balancing_workhorse.c
    balancing_file_source.c
    balancing_mpi_source.c)

if (HAVE_CURL)
set(balancing_sources ${balancing_sources} balancing_curl_source.c)
endif (HAVE_CURL)

set(bwc_base_sources
    ${balancing_sources}
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

add_library(bwc_base ${CADO_SHARED} ${bwc_base_sources})

if (HAVE_CURL)
target_link_libraries(bwc_base curl)
endif (HAVE_CURL)


set_target_properties(bwc_base PROPERTIES
        COMPILE_FLAGS "--mpi"
        LINK_FLAGS "--mpi")

if(ENABLE_SHARED)
    install(TARGETS bwc_base DESTINATION ${LIBSUFFIX})
endif(ENABLE_SHARED)

### end bwc_base


add_library(mf ${CADO_SHARED} mf.c)
target_link_libraries(mf bwc_base utils)
if(ENABLE_SHARED)
    install(TARGETS mf DESTINATION ${LIBSUFFIX})
endif(ENABLE_SHARED)
add_executable(mf_scan ${EXCLUDE_FROM_RELEASE_BUILD} mf_scan.c)
target_link_libraries(mf_scan mf)


include_directories(${CADO_NFS_SOURCE_DIR}/linalg)
add_executable(mf_bal mf_bal.c ../rowset_heap.cpp)
target_link_libraries(mf_bal mf m)
install(TARGETS mf_bal RUNTIME DESTINATION ${BINSUFFIX}/linalg/bwc)


configure_file(bwc.pl ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl ESCAPE_QUOTES @ONLY)

add_executable(random_matrix ${EXCLUDE_FROM_RELEASE_BUILD} random_matrix.c)

# balancing.c itself is in fact standalone, and lightweight. We prefer to
# avoid having the whole thing entangled with too many dependencies all
# around, since all that split.c actually needs is very, very little.
add_executable(split split.c bw-common.c balancing.c)

add_executable(acollect acollect.c bw-common.c)

install(TARGETS split acollect RUNTIME DESTINATION ${BINSUFFIX}/linalg/bwc)

target_link_libraries(random_matrix utils m)

target_link_libraries(split utils)
target_link_libraries(acollect utils)

add_executable(bcast-file ${EXCLUDE_FROM_RELEASE_BUILD} bcast-file.c)
set_target_properties(bcast-file PROPERTIES
        COMPILE_FLAGS "--mpi"
        LINK_FLAGS "--mpi")

### bwc_mpfq

set(bwc_abases u64k1 u64k2)

set(bwc_mpfq_sources mpfq/abase_vbase.c)
foreach(v ${bwc_abases})
    set(bwc_mpfq_sources ${bwc_mpfq_sources} mpfq/abase_${v}.c mpfq/abase_${v}_t.c )
endforeach(v)

foreach(v ${BWC_GFP_ARITHMETIC_BACKENDS})
    set(bwc_mpfq_sources ${bwc_mpfq_sources} mpfq/abase_${v}.c mpfq/abase_${v}_t.c )
endforeach(v)

add_library(bwc_mpfq ${CADO_SHARED} ${bwc_mpfq_sources})

set(bwc_mpfq_flags)
foreach(v ${BWC_GFP_ARITHMETIC_BACKENDS})
    set(bwc_mpfq_flags "${bwc_mpfq_flags} -DCOMPILE_MPFQ_PRIME_FIELD_${v}")
endforeach(v)

if(ENABLE_SHARED)
    install(TARGETS bwc_mpfq DESTINATION ${LIBSUFFIX})
endif(ENABLE_SHARED)

target_link_libraries(bwc_mpfq ${gmp_libname} utils)

set_target_properties(bwc_mpfq PROPERTIES COMPILE_FLAGS "--mpi ${bwc_mpfq_flags}")
set_target_properties(bwc_mpfq PROPERTIES LINK_FLAGS "--mpi")

### end bwc_mpfq


### link dependencies

IF(UNIX)
FIND_PROGRAM(CMAKE_UNAME uname /bin /usr/bin /usr/local/bin )
IF(CMAKE_UNAME)
EXEC_PROGRAM(uname ARGS -s OUTPUT_VARIABLE CMAKE_SYSTEM_NAME)
ENDIF(CMAKE_UNAME)
ENDIF(UNIX)
IF(CMAKE_SYSTEM_NAME MATCHES BSD)
target_link_libraries(bwc_base utils bwc_mpfq ${pthread_libs})
ELSE(CMAKE_SYSTEM_NAME MATCHES BSD)
if(BUILD_DYNAMICALLY_LINKABLE_BWC)
target_link_libraries(bwc_base utils bwc_mpfq ${pthread_libs} dl)
else(BUILD_DYNAMICALLY_LINKABLE_BWC)
target_link_libraries(bwc_base utils bwc_mpfq ${pthread_libs})
endif(BUILD_DYNAMICALLY_LINKABLE_BWC)
ENDIF(CMAKE_SYSTEM_NAME MATCHES BSD)

# matmul_common is a dependency of all the impl-specific libs. When doing
# dynamic linking, there's no real need to make matmul_common *also*
# shared, as embarking the code within the impl-specific .so is fine.
# Unfortunately, we do want matmul_common to be compiled with -fPIC,
# then. Which is easiest to trigger by making it a shared object.
add_library(matmul_common ${CADO_SHARED} matmul-common.c)
set_target_properties(matmul_common PROPERTIES
    COMPILE_FLAGS "--mpi"
    LINK_FLAGS "--mpi")
if(ENABLE_SHARED)
    install(TARGETS matmul_common DESTINATION ${LIBSUFFIX})
endif(ENABLE_SHARED)


set(matmul_basic_sources    matmul-basic.c)
set(matmul_basicp_sources   matmul-basicp.c)
set(matmul_sliced_sources   matmul-sliced.cpp)
set(matmul_bucket_sources   matmul-bucket.cpp)
set(matmul_threaded_sources matmul-threaded.c worker-threads.c)

macro(CONFIGURE_MATMUL_LIB v i)
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
    if(ENABLE_SHARED)
        install(TARGETS matmul_${v}_${i} DESTINATION ${LIBSUFFIX})
    endif(ENABLE_SHARED)
endmacro(CONFIGURE_MATMUL_LIB)

foreach(v ${BWC_GF2_ARITHMETIC_BACKENDS})
    foreach(w ${BWC_GF2_MATMUL_BACKENDS})
        CONFIGURE_MATMUL_LIB(${v} ${w})
    endforeach(w)
endforeach(v)
foreach(v ${BWC_GFP_ARITHMETIC_BACKENDS})
    foreach(w ${BWC_GFP_MATMUL_BACKENDS})
        CONFIGURE_MATMUL_LIB(${v} ${w})
    endforeach(w)
endforeach(v)


if(BUILD_DYNAMICALLY_LINKABLE_BWC)
set(MATMUL_LIBS_PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
set(MATMUL_LIBS_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
endif(BUILD_DYNAMICALLY_LINKABLE_BWC)

configure_file(matmul-libnames.h.in matmul-libnames.h ESCAPE_QUOTES @ONLY)
include_directories (${CMAKE_CURRENT_BINARY_DIR})

if(NOT BUILD_DYNAMICALLY_LINKABLE_BWC)
target_link_libraries(bwc_base ${all_matmul_libs})
endif(NOT BUILD_DYNAMICALLY_LINKABLE_BWC)

### matmul_mf
add_library(matmul_mf ${CADO_SHARED} matmul-mf.c)
    set_target_properties(matmul_mf PROPERTIES
            COMPILE_FLAGS "--mpi"
            LINK_FLAGS "--mpi")
target_link_libraries(matmul_mf mf)
if(ENABLE_SHARED)
    install(TARGETS matmul_mf DESTINATION ${LIBSUFFIX})
endif(ENABLE_SHARED)
### end matmul_mf

macro(bwc_program v)
    add_executable(${v} ${ARGN})
    install(TARGETS ${v} RUNTIME DESTINATION ${BINSUFFIX}/linalg/bwc)
    set_target_properties(${v} PROPERTIES
            LINKER_LANGUAGE CXX
            COMPILE_FLAGS "--mpi"
            LINK_FLAGS "--mpi")
        # -Wl,-rpath,${CMAKE_CURRENT_BINARY_DIR}")
    target_link_libraries(${v} bwc_base)
endmacro(bwc_program)

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

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl DESTINATION ${BINSUFFIX}/linalg/bwc)

# Also lingen.

link_directories (${CADO_NFS_BINARY_DIR}/gf2x/.libs)
link_directories (${CADO_NFS_BINARY_DIR}/gf2x/gf2x-fft/.libs)

include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/gf2x-fft)
add_executable(lingen lingen_binary.cpp)
install(TARGETS lingen RUNTIME DESTINATION ${BINSUFFIX}/linalg/bwc)
target_link_libraries(lingen bwc_base utils gf2x-fft gf2x cantor128)
add_dependencies(lingen gf2x-build)

set_target_properties(lingen PROPERTIES
        COMPILE_FLAGS "--mpi"
        LINK_FLAGS "--mpi")

add_executable(short_matmul ${EXCLUDE_FROM_RELEASE_BUILD} short_matmul.c)
target_link_libraries(short_matmul bwc_base)

add_executable(cleanup cleanup.c ../blockmatrix.c ../bit_matrices.c
    ../gauss.c)
target_link_libraries(cleanup utils m)

install(TARGETS cleanup RUNTIME DESTINATION ${BINSUFFIX}/linalg/bwc)

if(WITH_MPIR)
    # This is a convenience library which is used **ONLY** by one binary,
    # so there's absolutely no point in making it shared.
    file(GLOB flint_fft_files flint-fft/*.c)
    add_library(flint-fft ${flint_fft_files})
    add_executable(test-flint test-flint.c)
    target_link_libraries(test-flint flint-fft utils ${gmp_libname})
endif(WITH_MPIR)

foreach(v ${BWC_GFP_ARITHMETIC_BACKENDS})
    set(plingen_${v}_sources
        bw-common.c
        bw-common-mpi.c
        plingen.c
        plingen-tuning.cpp
        logline.c
        lingen-matpoly.c
        lingen-bigmatpoly.c
        lingen-polymat.c
        lingen-bigpolymat.c
        mpfq/abase_${v}.c
        )
    if(WITH_MPIR)
        set(plingen_${v}_sources
            ${plingen_${v}_sources}
            lingen-matpoly-ft.c
            lingen-bigmatpoly-ft.c)
    endif(WITH_MPIR)

    add_executable(plingen_${v} EXCLUDE_FROM_ALL ${plingen_${v}_sources})

    if(WITH_MPIR)
        target_link_libraries(plingen_${v} flint-fft ${gmp_libname} utils)
    else(WITH_MPIR)
        target_link_libraries(plingen_${v} ${gmp_libname} utils)
    endif(WITH_MPIR)
    set_target_properties(plingen_${v} PROPERTIES
        COMPILE_FLAGS "--mpi -DSELECT_ABASE_${v}"
        LINK_FLAGS "--mpi"
        )
endforeach(v)



add_custom_target(bwc_full)
add_dependencies(bwc_full mf_bal random_matrix split acollect lingen
    cleanup prep krylov dispatch mksol gather secure bwc_base
    all_matmul_libs)

