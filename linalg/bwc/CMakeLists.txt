
include_directories(..)

Include(FindMPI)

add_executable(test_hello EXCLUDE_FROM_ALL test_hello.c parallelizing_info.c)
add_executable(random random.c)
add_executable(split split.c bw-common.c)
add_executable(shell shell.c parallelizing_info.c bw-common.c bw-common-mpi.c)
add_executable(acollect acollect.c bw-common.c)

install(TARGETS random split shell acollect RUNTIME DESTINATION bin/linalg/bwc)

target_link_libraries(test_hello utils)
target_link_libraries(random utils m)
target_link_libraries(split utils)
target_link_libraries(shell utils)
target_link_libraries(acollect utils)

set(bwc_base_objs
    parallelizing_info.c
    intersections.cpp
    xvectors.c
    bw-common.c
    bw-common-mpi.c
    rusage.c
    async.c
    debug.c)

add_library(bwc_base ${bwc_base_objs})

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
        message(STATUS "Trying to compile a sample .S file")
        TRY_COMPILE(asm_compiles ${CADO_NFS_BINARY_DIR}/config ${CADO_NFS_SOURCE_DIR}/config/asm_sample.c
            COMPILE_DEFINITIONS "-x assembler"
            OUTPUT_VARIABLE asm_compile_output)
        if(asm_compiles)
            message(STATUS "Trying to compile a sample .S file -- Success")
            set(${v}_matmul_sources ${${v}_matmul_sources} matmul-sliced-asm.S)
            set_property(SOURCE matmul-sliced-asm.S PROPERTY LANGUAGE C)
        else(asm_compiles)
            message(STATUS "Trying to compile a sample .S file -- Failed")
            message(STATUS ${asm_compile_output})
        endif(asm_compiles)
    endif(doing_u64)

    add_library(${v}_matmul ${${v}_matmul_sources})

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

install(PROGRAMS bwc.pl DESTINATION bin/linalg/bwc)

# Also lingen.

include_directories(${CADO_NFS_SOURCE_DIR}/gf2x)
include_directories(${CADO_NFS_BINARY_DIR}/gf2x)
link_directories(${CADO_NFS_BINARY_DIR}/gf2x/.libs)
include_directories(${CADO_NFS_SOURCE_DIR}/cantor)
add_executable(lingen fake_fft.c lingen_binary.cpp)
install(TARGETS lingen RUNTIME DESTINATION bin/linalg/bwc)
target_link_libraries(lingen bwc_base gf2x cantor128 utils)

