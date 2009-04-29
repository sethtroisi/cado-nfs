
include_directories(..)

Include(FindMPI)

add_executable(test_hello test_hello.c parallelizing_info.c)
add_executable(random random.c)
add_executable(split split.c bw-common.c)
add_executable(shell shell.c parallelizing_info.c bw-common.c bw-common-mpi.c)
add_executable(acollect acollect.c bw-common.c)

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
    add_library(${v}_matmul
        matmul-common.c matmul-basic.c matmul-sliced.cpp
        matmul-threaded.c matmul-bucket.cpp matmul.c
        matmul-sliced-asm.S)

    set_property(SOURCE matmul-sliced-asm.S PROPERTY LANGUAGE C)

    push_abase_cflags(${v} matmul)
    target_link_libraries(${v}_matmul mat_io)

    add_library(${v}_bwc info_file.c matmul_top.c)
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
        target_link_libraries(${v}_${prg} ${v}_bwc)
        push_abase_cflags(${v} ${prg})
        set_target_properties(${v}_${prg} PROPERTIES LINKER_LANGUAGE CXX)
    endforeach(v)
endmacro(configure_program)

configure_program(krylov ${allvariants})
configure_program(mksol  ${allvariants})
configure_program(build  ${allvariants})
configure_program(bench  ${allvariants})
configure_program(prep   u64n)
configure_program(gather u64n)
configure_program(secure u64)
