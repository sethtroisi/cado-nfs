
include_directories(..)

if (${WITH_MPI})
        set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
        set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
        add_definitions(-DWITH_MPI)
endif (${WITH_MPI})

configure_file(bwc.pl ${CMAKE_CURRENT_BINARY_DIR}/bwc.pl ESCAPE_QUOTES @ONLY)

add_executable(test_hello EXCLUDE_FROM_ALL
    test_hello.c parallelizing_info.c)
add_executable(random random.c)
add_executable(split split.c bw-common.c)
add_executable(shell
    shell.c parallelizing_info.c bw-common.c bw-common-mpi.c)
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
