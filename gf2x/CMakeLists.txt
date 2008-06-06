
##########################################################################
# Some trivialities

include(CheckTypeSize)
check_type_size("unsigned long int" ulong_size)
math (EXPR w "8*${ulong_size}")
message("-- Size of ulong is ${w} bits")

set(o ${CMAKE_CURRENT_BINARY_DIR})
set(s ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -W -Wall -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -W -Wall")

macro(mycopy src dst)
    add_custom_command(OUTPUT ${dst}
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy ${src} ${dst}
        DEPENDS ${src})
endmacro(mycopy)

macro(copy_from_source src dst)
mycopy(${s}/${src} ${o}/${dst})
endmacro(copy_from_source)

macro(add_cflags f)
    get_source_file_property(tmpvar ${f} COMPILE_FLAGS)
    if(${tmpvar} STREQUAL NOTFOUND)
        set(tmpvar)
    endif(${tmpvar} STREQUAL NOTFOUND)
    foreach(v ${ARGN})
        set(tmpvar "${tmpvar} ${v}")
    endforeach(v)
    set_source_files_properties(${f} COMPILE_FLAGS ${tmpvar})
endmacro(add_cflags)

##########################################################################

add_executable(gen_bb_mul_code EXCLUDE_FROM_ALL gen_bb_mul_code.c)
get_target_property(gen_bb_mul_code genbb LOCATION)

include_directories(${s})
include_directories(${o})

set(small_targets)
set(small_sizes)
macro(handle_small_size k)
    add_executable(tuneup-${k} EXCLUDE_FROM_ALL ${o}/tuneup-${k}.c)
    copy_from_source(tuneup.c tuneup-${k}.c)
    add_cflags(${o}/tuneup-${k}.c -DTUNE_MUL1)
    add_cflags(${o}/tuneup-${k}.c -DINLINES_FILE=\\\\"mul-inlines-k${k}.c\\\\")
    set_source_files_properties(${o}/tuneup-${k}.c OBJECT_DEPENDS mul-inlines-k${k}.c)
    add_custom_command(
        OUTPUT ${o}/mul-inlines-k${k}.c
        DEPENDS gen_bb_mul_code
        COMMAND ${o}/gen_bb_mul_code ${w} ${k} > ${o}/mul-inlines-k${k}.c
        )
    set(small_targets ${small_targets} tuneup-${k})
    set(small_sizes ${small_sizes} ${k})
endmacro(handle_small_size)

foreach(v RANGE 1 8)
    handle_small_size (${v})
endforeach(v)

add_custom_command(
    OUTPUT ${o}/mul-inlines.c
    DEPENDS ${small_targets}
    COMMAND
    ${s}/wizard_mul1.sh ${o} tuneup- ${small_sizes}
)

# useful ?
set_source_files_properties(${o}/mul-inlines.c GENERATED)

add_custom_target(tune1 DEPENDS ${o}/mul-inlines.c)

set_source_files_properties(mul-toom.c OBJECT_DEPENDS ${o}/mul-inlines.c)

set(libgf2x_files mul-toom.c mul-gf2x.c mulfft-bit.c)

add_executable(tunetoom ${libgf2x_files} tunetoom.c timing.c replace.c)
add_executable(tunefft ${libgf2x_files} tunefft.c timing.c replace.c)
add_executable(test ${libgf2x_files} test.c timing.c)
add_library(gf2x ${libgf2x_files})

add_executable(factor factor.cpp)
target_link_libraries(factor gf2x ntl gmp)

add_custom_target(tune-all
    DEPENDS tune1 do-tunetoom do-tunefft
    EXCLUDE_FROM_ALL
)
add_custom_target(do-tunetoom DEPENDS tunetoom COMMAND
    ${o}/tunetoom 2048)
add_custom_target(do-tunefft DEPENDS tunefft COMMAND
    ${o}/tunefft 100000)
add_dependencies(do-tunefft do-tunetoom)




