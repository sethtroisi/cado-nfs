
# Some compilation defaults specific to this directory
# The following flags are specific to gcc, and fix a regression in gcc 4.3.0

set_property(DIRECTORY APPEND PROPERTY COMPILE_FLAGS fno-ivopts -fno-regmove -fno-tree-sra -fno-tree-ter)


# CFLAGS=-O0 -g
# CFLAGS+=-DPROFILE

# if you compile with GCC 4.2 or later, you can also try -march=native
# instead of the default -mtune=generic (cf man gcc and
# http://gcc.gnu.org/ml/gcc-help/2008-07/msg00060.html)

# using -DSSE_NORM_INIT seems to be slower (about 3% with gcc 4.3.0)
# CFLAGS+=-DSSE_NORM_INIT
# -fkeep-inline-functions does not yield any speedup, is it useful?
# CFLAGS+=-fkeep-inline-functions

include_directories (ecm)
include_directories (${CMAKE_CURRENT_BINARY_DIR})
link_directories (ecm)

# XXX TODO: sieve_config.h is dead -- remove it, place relevant #defines in
# the source file instead. Then remove the custom command below as well as
# references to sieve_config.h here and in the companion file
# CMakeLists-nodist.txt

add_custom_command(OUTPUT sieve_config.h
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/sieve_config_h.in ${CMAKE_CURRENT_BINARY_DIR}/sieve_config.h)

set_source_files_properties(sieve_config.h PROPERTIES GENERATED true)

add_executable (makefb makefb.c)
target_link_libraries (makefb utils m)

add_executable (las las.c fb.c trialdiv.c sieve_config.h)
target_link_libraries (las facul utils m pthread)

# Compiling stuff for las-mt: uncomment those two lines for testing
# add_executable (las-mt las-mt.c fb.c trialdiv.c sieve_config.h)
# target_link_libraries (las-mt facul utils m pthread)

add_executable(trialdiv EXCLUDE_FROM_ALL trialdiv.c)
target_link_libraries (trialdiv utils)
set_property(TARGET trialdiv APPEND PROPERTY COMPILE_FLAGS -DTESTDRIVE)

add_executable(freerel freerel.c fb.c sieve_config.h)
target_link_libraries(freerel utils m)

add_subdirectory(ecm)

install(TARGETS makefb las freerel RUNTIME DESTINATION bin/sieve)

include(CMakeLists-nodist.txt OPTIONAL)

