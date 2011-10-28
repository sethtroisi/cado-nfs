# Some compilation defaults specific to this directory

# The following is needed for some SSE stuff
# set (CMAKE_C_FLAGS   "-flax-vector-conversions ${CMAKE_C_FLAGS}")

include_directories (ecm)
include_directories (${CMAKE_CURRENT_BINARY_DIR})
link_directories (${CMAKE_CURRENT_BINARY_DIR}/ecm)

add_executable (makefb makefb.c)
target_link_libraries (makefb utils m)

add_executable (las las.c fb.c trialdiv.c
    las-coordinates.c las-debug.c las-report-stats.c
    las-config.c las-norms.c bucket.c las-smallsieve.c las-unsieve.c las-arith.c
    las-qlattice.c)
target_link_libraries (las facul utils m ${pthread_libs})

add_executable(trialdiv EXCLUDE_FROM_ALL trialdiv.c)
target_link_libraries (trialdiv utils)
set_property(TARGET trialdiv APPEND PROPERTY COMPILE_FLAGS -DTESTDRIVE)

add_executable(freerel freerel.c fb.c)
target_link_libraries(freerel utils m)

add_subdirectory(ecm)

install(TARGETS makefb las freerel RUNTIME DESTINATION bin/sieve)

include(CMakeLists-nodist.txt OPTIONAL)
