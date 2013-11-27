# Some compilation defaults specific to this directory

# The following is needed for some SSE stuff
# set (CMAKE_C_FLAGS   "-flax-vector-conversions ${CMAKE_C_FLAGS}")

include_directories (ecm)
include_directories (${CMAKE_CURRENT_BINARY_DIR})
link_directories (${CMAKE_CURRENT_BINARY_DIR}/ecm)

add_executable (makefb makefb.c)
target_link_libraries (makefb utils)

add_executable (las las.c fb.c trialdiv.c
    las-coordinates.c las-debug.c las-report-stats.c las-duplicate.c
    las-config.c las-norms.c bucket.c las-smallsieve.c las-unsieve.c las-arith.c
    las-qlattice.c las-descent-helpers.cpp)
target_link_libraries (las facul utils ${pthread_libs})

add_executable(test_trialdiv EXCLUDE_FROM_ALL test_trialdiv.c trialdiv.c)
target_link_libraries (test_trialdiv utils)

add_executable(freerel freerel.c fb.c)
target_link_libraries(freerel utils)

add_subdirectory(ecm)

install(TARGETS makefb las freerel RUNTIME DESTINATION bin/sieve)

include(CMakeLists-nodist.txt OPTIONAL)
