# Some compilation defaults specific to this directory

# The following is needed for some SSE stuff
# set (CMAKE_C_FLAGS   "-flax-vector-conversions ${CMAKE_C_FLAGS}")

include_directories (ecm)
include_directories (${CMAKE_CURRENT_BINARY_DIR})
link_directories (${CMAKE_CURRENT_BINARY_DIR}/ecm)

add_executable (makefb makefb.c)
target_link_libraries (makefb utils)

set (las_sources las.c fb.c
    las-coordinates.c las-debug.c las-report-stats.c las-duplicate.c
    las-config.c bucket.c las-smallsieve.c las-unsieve.c las-arith.c
    las-qlattice.c las-descent-helpers.cpp las-cofactor.c)
set (las_libraries las-norms trialdiv facul utils ${pthread_libs})

add_executable (las ${las_sources})
target_link_libraries (las ${las_libraries})

add_executable (las_tracek EXCLUDE_FROM_ALL ${las_sources})
set_property(TARGET las_tracek PROPERTY COMPILE_DEFINITIONS TRACE_K=1)
target_link_libraries (las_tracek ${las_libraries})

# Targets with the minimum acceptable thresholds for K-buckets and
# M-buckets, so that this code can be tested without requiring an
# exorbitant amount of memory
add_executable (las_kbucket EXCLUDE_FROM_ALL ${las_sources})
set_property(TARGET las_kbucket PROPERTY COMPILE_DEFINITIONS THRESHOLD_K_BUCKETS=16)
target_link_libraries (las_kbucket ${las_libraries})
add_executable (las_mbucket EXCLUDE_FROM_ALL ${las_sources})
set_property(TARGET las_mbucket PROPERTY COMPILE_DEFINITIONS THRESHOLD_K_BUCKETS=16 THRESHOLD_M_BUCKETS=64)
target_link_libraries (las_mbucket ${las_libraries})

add_library (trialdiv trialdiv.c)

add_library (las-norms las-norms.c)

add_executable(freerel freerel.c fb.c)
target_link_libraries(freerel utils)

add_subdirectory(ecm)

add_executable(dupsup dupsup.c las-duplicate.c las-norms.c las-qlattice.c las-coordinates.c las-cofactor.c fb.c)
target_link_libraries(dupsup utils_with_io facul)

install(TARGETS makefb las freerel RUNTIME DESTINATION ${BINSUFFIX}/sieve)

include(CMakeLists-nodist.txt OPTIONAL)
