# Some compilation defaults specific to this directory

# The following is needed for some SSE stuff
# set (CMAKE_C_FLAGS   "-flax-vector-conversions ${CMAKE_C_FLAGS}")

include_directories (ecm)
#include_directories (${CMAKE_CURRENT_BINARY_DIR})
link_directories (${CMAKE_CURRENT_BINARY_DIR}/ecm)

# ugly (copied from linalg/bwc/CMakeLists.txt...)
add_executable (makefb makefb.c)
target_link_libraries (makefb utils)

set (las_sources las.cpp fb.cpp las-plattice.cpp
  las-coordinates.cpp las-debug.cpp las-report-stats.cpp las-duplicate.cpp
  las-config.c bucket.cpp las-smallsieve.cpp las-unsieve.cpp
  las-unsieve-sse2.cpp las-arith.cpp
  las-trialdiv.cpp las-fb.cpp
  las-qlattice.cpp las-descent-trees.cpp las-cofactor.cpp
  las-fill-in-buckets.cpp las-threads.cpp)
set (las_libraries las-norms trialdiv facul utils ${pthread_libs})
# set (las_libraries trialdiv facul utils ${pthread_libs})

add_executable (las ${las_sources})
target_link_libraries (las ${las_libraries})
set_target_properties(las PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")

add_executable (las_tracek EXCLUDE_FROM_ALL ${las_sources})
set_property(TARGET las_tracek PROPERTY COMPILE_DEFINITIONS TRACE_K=1)
set_target_properties(las_tracek PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")
target_link_libraries (las_tracek ${las_libraries})

add_executable (fb_test EXCLUDE_FROM_ALL fb.cpp)
set_property(TARGET fb_test PROPERTY COMPILE_DEFINITIONS TESTDRIVE=1)
target_link_libraries (fb_test utils)

add_executable (las_descent ${las_sources} las-dlog-base.cpp)
set_property(TARGET las_descent PROPERTY COMPILE_DEFINITIONS DLP_DESCENT SUPPORT_LARGE_Q)
set_target_properties(las_descent PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")
target_link_libraries (las_descent ${las_libraries})

add_library (trialdiv trialdiv.c)
add_library (las-norms las-norms.cpp)

add_executable(freerel freerel.c fb.cpp las-plattice.cpp)
target_link_libraries(freerel utils)

add_subdirectory(ecm)
add_subdirectory(strategies)

add_executable(dupsup dupsup.cpp las-duplicate.cpp las-qlattice.cpp
    las-plattice.cpp las-coordinates.cpp las-cofactor.cpp fb.cpp
    las-norms.cpp
    las-unsieve.cpp
    las-unsieve-sse2.cpp
    las-smallsieve.cpp
    las-trialdiv.cpp las-fb.cpp
    trialdiv.c
    )
target_link_libraries(dupsup utils_with_io facul las-norms)

add_executable(fake_rels EXCLUDE_FROM_ALL fake_rels.cpp)
target_link_libraries(fake_rels utils_with_io)


install(TARGETS makefb las las_descent freerel RUNTIME DESTINATION ${LIBSUFFIX}/sieve) 

include(CMakeLists-nodist.txt OPTIONAL)
