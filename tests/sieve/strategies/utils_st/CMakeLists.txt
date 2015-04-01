set (SUF_SRC ${CADO_NFS_SOURCE_DIR}/sieve/strategies/utils_st/)

include_directories(${CADO_NFS_SOURCE_DIR}/sieve/strategies/utils_st)
include_directories(${CADO_NFS_SOURCE_DIR})
include_directories(${CADO_NFS_SOURCE_DIR}/sieve/ecm)
include_directories(${CADO_NFS_SOURCE_DIR}/sieve)

add_executable(test_fm test_fm.c ${SUF_SRC}fm.c ${SUF_SRC}tab_fm.c) 
run_test (test_fm test_fm)

add_executable(test_strategy
  test_strategy.c
  ${SUF_SRC}tab_strategy.c
  ${SUF_SRC}fm.c
  ${SUF_SRC}strategy.c
  ${SUF_SRC}fm.c
  ${SUF_SRC}tab_fm.c)

run_test (test_strategy test_strategy)

add_executable(test_decomp
  test_decomp.c
  ${SUF_SRC}tab_decomp.c
  ${SUF_SRC}decomp.c)

if(MINGW)
    target_link_libraries(test_decomp gnurx)
endif(MINGW)

run_test (test_decomp test_decomp)

add_executable(test_convex_hull
  test_convex_hull.c
  ${SUF_SRC}convex_hull.c
  ${SUF_SRC}point.c
  ${SUF_SRC}tab_point.c)

run_test (test_convex_hull1 test_convex_hull)
run_test (test_convex_hull2 test_convex_hull 1000 0.001)
run_test (test_convex_hull3 test_convex_hull 100000 0.00001)
run_test (test_convex_hull4 test_convex_hull 100000 -0.00001)
