set (SUF_SRC ${CADO_NFS_SOURCE_DIR}/sieve/strategies)

set (SUF_UTILS ${SUF_SRC}/utils_st/)
set (SRC_UTILS ${SUF_UTILS}point.c ${SUF_UTILS}tab_point.c
  ${SUF_UTILS}convex_hull.c ${SUF_UTILS}fm.c ${SUF_UTILS}tab_fm.c
  ${SUF_UTILS}strategy.c ${SUF_UTILS}tab_strategy.c)

set (LIB_EXEC m facul utils)

include_directories(${CADO_NFS_SOURCE_DIR})
include_directories(${CADO_NFS_SOURCE_DIR}/sieve)
include_directories(${CADO_NFS_SOURCE_DIR}/sieve/ecm)
include_directories(${CADO_NFS_SOURCE_DIR}/sieve/strategies)
include_directories(${CADO_NFS_SOURCE_DIR}/sieve/strategies/utils_st)

add_executable(test_generate_composite_numbers
  test_generate_composite_numbers.c
  ${SUF_SRC}/generate_factoring_method.c
  ${SRC_UTILS})
target_link_libraries(test_generate_composite_numbers ${LIB_EXEC})
run_test (test_generate_composite_numbers test_generate_composite_numbers)


add_subdirectory (utils_st)