cmake_minimum_required(VERSION 3.0)

include_directories (utils_st)
include_directories (${CMAKE_CURRENT_BINARY_DIR})


set (SUF_UTILS utils_st/)
set (SRC_UTILS ${SUF_UTILS}point.c ${SUF_UTILS}tab_point.c ${SUF_UTILS}convex_hull.c 
  ${SUF_UTILS}fm.c ${SUF_UTILS}tab_fm.c ${SUF_UTILS}strategy.c ${SUF_UTILS}tab_strategy.c)

set (LIB_EXEC m gmp facul utils)

set(EXEC_NAME1 gfm)
set(SRC_1 gfm_exec.c generate_factoring_method.c ${SRC_UTILS})
add_executable (${EXEC_NAME1} ${SRC_1})
target_link_libraries(${EXEC_NAME1} ${LIB_EXEC})

set(EXEC_NAME2 benchfm)
set(SRC_2 benchfm_exec.c generate_factoring_method.c ${SRC_UTILS})
add_executable (${EXEC_NAME2} ${SRC_2})
target_link_libraries(${EXEC_NAME2} ${LIB_EXEC})

