include_directories (utils_st)
include_directories (${CMAKE_CURRENT_BINARY_DIR})

set (SUF_UTILS utils_st/)
set (SRC_UTILS ${SUF_UTILS}point.c ${SUF_UTILS}tab_point.c
    ${SUF_UTILS}convex_hull.c ${SUF_UTILS}fm.c ${SUF_UTILS}tab_fm.c
    ${SUF_UTILS}strategy.c ${SUF_UTILS}tab_strategy.c)

set (LIB_EXEC m facul utils)

add_executable(gfm gfm_exec.c generate_factoring_method.c ${SRC_UTILS})
target_link_libraries(gfm ${LIB_EXEC})

add_executable(benchfm benchfm_exec.c generate_factoring_method.c ${SRC_UTILS})
target_link_libraries(benchfm ${LIB_EXEC})

add_executable(gst 
  gst_exec.c 
  generate_strategies.c 
  ${SUF_UTILS}/decomp.c
  ${SUF_UTILS}/tab_decomp.c
  ${SRC_UTILS})

target_link_libraries(gst ${LIB_EXEC})

add_executable(finalst 
  final_exec.c
  finding_good_strategy.c 
  ${SRC_UTILS})

target_link_libraries(finalst ${LIB_EXEC})

add_executable(gen_decomp gen_decomp.c)
target_link_libraries(gen_decomp m)
