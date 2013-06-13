# the utils library is already linked with -lgmp and -lm
add_executable(dup1 dup1.c)
target_link_libraries(dup1 utils)

add_executable(dup2 dup2.c filter_utils.c filter_io.c filter_memalloc.c)
target_link_libraries(dup2 utils)
set_target_properties(dup2 PROPERTIES COMPILE_FLAGS "-DFOR_DUP2")
add_dependencies(dup2 antebuffer)

add_executable(purge purge.c filter_utils.c filter_io.c filter_memalloc.c)
target_link_libraries(purge utils pthread)
add_dependencies(purge antebuffer)

add_executable(merge merge.c markowitz.c mst.c report.c merge_mono.c filter_matrix.c sparse.c filter_utils.c filter_io.c filter_memalloc.c)
target_link_libraries(merge utils)

add_executable(replay replay.c sparse.c)
target_link_libraries(replay utils)

install(TARGETS dup1 dup2 purge merge replay RUNTIME DESTINATION bin/filter)

include(CMakeLists-nodist.txt OPTIONAL)
