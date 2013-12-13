# the utils library is already linked with -lgmp and -lm
add_executable(dup1 dup1.c)
target_link_libraries(dup1 filter_with_io)

add_executable(dup2 dup2.c)
target_link_libraries(dup2 filter_with_io)

add_executable(purge filter_utils.c purge.c)
target_link_libraries(purge filter_with_io)

add_executable(merge merge.c filter_utils.c markowitz.c mst.c report.c merge_mono.c filter_matrix.c sparse.c)
target_link_libraries(merge filter_with_io)

add_executable(merge-dl merge.c filter_utils.c markowitz.c mst.c report.c merge_mono.c filter_matrix.c sparse.c)
target_link_libraries(merge-dl filter_with_io)
set_target_properties(merge-dl PROPERTIES COMPILE_FLAGS "-DFOR_DL")

add_executable(replay filter_utils.c replay.c sparse.c)
target_link_libraries(replay filter_with_io)

add_executable(replay-dl filter_utils.c replay.c sparse.c)
target_link_libraries(replay-dl filter_with_io)
set_target_properties(replay-dl PROPERTIES COMPILE_FLAGS "-DFOR_DL")

add_executable(reconstructlog-dl reconstructlog.c)
target_link_libraries(reconstructlog-dl filter_with_io)

install(TARGETS dup1 dup2 purge merge replay RUNTIME DESTINATION bin/filter)

include(CMakeLists-nodist.txt OPTIONAL)
