add_executable(duplicates duplicates.c files.c)
target_link_libraries(duplicates utils)

add_executable(purge purge.c files.c)
target_link_libraries(purge utils m)

add_executable(merge merge.c markowitz.c swar.c mst.c report.c dclist.c merge_mono.c sparse_mat.c sparse.c)
target_link_libraries(merge utils m)

add_executable(postprocess postprocess.c markowitz.c mst.c report.c merge_mono.c sparse_mat.c sparse.c)
set_property(TARGET postprocess APPEND PROPERTY COMPILE_FLAGS -DUSE_MARKOWITZ)
target_link_libraries(postprocess utils m)

add_executable(replay replay.c sparse.c)
target_link_libraries(replay utils)

install(TARGETS duplicates purge merge replay RUNTIME DESTINATION bin/merge)
install(TARGETS postprocess RUNTIME DESTINATION bin/merge)
