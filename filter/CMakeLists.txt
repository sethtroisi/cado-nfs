# the utils library is already linked with -lgmp and -lm
add_executable(dup1 dup1.c)
target_link_libraries(dup1 utils)

add_executable(dup2 dup2.c)
target_link_libraries(dup2 utils)

add_executable(purge purge.c)
target_link_libraries(purge utils pthread)

add_executable(purge-mt purge-mt.c)
target_link_libraries(purge-mt utils ${pthread_libs})

add_executable(merge merge.c markowitz.c mst.c report.c merge_mono.c filter_matrix.c sparse.c)
target_link_libraries(merge utils)

add_executable(replay replay.c sparse.c)
target_link_libraries(replay utils)

install(TARGETS dup1 dup2 purge merge replay RUNTIME DESTINATION bin/filter)

include(CMakeLists-nodist.txt OPTIONAL)
