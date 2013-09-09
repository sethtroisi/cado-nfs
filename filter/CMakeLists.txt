set (sources
        filter_utils.c filter_io.cpp filter_memalloc.c
)

add_library (filterutils ${CADO_SHARED} ${sources})
target_link_libraries(filterutils utils pthread m)
install(TARGETS filterutils DESTINATION lib/${package_basename})
add_dependencies(filterutils antebuffer)

# the utils library is already linked with -lgmp and -lm
add_executable(dup1 dup1.c)
target_link_libraries(dup1 filterutils)

add_executable(dup2 dup2.c)
target_link_libraries(dup2 filterutils)

add_executable(purge purge.c)
target_link_libraries(purge filterutils)

add_executable(merge merge.c markowitz.c mst.c report.c merge_mono.c filter_matrix.c sparse.c)
target_link_libraries(merge filterutils)

add_executable(merge-dl merge.c markowitz.c mst.c report.c merge_mono.c filter_matrix.c sparse.c)
target_link_libraries(merge-dl filterutils)
set_target_properties(merge-dl PROPERTIES COMPILE_FLAGS "-DFOR_DL")

add_executable(replay replay.c sparse.c)
target_link_libraries(replay utils)

add_executable(replay-dl replay.c sparse.c)
target_link_libraries(replay-dl filterutils)
set_target_properties(replay-dl PROPERTIES COMPILE_FLAGS "-DFOR_DL")

add_executable(reconstructlog-dl reconstructlog.c)
target_link_libraries(reconstructlog-dl filterutils)

install(TARGETS dup1 dup2 purge merge replay RUNTIME DESTINATION bin/filter)

include(CMakeLists-nodist.txt OPTIONAL)
