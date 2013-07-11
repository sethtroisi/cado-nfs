
include_directories(bwc) # for worker-threads.[ch]

add_executable(characters characters.c gauss.c blockmatrix.c
    bit_matrices.c bwc/worker-threads.c)
add_executable(sm sm.c)

target_link_libraries(characters utils pthread)
target_link_libraries(sm utils pthread)

install(TARGETS characters RUNTIME DESTINATION bin/linalg)
install(TARGETS sm RUNTIME DESTINATION bin/linalg)

set_source_files_properties (gauss.c COMPILE_FLAGS "-DMULTI_ROW=3 -DNO_MAIN")

include(CMakeLists-nodist.txt OPTIONAL)
