
include_directories(bwc) # for worker-threads.[ch]

add_executable(characters characters.c gauss.c blockmatrix.c
    bit_matrices.c bwc/worker-threads.c)

target_link_libraries(characters utils_with_io pthread)

install(TARGETS characters RUNTIME DESTINATION ${BINSUFFIX}/linalg)

include(CMakeLists-nodist.txt OPTIONAL)
