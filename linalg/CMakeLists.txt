
include_directories(bwc) # for worker-threads.[ch]

add_executable(characters characters.c gauss.c blockmatrix.c
    bit_matrices.c bwc/worker-threads.c)
target_link_libraries(characters utils pthread)
install(TARGETS characters RUNTIME DESTINATION bin/linalg)

# This library is used by bwc mostly, but since its files are relevant here as
# well... Note that because this library is required by the low-level
# matmul libraries in bwc, which are dlopen'ed, it is necessary to tag it
# as shared.
add_library(mat_io SHARED readbuffer.c)
install(TARGETS mat_io LIBRARY DESTINATION "lib/${package_basename}")

# add_executable(transpose EXCLUDE_FROM_ALL transpose.c readmat.c)
# target_link_libraries(transpose utils m)

set_source_files_properties (gauss.c COMPILE_FLAGS "-DMULTI_ROW=3 -DNO_MAIN")

include(CMakeLists-nodist.txt OPTIONAL)
