
add_executable(characters characters.c gauss.c blockmatrix.c bit_matrices.c)
target_link_libraries(characters utils)

# This library is used by bwc mostly, but since its files are relevant here as
# well...
add_library(mat_io readbuffer.c readmat.c readmat-easy.c)

add_executable(transpose EXCLUDE_FROM_ALL transpose.c readmat.c)
target_link_libraries(transpose utils m)

set_source_files_properties (gauss.c COMPILE_FLAGS "-DMULTI_ROW=3 -DNO_MAIN")

include(CMakeLists-nodist.txt OPTIONAL)
