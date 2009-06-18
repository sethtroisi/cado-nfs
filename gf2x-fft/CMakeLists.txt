
# include_directories(${CADO_NFS_SOURCE_DIR}/gf2x/gf2x)

add_subdirectory(cantor)
link_directories(cantor)

add_library(gf2x-fft fake_fft.c)

# It's a kludge. But putting a bare ``gf2x'' here won't work because cmake
# can't see whether libgf2x.a has changed (it has no idea of its location !)
message("gf2x lib is ${GF2X_LIB_FILE}")
target_link_libraries(gf2x-fft ${GF2X_LIB_FILE} cantor128)
add_dependencies(gf2x-fft gf2x-build)

include(CMakeLists-nodist.txt OPTIONAL)
