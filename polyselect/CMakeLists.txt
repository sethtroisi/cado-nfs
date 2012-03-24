# defining a library avoids compiling the same oject file once for each
# binary. The choice of setting it as a shared library is completely
# arbitrary here. A static would do as well.
add_library(polyselect_common SHARED auxiliary.c rootsieve.c rho.c murphyE.c)
add_executable (rotate rotate.c)
add_executable (rotate_all rotate_all.c)
add_executable (polyselect2 polyselect2.c ropt.c
                ropt_param.c ropt_tree.c ropt_arith.c
                ropt_io.c ropt_stage1.c ropt_stage2.c)
add_executable (polyselect2l polyselect2l.c ropt.c
                ropt_param.c ropt_tree.c ropt_arith.c
                ropt_io.c ropt_stage1.c ropt_stage2.c)
add_executable (twoquadratics twoquadratics.c)
add_executable (sopt sopt.c)

# the utils library already includes -lgmp and -lm
target_link_libraries (polyselect_common utils)
target_link_libraries (twoquadratics polyselect_common utils)
target_link_libraries (rotate polyselect_common utils)
target_link_libraries (rotate_all polyselect_common utils)
target_link_libraries (polyselect2 polyselect_common utils ${pthread_libs})
target_link_libraries (polyselect2l polyselect_common utils ${pthread_libs})
target_link_libraries (sopt polyselect_common utils)

install(TARGETS polyselect_common DESTINATION lib/${package_basename})
install(TARGETS rotate rotate_all polyselect2 RUNTIME DESTINATION bin/polyselect)
include(CMakeLists-nodist.txt OPTIONAL)
