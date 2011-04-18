# defining a library avoids compiling the same oject file once for each
# binary. The choice of setting it as a shared library is completely
# arbitrary here. A static would do as well.
add_library(polyselect_common SHARED auxiliary.c rootsieve.c rho.c murphyE.c)
add_executable (polyselect polyselect.c)
add_executable (rotate rotate.c)
add_executable (rotate_all rotate_all.c)
add_executable (polyselect2 polyselect2.c rootsieve5.c)
add_executable (twoquadratics twoquadratics.c)
add_executable (polyopt polyopt.c)

target_link_libraries (polyselect polyselect_common utils m)
target_link_libraries (twoquadratics polyselect_common utils m gmp)
target_link_libraries (rotate polyselect_common utils m)
target_link_libraries (rotate_all polyselect_common utils m)
target_link_libraries (polyselect2 polyselect_common utils m ${pthread_libs})
target_link_libraries (polyopt polyselect_common utils m gmp)

install(TARGETS polyselect_common DESTINATION lib/${package_basename})
install(TARGETS polyselect rotate rotate_all polyselect2 RUNTIME
    DESTINATION bin/polyselect)
include(CMakeLists-nodist.txt OPTIONAL)
