# defining a library avoids compiling the same oject file once for each
# binary. The choice of setting it as a shared library is completely
# arbitrary here. A static would do as well.

# Root OPTimization sources.
set (ropt_sources 
  ropt.c ropt_linear.c ropt_quadratic.c ropt_str.c ropt_param.c
  ropt_tree.c ropt_arith.c ropt_io.c ropt_stage1.c ropt_stage2.c
  )

add_library(polyselect_common ${CADO_SHARED} auxiliary.c rho.c murphyE.c area.c size_optimization.c)
add_executable (polyselect2l polyselect2l.c polyselect2l_str.c polyselect2l_arith.c ${ropt_sources})
add_executable (polyselect_ropt polyselect_ropt.c ${ropt_sources})
add_executable (polyselect3 polyselect3.c ${ropt_sources})

# the utils library already includes -lgmp and -lm
target_link_libraries (polyselect_common utils)
target_link_libraries (polyselect2l polyselect_common utils ${pthread_libs})
target_link_libraries (polyselect_ropt polyselect_common utils ${pthread_libs})
target_link_libraries (polyselect3 polyselect_common utils ${pthread_libs})

if(ENABLE_SHARED)
    install(TARGETS polyselect_common DESTINATION ${LIBSUFFIX})
endif(ENABLE_SHARED)
install(TARGETS polyselect2l RUNTIME DESTINATION ${BINSUFFIX}/polyselect)
install(TARGETS polyselect_ropt RUNTIME DESTINATION ${BINSUFFIX}/polyselect)
include(CMakeLists-nodist.txt OPTIONAL)
