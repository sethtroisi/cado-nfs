add_executable (polyselect polyselect.c auxiliary.c rootsieve.c rho.c murphyE.c)
add_executable (rotate rotate.c auxiliary.c rootsieve.c murphyE.c rho.c)
add_executable (rotate_all rotate_all.c auxiliary.c rootsieve.c rho.c murphyE.c)
add_executable (kleinjung2a kleinjung2a.c auxiliary.c murphyE.c rho.c)

target_link_libraries (polyselect utils m)
target_link_libraries (rotate utils m)
target_link_libraries (rotate_all utils m)
target_link_libraries (kleinjung2a utils m ${pthread_libs})

install(TARGETS polyselect rotate rotate_all RUNTIME DESTINATION bin/polyselect)
include(CMakeLists-nodist.txt OPTIONAL)
