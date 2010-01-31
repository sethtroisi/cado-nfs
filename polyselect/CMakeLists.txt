
add_executable (polyselect polyselect.c auxiliary.c rootsieve.c)
add_executable (rotate rotate.c auxiliary.c rootsieve.c)
add_executable (rotate_all rotate_all.c auxiliary.c rootsieve.c)

target_link_libraries (polyselect utils m)
target_link_libraries (rotate utils m)
target_link_libraries (rotate_all utils m)

install(TARGETS polyselect rotate rotate_all RUNTIME DESTINATION bin/polyselect)
include(CMakeLists-nodist.txt OPTIONAL)
