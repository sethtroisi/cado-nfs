
add_executable (polyselect polyselect.c auxiliary.c)
add_executable (rotate rotate.c auxiliary.c)

target_link_libraries (polyselect utils m)
target_link_libraries (rotate utils m)

install(TARGETS polyselect rotate RUNTIME DESTINATION bin/polyselect)
