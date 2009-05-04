
add_executable (polyselect polyselect.c aux.c)
add_executable (rotate rotate.c aux.c)

target_link_libraries (polyselect utils m)
target_link_libraries (rotate utils m)

install(TARGETS polyselect rotate RUNTIME DESTINATION bin/polyselect)
