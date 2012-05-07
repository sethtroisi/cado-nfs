add_executable (convert_poly convert_poly.c)

target_link_libraries (convert_poly utils)

install(TARGETS convert_poly RUNTIME DESTINATION bin/misc)
