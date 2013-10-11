add_executable (convert_poly convert_poly.c)
target_link_libraries (convert_poly utils)
install(TARGETS convert_poly RUNTIME DESTINATION bin/misc)

add_executable (convert_rels convert_rels.c)
target_link_libraries (convert_rels utils)
install(TARGETS convert_rels RUNTIME DESTINATION bin/misc)

add_executable(complete_rels complete_rels.c)
target_link_libraries(complete_rels filterutils)
install(TARGETS complete_rels RUNTIME DESTINATION bin/misc)

add_executable(check_rels check_rels.c)
target_link_libraries(check_rels filterutils)
install(TARGETS check_rels RUNTIME DESTINATION bin/misc)
