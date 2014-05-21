add_executable (convert_poly convert_poly.c)
target_link_libraries (convert_poly utils)
install(TARGETS convert_poly RUNTIME DESTINATION ${BINSUFFIX}/misc)

add_executable (convert_rels convert_rels.c)
target_link_libraries (convert_rels utils)
install(TARGETS convert_rels RUNTIME DESTINATION ${BINSUFFIX}/misc)

add_executable(check_rels check_rels.c)
target_link_libraries(check_rels utils_with_io)
install(TARGETS check_rels RUNTIME DESTINATION ${BINSUFFIX}/misc)

add_executable(debug_renumber EXCLUDE_FROM_ALL debug_renumber.c)
target_link_libraries(debug_renumber utils)

include(CMakeLists-nodist.txt OPTIONAL)
