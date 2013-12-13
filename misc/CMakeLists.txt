link_directories (${CADO_NFS_BINARY_DIR}/ffs/f2)
link_directories (${CADO_NFS_BINARY_DIR}/ffs/f3)
include_directories (${CADO_NFS_SOURCE_DIR}/ffs)

add_executable (convert_poly convert_poly.c)
target_link_libraries (convert_poly utils)
install(TARGETS convert_poly RUNTIME DESTINATION bin/misc)

add_executable (convert_rels convert_rels.c)
target_link_libraries (convert_rels utils)
install(TARGETS convert_rels RUNTIME DESTINATION bin/misc)

add_executable(check_rels check_rels.c)
target_link_libraries(check_rels filterutils)
install(TARGETS check_rels RUNTIME DESTINATION bin/misc)

add_executable(debug_renumber EXCLUDE_FROM_ALL debug_renumber.c)
target_link_libraries(debug_renumber utils)

add_executable(debug_renumber-ffs-f2 EXCLUDE_FROM_ALL debug_renumber.c)
target_link_libraries(debug_renumber-ffs-f2 ffs-f2)
set_target_properties(debug_renumber-ffs-f2 PROPERTIES COMPILE_FLAGS "-DFOR_FFS -DUSE_F2")
add_executable(debug_renumber-ffs-f3 EXCLUDE_FROM_ALL debug_renumber.c)
target_link_libraries(debug_renumber-ffs-f3 ffs-f3)
set_target_properties(debug_renumber-ffs-f3 PROPERTIES COMPILE_FLAGS "-DFOR_FFS -DUSE_F3")
