
# All executables are somewhat optional ; testbench is the one which is
# used for tests.

# Would compile, but implies a dependency on gmp-ecm, or so it seems.
# add_executable(find_primes_hit find_primes_hit.c)
# target_link_libraries(find_primes_hit utils m)

# DNC
# add_executable(addchain addchain.c prac_bc.c)
# target_link_libraries(addchain utils m)
add_executable(testbench testbench.c)
target_link_libraries(testbench facul utils m)
# DNC
# add_executable(ecm_check ecm_check.c)
# target_link_libraries(ecm_check utils m facul)

macro(variant_facul base typewidth)
string(TOUPPER ${typewidth} TW)
    SET(src ${CMAKE_CURRENT_SOURCE_DIR}/${base}.c)
    SET(dst ${CMAKE_CURRENT_BINARY_DIR}/${base}_${typewidth}.c)
add_custom_command(OUTPUT  ${dst}
    COMMAND ln -s ${src} ${dst}
    DEPENDS  ${src})
set_source_files_properties(${dst} PROPERTIES COMPILE_FLAGS -DMODREDC${TW})
endmacro(variant_facul)

variant_facul(pm1 15ul)
variant_facul(pp1 15ul)
variant_facul(ecm 15ul)
variant_facul(facul_doit 15ul)
variant_facul(pm1 ul)
variant_facul(pp1 ul)
variant_facul(ecm ul)
variant_facul(facul_doit ul)

add_library(facul
        ${CMAKE_CURRENT_BINARY_DIR}/ecm_ul.c
        ${CMAKE_CURRENT_BINARY_DIR}/pm1_ul.c
        ${CMAKE_CURRENT_BINARY_DIR}/pp1_ul.c
        ${CMAKE_CURRENT_BINARY_DIR}/facul_doit_ul.c
        ${CMAKE_CURRENT_BINARY_DIR}/ecm_15ul.c
        ${CMAKE_CURRENT_BINARY_DIR}/pm1_15ul.c
        ${CMAKE_CURRENT_BINARY_DIR}/pp1_15ul.c
        ${CMAKE_CURRENT_BINARY_DIR}/facul_doit_15ul.c
        stage2.c prac_bc.c
	makeplan.c facul.c)

