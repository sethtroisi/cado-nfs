# CMake 2.8.4 and newer has a WORKING_DIRECTORY property, but we may want
# to stay compatible with older versions
add_test(F9_sievetest "${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh" "${CADO_NFS_SOURCE_DIR}" "${CADO_NFS_BINARY_DIR}")

# Try with default TRIALDIV_MAXLEN, whatever it may be
add_executable(test_trialdiv test_trialdiv.c "${CADO_NFS_SOURCE_DIR}/sieve/trialdiv.c")
target_link_libraries(test_trialdiv utils tests)
run_test(test_trialdiv test_trialdiv)

macro (add_trialdiv_test MAXLEN)
# add_executable(test_trialdiv${MAXLEN} test_trialdiv.c "${CADO_NFS_SOURCE_DIR}/sieve/trialdiv.c")
add_executable(test_trialdiv${MAXLEN} test_trialdiv.c)
set_property(TARGET test_trialdiv${MAXLEN} PROPERTY COMPILE_DEFINITIONS TRIALDIV_MAXLEN=${MAXLEN})
target_link_libraries(test_trialdiv${MAXLEN} libtrialdiv utils tests)
run_test(test_trialdiv${MAXLEN} test_trialdiv${MAXLEN})
endmacro(add_trialdiv_test)

# Try specific TRIALDIV_MAXLEN from 1 to 6
add_trialdiv_test(1)
add_trialdiv_test(2)
add_trialdiv_test(3)
add_trialdiv_test(4)
add_trialdiv_test(5)
add_trialdiv_test(6)

add_subdirectory (ecm)
