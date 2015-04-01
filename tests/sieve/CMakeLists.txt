# CMake 2.8.4 and newer has a WORKING_DIRECTORY property, but we may want
# to stay compatible with older versions

# Build makefb and las before running check
include_directories(${PROJECT_SOURCE_DIR}/sieve)

# Get paths to the makefb and las binaries and store in variables
# makefb_path and las_path, resp.
# for cmake pre-3.0, we did:
# get_property(makefb_path TARGET makefb PROPERTY LOCATION)

# now it seems to be forbidden with cmake 3.0+ (see cmake --help-policy
# CMP0026). The following is probably safe enough.
set(makefb_path         ${CADO_NFS_BINARY_DIR}/sieve/makefb)
set(las_path            ${CADO_NFS_BINARY_DIR}/sieve/las)
set(las_tracek_path     ${CADO_NFS_BINARY_DIR}/sieve/las_tracek)
# FAILING # set(las_kbucket_path    ${CADO_NFS_BINARY_DIR}/sieve/las_kbucket)
# FAILING # set(las_mbucket_path    ${CADO_NFS_BINARY_DIR}/sieve/las_mbucket)
set(las_descent_path    ${CADO_NFS_BINARY_DIR}/sieve/las_descent)

set(F9_roots "${CMAKE_CURRENT_BINARY_DIR}/F9.roots")
add_test(F9_makefbtest "${makefb_path}" -poly "${CADO_NFS_SOURCE_DIR}/params/F9.poly" -alim 1200000 -maxbits 10 -out "${F9_roots}")
set_test_build_depends(F9_makefbtest makefb)

add_test(F9_sievetest bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh\" \"${F9_roots}\" \"${las_path}\" \"${CADO_NFS_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}/F9.checksums\"")
# Before running F9_sievetest, remove the file with the captured checksums
# We have to make this a "test" so that ctest runs it, even though it does
# not actually test anything
add_test(F9_sievetest_cleanup "${CMAKE_COMMAND}" -E remove -f "${CMAKE_CURRENT_BINARY_DIR}/F9.checksums")
set_test_build_depends(F9_sievetest las)
set_test_dependencies(F9_sievetest F9_sievetest_cleanup F9_makefbtest)

# We make these depend on sievetest, for checking the checksum file.
# Since sievetest depends on makefbtest, we don't have to specify it
# again.
add_test(F9_sievetest_unsieve1 bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh\" \"${F9_roots}\" \"${las_path}\" \"${CADO_NFS_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}/F9.checksums\" -unsievethresh 1")
set_test_build_depends(F9_sievetest_unsieve1 las)
set_test_dependencies(F9_sievetest_unsieve1 F9_sievetest)

# FAILING # add_test(F9_kbucket_test bash -c "env EXPECTED_FAIL=1 \"${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh\" \"${F9_roots}\" \"${las_kbucket_path}\" \"${CADO_NFS_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}/F9.checksums\"")
# FAILING # SET_TESTS_PROPERTIES(F9_kbucket_test PROPERTIES DEPENDS F9_sievetest WILL_FAIL True)
# FAILING # add_test(F9_mbucket_test bash -c "env EXPECTED_FAIL=1 \"${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh\" \"${F9_roots}\" \"${las_mbucket_path}\" \"${CADO_NFS_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}/F9.checksums\"")
# FAILING # SET_TESTS_PROPERTIES(F9_mbucket_test PROPERTIES DEPENDS F9_sievetest WILL_FAIL True)
add_test(F9_cofactest bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.cofactest.sh\" \"${F9_roots}\" \"${las_path}\" \"${CADO_NFS_SOURCE_DIR}\" \"${CMAKE_CURRENT_BINARY_DIR}/F9.checksums\" -file-cofact \"${CADO_NFS_SOURCE_DIR}/params/F9.st\"")
set_test_build_depends(F9_cofactest las)
set_test_dependencies(F9_cofactest F9_sievetest)

# The TRACE_K test does only a single special-q, so we don't compare output
# relations or checksums
add_test(F9_tracektest bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.tracektest.sh\" \"${F9_roots}\" \"${las_tracek_path}\" \"${CADO_NFS_SOURCE_DIR}\" -traceab -8517,584707")
set_test_build_depends(F9_tracektest las_tracek)
set_test_dependencies(F9_tracektest F9_makefbtest)

# Try a too-large lambda value
add_test(F9_rlambdatest bash -c "rlambda=5 \"${CMAKE_CURRENT_SOURCE_DIR}/F9.tracektest.sh\" \"${F9_roots}\" \"${las_path}\" \"${CADO_NFS_SOURCE_DIR}\" -regex \"rlambda>.* does not make sense\"")
set_test_build_depends(F9_rlambdatest las)
set_test_dependencies(F9_rlambdatest F9_makefbtest)

add_test(F9_alambdatest bash -c "alambda=5 \"${CMAKE_CURRENT_SOURCE_DIR}/F9.tracektest.sh\" \"${F9_roots}\" \"${las_path}\" \"${CADO_NFS_SOURCE_DIR}\" -regex \"alambda>.* does not make sense\"")
set_test_build_depends(F9_alambdatest las)
set_test_dependencies(F9_alambdatest F9_makefbtest)

# Try the todo-list mode of las
add_test(F9_todotest bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.todotest.sh\" \"${F9_roots}\" \"${las_path}\" \"${CADO_NFS_SOURCE_DIR}\" ${CADO_NFS_SOURCE_DIR}/tests/sieve/F9.todo")
set_test_build_depends(F9_todotest las)
set_test_dependencies(F9_todotest F9_makefbtest)

# Try the descent 
add_test(F9_descenttest bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.descenttest.sh\" \"${F9_roots}\" \"${las_descent_path}\" \"${CADO_NFS_SOURCE_DIR}\" ${CADO_NFS_SOURCE_DIR}/tests/sieve/F9.hint")
set_test_build_depends(F9_descenttest las_descent)
set_test_dependencies(F9_descenttest F9_makefbtest)

# Try with default TRIALDIV_MAXLEN, whatever it may be
add_executable(test_trialdiv test_trialdiv.c "${CADO_NFS_SOURCE_DIR}/sieve/trialdiv.c")
target_link_libraries(test_trialdiv utils tests)
run_test(test_trialdiv test_trialdiv)

build_and_run_test(test_memset_things "las-norms;utils")

macro (add_trialdiv_test MAXLEN)
  add_executable(test_trialdiv${MAXLEN} test_trialdiv.c)
  set_property(TARGET test_trialdiv${MAXLEN} PROPERTY COMPILE_DEFINITIONS TRIALDIV_MAXLEN=${MAXLEN})
  target_link_libraries(test_trialdiv${MAXLEN} trialdiv utils tests)
  run_test(test_trialdiv${MAXLEN} test_trialdiv${MAXLEN})
endmacro(add_trialdiv_test)

macro (add_init_norms_bucket_region)
  add_executable(test_init_norms_bucket_region test_init_norms_bucket_region.cpp )
  set_property(TARGET test_init_norms_bucket_region PROPERTY COMPILE_DEFINITIONS)
  target_link_libraries(test_init_norms_bucket_region las-norms utils tests)
  run_test(test_init_norms_bucket_region test_init_norms_bucket_region)
endmacro(add_init_norms_bucket_region)

# Try specific TRIALDIV_MAXLEN from 1 to 6
add_trialdiv_test(1)
add_trialdiv_test(2)
add_trialdiv_test(3)
add_trialdiv_test(4)
add_trialdiv_test(5)
add_trialdiv_test(6)

add_init_norms_bucket_region()

add_executable(test_fb_make_step test_fb_make_step.c "${CADO_NFS_SOURCE_DIR}/sieve/fb.cpp")
target_link_libraries(test_fb_make_step utils tests)

add_subdirectory (ecm)
add_subdirectory (strategies)
