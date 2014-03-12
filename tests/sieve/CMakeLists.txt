# CMake 2.8.4 and newer has a WORKING_DIRECTORY property, but we may want
# to stay compatible with older versions

# Build makefb and las before running check
add_dependencies(check makefb las las_tracek)

# Get paths to the makefb and las binaries and store in variables
# makefb_path and las_path, resp.
get_property(makefb_path TARGET makefb PROPERTY LOCATION)
get_property(las_path TARGET las PROPERTY LOCATION)
get_property(las_tracek_path TARGET las_tracek PROPERTY LOCATION)

add_test(F9_sievetest bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.sievetest.sh\" \"${makefb_path}\" \"${las_path}\" \"${CADO_NFS_SOURCE_DIR}\"")
add_test(F9_tracektest bash -c "\"${CMAKE_CURRENT_SOURCE_DIR}/F9.tracektest.sh\" \"${makefb_path}\" \"${las_tracek_path}\" \"${CADO_NFS_SOURCE_DIR}\"")

# Try with default TRIALDIV_MAXLEN, whatever it may be
add_executable(test_trialdiv test_trialdiv.c "${CADO_NFS_SOURCE_DIR}/sieve/trialdiv.c")
target_link_libraries(test_trialdiv utils tests)
run_test(test_trialdiv test_trialdiv)

macro (add_trialdiv_test MAXLEN)
  add_executable(test_trialdiv${MAXLEN} test_trialdiv.c)
  set_property(TARGET test_trialdiv${MAXLEN} PROPERTY COMPILE_DEFINITIONS TRIALDIV_MAXLEN=${MAXLEN})
  target_link_libraries(test_trialdiv${MAXLEN} libtrialdiv utils tests)
  run_test(test_trialdiv${MAXLEN} test_trialdiv${MAXLEN})
endmacro(add_trialdiv_test)

macro (add_init_alg_norms_bucket_region)
  add_executable(test_init_alg_norms_bucket_region test_init_alg_norms_bucket_region.c)
  set_property(TARGET test_init_alg_norms_bucket_region PROPERTY COMPILE_DEFINITIONS)
  target_link_libraries(test_init_alg_norms_bucket_region las-norms utils tests)
  run_test(test_init_alg_norms_bucket_region test_init_alg_norms_bucket_region)
endmacro(add_init_alg_norms_bucket_region)

macro (add_init_rat_norms_bucket_region)
  add_executable(test_init_rat_norms_bucket_region test_init_rat_norms_bucket_region.c)
  set_property(TARGET test_init_rat_norms_bucket_region PROPERTY COMPILE_DEFINITIONS)
  target_link_libraries(test_init_rat_norms_bucket_region las-norms utils tests)
  run_test(test_init_rat_norms_bucket_region test_init_rat_norms_bucket_region)
endmacro(add_init_rat_norms_bucket_region)

# Try specific TRIALDIV_MAXLEN from 1 to 6
add_trialdiv_test(1)
add_trialdiv_test(2)
add_trialdiv_test(3)
add_trialdiv_test(4)
add_trialdiv_test(5)
add_trialdiv_test(6)

add_init_alg_norms_bucket_region()
add_init_rat_norms_bucket_region()

add_subdirectory (ecm)
