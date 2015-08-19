add_custom_target(all_test_dependencies)

# cmake understands tests and targets as being two different things. It
# is perfectly common to have both a target and a test having the same
# name X. Yet testing X will not trigger building X, and that is
# admittedly unfortunate.

# INTERNALS:
#
# we define several macros for dealing with this in a better way.
#
#  For each test X, we define another test builddep_X.
#
#  The test builddep_X which only takes the necessary steps for building
#  the dependencies of the test X. Those dependencies form a target
#  named X_dependencies.
#
#  The target X_dependencies may perhaps include a target named X.
#
#  The test X is made to depend on builddep_X
#
#  The target all_test_dependencies is made to depend on X_dependencies
#
#  Tests may have other dependency relationships, which should be set
#  in a manner which is compatible with our dependency X <- builddep_X.
#  This is taken care of by update_test_build_depends

# WHICH CMAKE COMMANDS SHOULD BE USED ?
#
# only one of the following:
#
#       build_and_run_test
#       run_target_as_a_test
#       run_target_as_a_parametric_test
#       run_test
#
# the naming is admittedly wholly inconsistent, and would deserve some
# cleaning.
# 
# 1 Typical case: you have one binary foo, which makes one test
# 
#  1a  foo is compiled only from foo.c, and perhaps links with library
#      baz and qux. It may need arguments 17 and 42
#
#      build_and_run_test(foo.c "baz;qux" 17 42)
#
#      note that if no libraries dependencies exist, those should be
#      specified as "".
#
#      [NOTE: this is equivalent to:
#          add_executable(foo foo.c)
#          target_link_libraries(foo baz qux)
#          run_test(foo foo 17 42)               # test foo versus target foo
#          set_test_build_depends(foo foo)
#      ]
#
#  1b  foo is compiled from foo.c and bar.c, and/or the binary foo needs
#      to be invoked with arguments 17 and 42
#      
#          add_executable(foo foo.c bar.c)
#          target_link_libraries(foo baz qux)
#          run_target_as_a_test(foo 17 42)
# 
# 2 you have one binary, which makes several tests depending on
#   arguments.
#
#          add_executable(foo foo.c bar.c)
#          target_link_libraries(foo baz qux)
#          run_target_as_a_parametric_test(1 foo 17)
#          run_target_as_a_parametric_test(2 foo 42)
# 
#   this creates tests foo_1 and foo_2
#
# 3 you have one script, which needs several binaries. In this case, you
#   need to use run_test and set the dependencies by yourself.
#
#          run_test(testname ${CMAKE_CURRENT_SOURCE_DIR}/test_main_script.sh ${CADO_NFS_SOURCE_DIR} ${CADO_NFS_BINARY_DIR})
#          set_test_build_depends(testname binary1 binary2)
#
#         (arguments passed to the script may vary of course, but needing
#         the bindir and/or srcdir can be considered typical)

macro(set_test_build_depends test_name)
    add_custom_target(${test_name}_dependencies)
    add_dependencies(${test_name}_dependencies ${ARGN})
    add_dependencies(all_test_dependencies ${test_name}_dependencies)
    add_test("builddep_${test_name}" "${CMAKE_COMMAND}" --build
        ${CMAKE_BINARY_DIR} --target ${test_name}_dependencies)
    set("${test_name}_has_build_depends" 1)
    # message(STATUS "Creating dependencies for ${test_name} -- ${ARGN}")
    set_tests_properties (${test_name} PROPERTIES DEPENDS "builddep_${test_name}")
endmacro(set_test_build_depends test_name ${ARGN})

macro(update_test_build_depends test_name)
    add_dependencies(${test_name}_dependencies ${ARGN})
endmacro(update_test_build_depends test_name)

if(CMAKE_VERSION VERSION_GREATER 3.1)
    cmake_policy(SET CMP0054 OLD)
endif(CMAKE_VERSION VERSION_GREATER 3.1)

macro(set_test_dependencies test_name)
    set(deps ${ARGN})
    if("${test_name}_has_build_depends")
        list(APPEND deps "builddep_${test_name}")
    endif("${test_name}_has_build_depends")
    # cmake has very quirky syntax...
    string_join(alldeps ";" ${deps})
    # message(STATUS "Updating dependencies for ${test_name} -- ${alldeps}")
    set_tests_properties (${test_name} PROPERTIES DEPENDS "${alldeps}")
endmacro(set_test_dependencies test_name ${ARGN})

macro(run_test NAME DSTFILE)
    # this is like add_test, but copes with scripts better.
    set(binary ${DSTFILE})
    if (${DSTFILE} MATCHES ".sh$")
        set(binary env bash ${DSTFILE})
    endif (${DSTFILE} MATCHES ".sh$")
    add_test (${NAME} ${binary} ${ARGN})
# old version.
#  # This add_test() call causes the ctest program to be called on the
#  # current source and binary directories. The ctest program then builds
#  # the DSTFILE binary by calling cmake (which we suppress here) and make,
#  # and then runs DSTFILE
#     add_test ("${NAME}" "${CMAKE_CTEST_COMMAND}"
#               # Specify the source and build directories for the Cmake/make
#               # steps that gets executed by CTest (except we suppress the
#               # CMake step in this CTest run, but we still have to specify the
#               # directories)
#               --build-and-test "${CMAKE_CURRENT_SOURCE_DIR}"
#                                "${CMAKE_CURRENT_BINARY_DIR}"
# 
#               # If CTest were to run Cmake, it should use same CMake generator
#               # as we use in the parent Cmake. We suppress the CMake call, but
#               # the "--build-generator" parameter is mandatory so we have to
#               # specify it anyway. Without it, we get this error:
#               # "--build-and-test requires that both the generator and 
#               #  makeprogram be provided using the --build-generator and 
#               #  --build-makeprogram command line options."
#               --build-generator "${CMAKE_GENERATOR}"
# 
#               # Use same make program as the parent CMake.
#               --build-makeprogram "${CMAKE_MAKE_PROGRAM}"
# 
#               # Do NOT run CMake. The CMakeCache.txt file generated by having
#               # add_subdirectory() in the parent directories is fine, and
#               # if we were to let CTest re-run CMake, we'd have to figure out
#               # directory locations all over again - this way, we inherit
#               # everything from the CMakeCache.txt generated by the parent
#               # CMake.
#               --build-nocmake
# 
#               # Normally, CTest runs "make clean" before building any target,
#               # so that only the most recently built target would remain in
#               # the build directory. I assume we don't want that, but want to
#               # be able to run them manually as well. Do tell CTest not to do
#               # "make clean."
#               --build-noclean
# 
#               # When we build and run this test, build only the binary
#               # required for this test
#               --build-target "${DSTFILE}"
# 
#               # If any test fails, print its output to stdout.
#               # ... except that does not actually work. When running the
#               # test suite with
#               # /usr/bin/ctest --force-new-ctest-process -V
#               # then the debugging output shows --output-on-failure being
#               # passed to the spawned ctest processes, but without -V, they
#               # still don't print anything in case of failure. Maybe a bug.
#               # To get output, we have to pass it via ARGS to "make test".
#               --output-on-failure
# 
#               # And finally the binary that should be run as the actual test
#               --test-command ${binary} ${ARGN})
endmacro(run_test NAME DSTFILE)

# The SRCFILE parameter to this macro is the file name of the source file to
# build. The LIBRARIES parameter is any libraries that should be linked in.
# Any additional parameters are added to the command line of the test
# executable that gets run.
macro(build_and_run_test SRCFILE LIBRARIES)
    # The name of the binary executable is equal to the name of the source
    # file, minus the extension
    get_filename_component (DSTFILE ${SRCFILE} NAME_WE)

    add_executable (${DSTFILE} ${SRCFILE})
    target_link_libraries (${DSTFILE} ${LIBRARIES})
    # The name of the test here is the same as the name of the test
    # executable file
    add_test (${DSTFILE} ${DSTFILE} ${ARGN})
    set_test_build_depends(${DSTFILE} ${DSTFILE})
endmacro(build_and_run_test)

macro(run_target_as_a_test x)
    add_test(${x} ${x} ${ARGN})
    set_test_build_depends(${x} ${x})
endmacro(run_target_as_a_test x)

macro(run_target_as_a_parametric_test n x)
    add_test(${x}_${n} ${x} ${ARGN})
    set_test_build_depends(${x}_${n} ${x})
endmacro(run_target_as_a_parametric_test n x)

if(DEFINED ENV{CHECKS_EXPENSIVE})
    run_test(full_c59 ${CADO_NFS_SOURCE_DIR}/factor.sh 90377629292003121684002147101760858109247336549001090677693 -t 2)
    set_test_build_depends(full_c59 polyselect2l polyselect_ropt freerel makefb las dup1 dup2 purge merge merge-dl replay replay-dl sm reconstructlog-dl characters bwc_full plingen_pz sqrt)
    find_program(MAGMA magma)
    if (MAGMA)
        run_test(full_p59 ${CMAKE_CURRENT_SOURCE_DIR}/test_full_p59.sh ${CADO_NFS_SOURCE_DIR})
        set_test_build_depends(full_p59 full_c59_dependencies debug_renumber)
        run_test(test_gfp3 ${CMAKE_CURRENT_SOURCE_DIR}/test_gfp3.sh ${CADO_NFS_SOURCE_DIR} ${CADO_NFS_BINARY_DIR})
        set_test_build_depends(test_gfp3 full_c59_dependencies debug_renumber)
        run_test(test_gfp2 ${CMAKE_CURRENT_SOURCE_DIR}/test_gfp2.sh ${CADO_NFS_SOURCE_DIR} ${CADO_NFS_BINARY_DIR})
        set_test_build_depends(test_gfp2 full_c59_dependencies filter_galois polyselect_gfpn)
    endif(MAGMA)
endif(DEFINED ENV{CHECKS_EXPENSIVE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library (tests STATIC test_iter.c tests_common.c)
target_link_libraries (tests ${gmp_libname})

add_executable(test_tests_common test_tests_common.c)
target_link_libraries(test_tests_common tests)
run_target_as_a_parametric_test(1 test_tests_common)
run_target_as_a_parametric_test(2 test_tests_common -seed)
run_target_as_a_parametric_test(3 test_tests_common -seed 1)
set_tests_properties(test_tests_common_3 PROPERTIES PASS_REGULAR_EXPRESSION "Using random seed=1")
run_target_as_a_parametric_test(4 test_tests_common -seed a)
run_target_as_a_parametric_test(5 test_tests_common -seed "")
set_tests_properties(test_tests_common_2 test_tests_common_4 test_tests_common_5 PROPERTIES WILL_FAIL 1)
run_target_as_a_parametric_test(6 test_tests_common -iter 10)
set_tests_properties(test_tests_common_6 PROPERTIES PASS_REGULAR_EXPRESSION "Using 10 iterations")

build_and_run_test(test_iceildiv.c "")

add_subdirectory (utils)
add_subdirectory (sieve)
add_subdirectory (polyselect)
add_subdirectory (filter)
add_subdirectory (linalg)
add_subdirectory (gf2x)
