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
# only the following:
#
#       cado_define_test()
#       cado_divert_test()
#
# those are documented below within their source definitions.
#
# Note that for test X, there are some use cases where it is nice to know
# about X_dependencies (for adding extra dependencies late), or the
# targets X (for setting compile time definitions).
if(CMAKE_VERSION VERSION_GREATER 3.1)
    cmake_policy(SET CMP0054 OLD)
endif(CMAKE_VERSION VERSION_GREATER 3.1)

macro(cado_define_test)
    # This macro is appropriate for a test which consists in running a
    # binary. The different ways to call this macro are:
    # 
    # cado_define_test(SOURCES foo.c bar.c
    #                        [TEST_NAME x]
    #                        [LIBRARIES x y z]
    #                        [TARGET_DEPENDENCIES x y z]
    #                        [TEST_DEPENDENCIES x y z]
    #                        [ARGUMENTS x y z]
    #                        [NO_DEFAULT_RUN])
    # cado_define_test(foo.c bar.c
    #                        [TEST_NAME x]
    #                        [LIBRARIES x y z]
    #                        [TARGET_DEPENDENCIES x y z]
    #                        [TEST_DEPENDENCIES x y z]
    #                        [ARGUMENTS x y z]
    #                        [NO_DEFAULT_RUN])
    # cado_define_test(foo
    #                        SOURCES foo.c bar.c
    #                        [LIBRARIES x y z]
    #                        [TARGET_DEPENDENCIES x y z]
    #                        [TEST_DEPENDENCIES x y z]
    #                        [ARGUMENTS x y z]
    #                        [NO_DEFAULT_RUN])
    # 
    # cado_define_test(foo
    #                        SCRIPT test.sh
    #                        [TARGET_DEPENDENCIES x y z]
    #                        [TEST_DEPENDENCIES x y z]
    #                        [ARGUMENTS x y z]
    #                        [NO_DEFAULT_RUN])
    # 
    # we define the "implicit" parameters as those occuring before any
    # tag of the form [SOURCES|LIBRARIES|...|...] (e.g foo.c bar.c in the
    # second example above).
    #
    # The test is named in either of the following ways:
    #   - after the TEST_NAME parameter if one is given,
    #   - as the implicit parameter if there is an explicit SOURCES or
    #     SCRIPT parameter and implicit parameters (there must not be
    #     more than one, then).
    #   - after the basename (with extension stripped) of the first source
    #     file or the SCRIPT parameter otherwise.
    #
    # The source files to be compiled are those defined by the SOURCES
    # parameter, or the implicit parameters otherwise.
    #
    # If the test consists in running a script (or any arbitrary external
    # command, really), then an explicit SCRIPT parameter must be
    # present, and explicit or implicit SOURCES parameters are forbidden.
    # If the SCRIPT parameter consists of more than one parameter and no
    # explicit ARGUMENTS parameter list is given, then the tail of the
    # SCRIPT parameter list is taken as an argument list.
    #
    # the PROGRAM tag has the exact same meaning as SCRIPT (handling is
    # identical)
    #
    # LIBRARIES, TARGET_DEPENDENCIES, TEST_DEPENDENCIES should be self
    # explanatory. LIBRARIES is incompatible with SCRIPT.
    #
    # ARGUMENTS defines which arguments are to be passed to the binary.
    # Optional.
    #
    # NO_DEFAULT_RUN (incompatible with ARGUMENTS) indicates that no test
    # is really defined, only a baseline for further calls to
    # cado_divert_test. (note that despite the fact that no test is
    # created in cmake's mind, we do need a test name for further
    # reference in cado_divert_test).
    set(IMPLICIT)
    set(TEST_DEPENDENCIES)
    set(TARGET_DEPENDENCIES)
    set(LIBRARIES)
    set(SOURCES)
    set(SCRIPT)
    set(ARGUMENTS)
    set(TEST_NAME)
    set(NO_DEFAULT_RUN)
    set(current IMPLICIT)
    foreach(x ${ARGN})
        if (17 EQUAL 42)        # cmake syntax sucks
            # dead branch, I know. But since the 1st test has privileged
            # status in cmake syntax, let's make it a straw man.
        elseif (x STREQUAL "TEST_DEPENDENCIES")
            SET(current TEST_DEPENDENCIES)
        elseif (x STREQUAL "TARGET_DEPENDENCIES")
            SET(current TARGET_DEPENDENCIES)
        elseif (x STREQUAL "LIBRARIES")
            SET(current LIBRARIES)
        elseif (x STREQUAL "SOURCES")
            SET(current SOURCES)
        elseif ((x STREQUAL "SCRIPT") OR (x STREQUAL "PROGRAM"))
            SET(current SCRIPT)
        elseif (x STREQUAL "ARGUMENTS")
            SET(current ARGUMENTS)
        elseif (x STREQUAL "TEST_NAME")
            SET(current TEST_NAME)
        elseif (x STREQUAL "NO_DEFAULT_RUN")
            SET(current NO_DEFAULT_RUN)
            # special, because it's a flag, really
            set(NO_DEFAULT_RUN 1)
        else(17 EQUAL 42)       # cmake syntax sucks
            list(APPEND ${current} ${x})
        endif(17 EQUAL 42)      # cmake syntax sucks
    endforeach(x)
    # do some checking.
	#message(STATUS "IMPLICIT is ${IMPLICIT}")
	#message(STATUS "TEST_DEPENDENCIES is ${TEST_DEPENDENCIES}")
	#message(STATUS "TARGET_DEPENDENCIES is ${TARGET_DEPENDENCIES}")
	#message(STATUS "LIBRARIES is ${LIBRARIES}")
	#message(STATUS "SOURCES is ${SOURCES}")
	#message(STATUS "SCRIPT is ${SCRIPT}")
	#message(STATUS "ARGUMENTS is ${ARGUMENTS}")
	#message(STATUS "TEST_NAME is ${TEST_NAME}")
	#message(STATUS "NO_DEFAULT_RUN is ${NO_DEFAULT_RUN}")
    list(LENGTH TEST_NAME nTEST_NAME)
    list(LENGTH SOURCES nSOURCES)
    list(LENGTH SCRIPT nSCRIPT)
    list(LENGTH LIBRARIES nLIBRARIES)
    list(LENGTH IMPLICIT nIMPLICIT)
    list(LENGTH NO_DEFAULT_RUN nNO_DEFAULT_RUN)
    list(LENGTH TARGET_DEPENDENCIES nTARGET_DEPENDENCIES)
    list(LENGTH ARGUMENTS nARGUMENTS)
    # direct the implicit parameters somewhere
    if(nIMPLICIT GREATER 0) 
        if (17 EQUAL 42)        # cmake syntax sucks
        elseif (((nSOURCES GREATER 0) OR (nSCRIPT GREATER 0)) AND (nTEST_NAME GREATER 0))
            message(FATAL_ERROR "bad syntax with implicit parameter list and both (SOURCES or SCRIPT) and TEST_NAME defined")
        elseif (((nSOURCES GREATER 0) OR (nSCRIPT GREATER 0)) AND (nTEST_NAME EQUAL 0))
            set(TEST_NAME ${IMPLICIT})
            list(LENGTH TEST_NAME nTEST_NAME)
        # from then on we know that SOURCES and SCRIPT are empty
        else(17 EQUAL 42)       # cmake syntax sucks
            set(SOURCES ${IMPLICIT})
            list(LENGTH SOURCES nSOURCES)
        endif (17 EQUAL 42)        # cmake syntax sucks
    endif(nIMPLICIT GREATER 0) 
    if(nSCRIPT GREATER 1)
        if (nARGUMENTS EQUAL 0)
            # Then use the tail as an argument list
            LIST(GET SCRIPT 0 x)
            LIST(REMOVE_AT SCRIPT 0)
            SET(ARGUMENTS ${SCRIPT})
            SET(SCRIPT ${x})
            list(LENGTH SCRIPT nSCRIPT)
        else (nARGUMENTS EQUAL 0)
            message(WARNING "too many script names, retaining only first")
            LIST(GET SCRIPT 0 x)
            SET(SCRIPT ${x})
            list(LENGTH SCRIPT nSCRIPT)
        endif (nARGUMENTS EQUAL 0)
    endif(nSCRIPT GREATER 1)
    if(nSCRIPT GREATER 0)
        if (nSOURCES GREATER 0)
            message(FATAL "SCRIPT and SOURCES incompatible")
        endif (nSOURCES GREATER 0)
        if (nLIBRARIES GREATER 0)
            message(FATAL "SCRIPT and LIBRARIES incompatible")
        endif (nLIBRARIES GREATER 0)
    endif(nSCRIPT GREATER 0)
    if(nTEST_NAME GREATER 1)
        message(WARNING "too many test names ${TEST_NAME}, retaining only first")
        LIST(GET TEST_NAME 0 x)
        SET(TEST_NAME ${x})
    elseif(nTEST_NAME EQUAL 0)
        # # message(STATUS "trying to find a test name")
        # # message(STATUS "SOURCES = ${SOURCES}")
        # # message(STATUS "SCRIPT = ${SCRIPT}")
        if (17 EQUAL 42)        # cmake syntax sucks
        elseif(nSOURCES GREATER 0)
            # define the test name as the base name without extension of the
            # first source file.
            LIST(GET SOURCES 0 x)
            get_filename_component (TEST_NAME ${x} NAME_WE)
        elseif(nSCRIPT GREATER 0)
            # define the test name as the base name without extension of the
            # first source file.
            LIST(GET SCRIPT 0 x)
            get_filename_component (TEST_NAME ${x} NAME_WE)
        else (17 EQUAL 42)        # cmake syntax sucks
            message(FATAL_ERROR "cannot find a name for test")
        endif (17 EQUAL 42)        # cmake syntax sucks
    endif(nTEST_NAME GREATER 1)
    if(nNO_DEFAULT_RUN GREATER 1)
        message(FATAL_ERROR "discarded arguments after NO_DEFAULT_RUN for test ${TEST_NAME}")
        set(NO_DEFAULT_RUN 1)
    elseif(nNO_DEFAULT_RUN GREATER 0)
        list(LENGTH ARGUMENTS nARGUMENTS)
        if (nARGUMENTS GREATER 0)
            message(FATAL_ERROR "NO_DEFAULT_RUN and ARGUMENTS are incomaptible for test ${TEST_NAME}")
        endif (nARGUMENTS GREATER 0)
    endif(nNO_DEFAULT_RUN GREATER 1)

    if (nSCRIPT EQUAL 0)
        add_executable(${TEST_NAME} ${SOURCES})
        target_link_libraries(${TEST_NAME} ${LIBRARIES})
        list(INSERT TARGET_DEPENDENCIES 0 ${TEST_NAME})
        list(LENGTH TARGET_DEPENDENCIES nTARGET_DEPENDENCIES)
    endif (nSCRIPT EQUAL 0)

    if(NOT TEST_NAME)
        message(FATAL_ERROR "could not find a test name")
        endif(NOT TEST_NAME)
    # even script tests may define target dependencies, even though it's
    # probably of little use.
    add_custom_target(${TEST_NAME}_dependencies)
    add_dependencies(all_test_dependencies ${TEST_NAME}_dependencies)

    if(nTARGET_DEPENDENCIES GREATER 0)
        add_dependencies(${TEST_NAME}_dependencies ${TARGET_DEPENDENCIES})
        add_test(builddep_${TEST_NAME} ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${TEST_NAME}_dependencies)
        # set(${TEST_NAME}_has_build_depends 1)
    endif(nTARGET_DEPENDENCIES GREATER 0)

    if (nSCRIPT EQUAL 0)
        set(${TEST_NAME}_EXECUTABLE ${TEST_NAME})
    else (nSCRIPT EQUAL 0)
        set(binary ${SCRIPT})
        if (${SCRIPT} MATCHES ".sh$")
            set(binary env bash ${SCRIPT})
        endif (${SCRIPT} MATCHES ".sh$")
        set(${TEST_NAME}_EXECUTABLE ${binary})
    endif(nSCRIPT EQUAL 0)

    if(NOT NO_DEFAULT_RUN)
        add_test (${TEST_NAME} ${${TEST_NAME}_EXECUTABLE} ${ARGUMENTS})
        set(${TEST_NAME}_ARGUMENTS ${ARGUMENTS})
        set(${TEST_NAME}_TEST_EXISTS 1)
        if(nTARGET_DEPENDENCIES GREATER 0)
            # seems that we need set_tests_properties to take only a
            # ;-separated list here.
            set(deps builddep_${TEST_NAME} ${TEST_DEPENDENCIES})
            string_join(all_deps ";" ${deps})
            set_tests_properties (${TEST_NAME} PROPERTIES DEPENDS "${all_deps}")
        endif(nTARGET_DEPENDENCIES GREATER 0)
    endif(NOT NO_DEFAULT_RUN)
endmacro(cado_define_test)

macro(cado_divert_test TEST_BASE DIVERSION)
    # This macro is appropriate when a test defined by
    # cado_define_test (possibly with NO_DEFAULT_RUN) needs to be
    # called several times with various parameter sets.
    # binary. The different ways to call this macro are:
    # 
    # optional arguments:
    #           TEST_DEPENDENCIES
    #           TARGET_DEPENDENCIES
    #           ARGUMENTS
    #           APPEND_ARGUMENTS
    #
    # example use:
    #
    #     cado_define_test(TEST_NAME foo .....)
    #
    #     cado_divert_test(foo 1 --fast)
    #     cado_divert_test(foo 2 --thorough)
    #     cado_divert_test(foo extra --fallback TEST_DEPENDENCIES test3)
    #     cado_divert_test(foo extra TEST_DEPENDENCIES test3
    #                                ARGUMENTS --fallback)
    #
    # this defines tests foo_1, foo_2, and foo_extra, with the specific
    # arguments given. The two last syntaxes are equivalent.
    #
    # NOTE that if the base test name does indeed correspond to a test
    # (that is, it was not defined with NO_DEFAULT_RUN), then the
    # diverted test is made to depend on this base test.
    set(ARGUMENTS)
    set(APPEND_ARGUMENTS)
    set(TEST_DEPENDENCIES)
    set(TARGET_DEPENDENCIES)
    set(current ARGUMENTS)
    foreach(x ${ARGN})
        if (17 EQUAL 42)        # cmake syntax sucks
            # dead branch, I know. But since the 1st test has privileged
            # status in cmake syntax, let's make it a straw man.
        elseif (x STREQUAL "APPEND_ARGUMENTS")
            SET(current APPEND_ARGUMENTS)
        elseif (x STREQUAL "TEST_DEPENDENCIES")
            SET(current TEST_DEPENDENCIES)
        elseif (x STREQUAL "TARGET_DEPENDENCIES")
            SET(current TARGET_DEPENDENCIES)
        else(17 EQUAL 42)       # cmake syntax sucks
            list(APPEND ${current} ${x})
        endif(17 EQUAL 42)      # cmake syntax sucks
    endforeach(x)
    # do some checking.
    list(LENGTH TARGET_DEPENDENCIES nTARGET_DEPENDENCIES)
    list(LENGTH ARGUMENTS nARGUMENTS)
    list(LENGTH APPEND_ARGUMENTS nAPPEND_ARGUMENTS)
    if ((nAPPEND_ARGUMENTS GREATER 0) AND (nARGUMENTS GREATER 0))
        message(FATAL_ERROR "ARGUMENTS and APPEND_ARGUMENTS are incompatible")
    endif ((nAPPEND_ARGUMENTS GREATER 0) AND (nARGUMENTS GREATER 0))

    set(TEST_NAME "${TEST_BASE}_${DIVERSION}")
    if (nAPPEND_ARGUMENTS GREATER 0)
        set(ARGUMENTS ${${TEST_BASE}_ARGUMENTS} ${APPEND_ARGUMENTS})
        list(LENGTH ARGUMENTS nARGUMENTS)
    endif (nAPPEND_ARGUMENTS GREATER 0)
    add_test(${TEST_NAME} ${${TEST_BASE}_EXECUTABLE} ${ARGUMENTS})

    # we want to define builddep_X unconditionally, depending at least on
    # the dependencies of the base test.
    add_custom_target(${TEST_NAME}_dependencies)
    add_test(builddep_${TEST_NAME} ${CMAKE_COMMAND} --build
        ${CMAKE_BINARY_DIR} --target ${TEST_NAME}_dependencies)

    add_dependencies(${TEST_NAME}_dependencies ${TEST_BASE}_dependencies ${TARGET_DEPENDENCIES})
    # set("${TEST_NAME}_has_build_depends" 1)

    # prepare a ;-separated list for set_tests_properties
    # note that we do not need to make the test depend on
    # builddep_${TEST_BASE}, since we have made ${TEST_NAME}_dependencies
    # depend on ${TEST_BASE}_dependencies. It's in fact relieving, given
    # that builddep_${TEST_BASE} does not necessarily exist.
    set(deps builddep_${TEST_NAME} ${TEST_DEPENDENCIES})
    if(${TEST_BASE}_TEST_EXISTS)
        list(APPEND deps ${TEST_BASE})
    endif(${TEST_BASE}_TEST_EXISTS)
    string_join(all_deps ";" ${deps})
    set_tests_properties (${TEST_NAME} PROPERTIES DEPENDS "${all_deps}")
endmacro(cado_divert_test TEST_BASE DIVERSION ARGUMENTS)



if(DEFINED ENV{CHECKS_EXPENSIVE})
    cado_define_test(full_c59
        SCRIPT
            ${CADO_NFS_SOURCE_DIR}/factor.sh
            90377629292003121684002147101760858109247336549001090677693 -t 2
        TARGET_DEPENDENCIES
            polyselect2l polyselect_ropt
            freerel makefb
            las
            dup1 dup2 purge merge replay
            characters
            bwc_full
            sqrt)
    find_program(MAGMA magma)
    if (MAGMA)
        cado_define_test(full_p59
            SCRIPT
                ${CMAKE_CURRENT_SOURCE_DIR}/test_full_p59.sh
                ${CADO_NFS_SOURCE_DIR}
            TARGET_DEPENDENCIES full_c59_dependencies debug_renumber
            plingen_pz sm reconstructlog-dl merge-dl replay-dl
        )
        cado_define_test(full_p30
            SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/test_full_p30.sh
            ARGUMENTS ${CADO_NFS_SOURCE_DIR}
            TARGET_DEPENDENCIES full_c59_dependencies debug_renumber
            sm_simple las_descent)
        cado_define_test(test_gfp3
            SCRIPT
                ${CMAKE_CURRENT_SOURCE_DIR}/test_gfp3.sh
                ${CADO_NFS_SOURCE_DIR} ${CADO_NFS_BINARY_DIR}
            TARGET_DEPENDENCIES full_p59_dependencies)
        cado_define_test(test_gfp2
            SCRIPT
                ${CMAKE_CURRENT_SOURCE_DIR}/test_gfp2.sh
                ${CADO_NFS_SOURCE_DIR} ${CADO_NFS_BINARY_DIR}
            TARGET_DEPENDENCIES full_p59_dependencies filter_galois
            polyselect_gfpn)
    endif(MAGMA)
endif(DEFINED ENV{CHECKS_EXPENSIVE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library (tests STATIC test_iter.c tests_common.c)
target_link_libraries (tests ${gmp_libname})

cado_define_test(test_tests_common.c NO_DEFAULT_RUN LIBRARIES tests)
cado_divert_test(test_tests_common 1)
cado_divert_test(test_tests_common 2 -seed)
cado_divert_test(test_tests_common 3 -seed 1)
set_tests_properties(test_tests_common_3 PROPERTIES PASS_REGULAR_EXPRESSION "Using random seed=1")
cado_divert_test(test_tests_common 4 -seed a)
cado_divert_test(test_tests_common 5 -seed "")
set_tests_properties(test_tests_common_2 test_tests_common_4 test_tests_common_5 PROPERTIES WILL_FAIL 1)
cado_divert_test(test_tests_common 6 -iter 10)
set_tests_properties(test_tests_common_6 PROPERTIES PASS_REGULAR_EXPRESSION "Using 10 iterations")

cado_define_test(test_iceildiv.c)

add_subdirectory (utils)
add_subdirectory (sieve)
add_subdirectory (polyselect)
add_subdirectory (filter)
add_subdirectory (linalg)
add_subdirectory (gf2x)
