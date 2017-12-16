
# The Jevents library is useful only for developers. We make little effort
# here of finding it. The developers can figure it out by themselves.

message(STATUS "Testing whether Jevents library can be used")
try_run(jevents_runs jevents_compiles
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/jevents.c
    COMPILE_OUTPUT_VARIABLE jevents_compile_output
    RUN_OUTPUT_VARIABLE jevents_run_output
    LINK_LIBRARIES jevents
)

if(jevents_compiles)
    if (jevents_runs MATCHES FAILED_TO_RUN)
        message(STATUS "Testing whether Jevents library can be used -- No")
	# message(STATUS "${jevents_run_output}")
        set (HAVE_JEVENTS 0)
    else()
        message(STATUS "Testing whether Jevents library can be used -- Yes")
        set (HAVE_JEVENTS 1)
    endif()
else()
    message(STATUS "Testing whether Jevents library can be used -- No")
    # message(STATUS "${jevents_compile_output}")
    set (HAVE_JEVENTS 0)
endif()
