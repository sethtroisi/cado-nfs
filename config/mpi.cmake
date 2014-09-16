
#############################################################
# mpi
# Don't use the FindMPI module, it's buggy.

# I can't make if($ENV{MPI}) evaluate to true as I want. In particular, I
# want all of these yield true:
# MPI=1 , MPI=on, MPI=yes, etc.
# MPI=/opt/openmpi-1.7/
# MPI=openmpi
# The following excerpt from the doc seems to be just plain wrong:
# #         if(variable)
# #
# #       True if the variable's value is not empty, 0, N, NO, OFF, FALSE,
# #       NOTFOUND, or <variable>-NOTFOUND.
if("$ENV{MPI}" MATCHES "^(0|NO|no|OFF|off|)$")
    message(STATUS "MPI is not enabled")
    set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
    set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
else("$ENV{MPI}" MATCHES "^(0|NO|no|OFF|off|)$")
    set(findprog_flags)
    set(mpicc_names)
    set(mpicxx_names)
    set(mpiexec_names)
    if("$ENV{MPI}" MATCHES "^(1|YES|yes|ON|on|)$")
        list(APPEND mpicc_names "mpicc")
        list(APPEND mpicxx_names "mpic++" "mpicxx" "mpiCC")
        list(APPEND mpiexec_names "mpiexec")
    else("$ENV{MPI}" MATCHES "^(1|YES|yes|ON|on|)$")
        if("$ENV{MPI}" MATCHES "/")
            # If MPI contains a /, then we assume it should be a path
            list(APPEND findprog_flags
                HINTS "$ENV{MPI}" "$ENV{MPI}/bin"
                NO_DEFAULT_PATH
                NO_CMAKE_ENVIRONMENT_PATH
                NO_CMAKE_PATH
                NO_SYSTEM_ENVIRONMENT_PATH
                NO_CMAKE_SYSTEM_PATH)
            # These are the standard names.
            list(APPEND mpicc_names "mpicc")
            list(APPEND mpicxx_names "mpic++" "mpicxx" "mpiCC")
            list(APPEND mpiexec_names "mpiexec")
        else("$ENV{MPI}" MATCHES "/")
            # otherwise we make the .<variant> binary names higher
            # priority than others This is for finding things such as
            # mpicc.mpich2 which get installed by the alternatives
            # mechanism on debian-like systems.
            list(APPEND mpicc_names "mpicc.${MPI}")
            list(APPEND mpicxx_names "mpic++.${MPI}" "mpicxx.${MPI}" "mpiCC.${MPI}")
            list(APPEND mpiexec_names "mpiexec.${MPI}")
            # Well. Presently we're in fact *not* pushing the standard
            # names in the search list. Should we ?
        endif("$ENV{MPI}" MATCHES "/")
    endif("$ENV{MPI}" MATCHES "^(1|YES|yes|ON|on|)$")

    find_program(MPI_C_COMPILER ${mpicc_names} ${findprog_flags})
    find_program(MPI_CXX_COMPILER ${mpicxx_names} ${findprog_flags})
    find_program(MPIEXEC ${mpiexec_names} ${findprog_flags})

    if (MPI_C_COMPILER AND MPI_CXX_COMPILER AND MPIEXEC)
        message(STATUS "Using MPI C compiler ${MPI_C_COMPILER}")
        message(STATUS "Using MPI C++ compiler ${MPI_CXX_COMPILER}")
        message(STATUS "Using MPI driver ${MPIEXEC}")
        get_filename_component(HAVE_MPI ${MPIEXEC} PATH)
        # We're using this variable in the top-level substitution, so it needs
        # to escape its scope and go into the cache right now.
        set(WITH_MPI 1 CACHE INTERNAL "MPI is being used (for relevant code parts)")

        # Now check for the MPI API version.
        macro(my_try_compile_mpicc SOURCE VAR)
            string(RANDOM LENGTH 8 ALPHABET "0123456789abcdef" uuid)
            string(TOLOWER "${VAR}" lcvar)
            set(checkfile
                "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/check_${lcvar}_${UUID}.c")

            FILE(WRITE "${checkfile}" "${SOURCE}\n")
            MESSAGE(STATUS "Performing Test ${VAR}")
            execute_process(
                COMMAND ${MPI_C_COMPILER} -c ${checkfile}
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                RESULT_VARIABLE test_return_code
                OUTPUT_VARIABLE test_stdout
                ERROR_VARIABLE test_stderr
                )
            FILE(APPEND
                ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
                "Summary for ${VAR} check:\n"
                "Source file was ${checkfile}\n"
                "Return code: ${test_return_code}\n"
                "stdout: ${test_stdout}\n"
                "stderr: ${test_stderr}\n"
                "\n"
                )
            if (test_return_code)
                SET(${VAR} 0 CACHE INTERNAL "Test ${VAR}")
                MESSAGE(STATUS "Performing Test ${VAR} -- Failed")
            else (test_return_code)
                SET(${VAR} 1 CACHE INTERNAL "Test ${VAR}")
                MESSAGE(STATUS "Performing Test ${VAR} -- Success")
            endif (test_return_code)
        endmacro(my_try_compile_mpicc SOURCE VAR)

        my_try_compile_mpicc("
        #include <mpi.h>
        #if !((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1))
        #error \"MPI version 2.1 not supported\"
        #endif
        int main(int argc, char** argv)
        {
        MPI_Init(&argc,&argv);
        MPI_Finalize();
        }
        " HAVE_MPI2_API)

        my_try_compile_mpicc("
        #include <mpi.h>
        #if !((MPI_VERSION > 3) || (MPI_VERSION == 3 && MPI_SUBVERSION >= 0))
        #error \"MPI version 3.0 not supported\"
        #endif
        int main(int argc, char** argv)
        {
        MPI_Init(&argc,&argv);
        MPI_Finalize();
        }
        " HAVE_MPI3_API)

    else(MPI_C_COMPILER AND MPI_CXX_COMPILER AND MPIEXEC)
        message(FATAL_ERROR "Cannot find all of mpicc/mpic++/mpiexec with MPI=$ENV{MPI}")
    endif(MPI_C_COMPILER AND MPI_CXX_COMPILER AND MPIEXEC)
endif("$ENV{MPI}" MATCHES "^(0|NO|no|OFF|off|)$")
