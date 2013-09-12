# This is a variation around CHECK_FUNCTION_EXISTS.

# search_for_function(FUNCTION VARIABLE [lib1 [lib2]])
# - macro which checks if the function exists
#  FUNCTION - the name of the function
#  VARIABLE - variable to store the result

# The following variables may be set before calling this macro to
# modify the way the check is run:

#  CMAKE_REQUIRED_FLAGS = string of compile command line flags
#  CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
#  CMAKE_REQUIRED_INCLUDES = list of include directories
#  CMAKE_REQUIRED_LIBRARIES = list of libraries to link

# Additionally, extra libraries may be suggested, and will be tried in
# order as supplementary arguments to CMAKE_REQUIRED_LIBRARIES -- note
# that each library that has been found will be tried and kept for later
# uses.

# If successful, VARIABLE is set to 1 and CMAKE_REQUIRED_LIBRARIES
# contains the winning library set. CMAKE_REQUIRED_LIBRARIES_EXTRA
# contains the libraries which had to be added to the initial library
# set.

macro(search_for_function FUNCTION VARIABLE)
    message(STATUS "Looking for ${FUNCTION}")
    set(sff_rl ${CMAKE_REQUIRED_LIBRARIES})
    set(atwork 1)
    foreach(extra "" ${ARGN})
        set(testing_this ${atwork})
        if(extra AND atwork)
            find_library(testing_this ${extra})
        endif(extra AND atwork)
        if(testing_this)
            set(sff_rl ${sff_rl} ${extra})
            set(MACRO_CHECK_FUNCTION_DEFINITIONS
              "-DCHECK_FUNCTION_EXISTS=${FUNCTION} ${CMAKE_REQUIRED_FLAGS}")
            if(sff_rl)
              set(sff_ADD_LIBRARIES "-DLINK_LIBRARIES:STRING=${sff_rl}")
            else(sff_rl)
              set(sff_ADD_LIBRARIES)
            endif(sff_rl)
            if(CMAKE_REQUIRED_INCLUDES)
              set(sff_ADD_INCLUDES
                "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}")
            else(CMAKE_REQUIRED_INCLUDES)
              set(sff_ADD_INCLUDES)
            endif(CMAKE_REQUIRED_INCLUDES)
            # message(STATUS "Trying to compile with ${CMAKE_REQUIRED_DEFINITIONS} -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS} and also ${sff_ADD_LIBRARIES} ${sff_ADD_INCLUDES}")
            try_compile(${VARIABLE}
              ${CMAKE_BINARY_DIR}
              ${CMAKE_ROOT}/Modules/CheckFunctionExists.c
              COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
              CMAKE_FLAGS
              -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
              "${sff_ADD_LIBRARIES}"
              "${sff_ADD_INCLUDES}"
              OUTPUT_VARIABLE OUTPUT)

            if(${VARIABLE})
              set(${VARIABLE} 1 CACHE INTERNAL "Have function ${FUNCTION}")
              if(sff_rl)
                  message(STATUS "Looking for ${FUNCTION} - requires ${sff_rl}")
              else(sff_rl)
                  message(STATUS "Looking for ${FUNCTION} - found")
              endif(sff_rl)
              file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
                "Determining if the function ${FUNCTION} exists passed with the following output:\n"
                "${OUTPUT}\n\n")
              set(CMAKE_REQUIRED_LIBRARIES ${sff_rl})
              set(CMAKE_REQUIRED_LIBRARIES_EXTRA ${extra})
              set(atwork 0)
            else(${VARIABLE})
              file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
                "Determining if the function ${FUNCTION} exists failed with the following output:\n"
                "${OUTPUT}\n\n")
            endif(${VARIABLE})
         endif(testing_this)
    endforeach(extra)

    if(NOT ${VARIABLE})
          message(STATUS "Looking for ${FUNCTION} - not found")
          set(${VARIABLE} "" CACHE INTERNAL "Have function ${FUNCTION}")
    endif(NOT ${VARIABLE})
endmacro(search_for_function FUNCTION VARIABLE)
