
find_package(OpenMP)
if (OPENMP_FOUND)
    set (HAVE_OPENMP 1)
    # Don't set unconditionally !
    # set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
