add_executable(sqrt sqrt.c)
target_link_libraries(sqrt utils m)

install(TARGETS sqrt RUNTIME DESTINATION bin/sqrt)

add_executable(checkdep checkdep.c)
target_link_libraries(checkdep utils)

install(TARGETS sqrt checkdep RUNTIME DESTINATION bin/sqrt)

add_executable(crtalgsqrt crtalgsqrt.c powers_of_p.cpp polyroots.c)
target_link_libraries(crtalgsqrt utils m)

set_target_properties(crtalgsqrt PROPERTIES
    COMPILE_FLAGS "--mpi" LINK_FLAGS "--mpi")

include(CMakeLists-nodist.txt OPTIONAL)
