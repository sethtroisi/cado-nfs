add_executable(sqrt sqrt.c)
target_link_libraries(sqrt utils m)

install(TARGETS sqrt RUNTIME DESTINATION bin/sqrt)

add_executable(checkdep checkdep.c)
target_link_libraries(checkdep utils)

install(TARGETS sqrt checkdep RUNTIME DESTINATION bin/sqrt)

include(CMakeLists-nodist.txt OPTIONAL)
