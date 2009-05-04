add_executable(allsqrt allsqrt.c)
target_link_libraries(allsqrt utils)

install(TARGETS allsqrt RUNTIME DESTINATION bin/sqrt)

add_executable(algsqrt algsqrt.c poly.c)
target_link_libraries(algsqrt utils m)

add_executable(checkdep checkdep.c)
target_link_libraries(checkdep utils)

install(TARGETS algsqrt allsqrt checkdep RUNTIME DESTINATION bin/sqrt)

include(CMakeLists-nodist.txt OPTIONAL)
