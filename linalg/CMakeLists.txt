
# checkdep: checkdep.o ;
# linalg: linalg.o gauss.o readmat.o ; target_link_libraries( utils
# mkbitstrings: mkbitstrings.c ; target_link_libraries( utils
# hismerge: hismerge.c ; target_link_libraries( utils
# gauss.c: gauss.c ; $(CC) $(MY_CFLAGS) -DMULTI_ROW=3 -DNO_MAIN $< -c -o $@

add_executable(characters characters.c gauss.c)
target_link_libraries(characters utils)

add_executable(apply_perm apply_perm.c)
target_link_libraries(apply_perm utils)

add_executable(mkbitstrings mkbitstrings.c)
target_link_libraries(mkbitstrings utils)

# This library is used by bwc mostly, but since its files are relevant here as
# well...
add_library(mat_io readbuffer.c readmat.c readmat-easy.c)

add_executable(balance balance.c readbuffer.c rowset_heap.cpp)
target_link_libraries(balance utils m)

add_executable(transpose transpose.c readmat.c)
target_link_libraries(transpose utils m)


install(TARGETS mkbitstrings characters apply_perm balance transpose RUNTIME DESTINATION bin/linalg)

set_source_files_properties (gauss.c COMPILE_FLAGS "-DMULTI_ROW=3 -DNO_MAIN")

include(CMakeLists-nodist.txt OPTIONAL)
