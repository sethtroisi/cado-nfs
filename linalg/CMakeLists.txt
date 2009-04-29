
set_source_files_properties (gauss.c COMPILE_FLAGS "-DMULTI_ROW=3 -DNO_MAIN")

add_executable(freerel freerel.c ../sieve/fb.c hashpair.c)
target_link_libraries(freerel utils m)

add_executable(duplicates duplicates.c hashpair.c files.c)
target_link_libraries(duplicates utils)

add_executable(purge purge.c hashpair.c files.c)
target_link_libraries(purge utils m)

add_executable(merge merge.c markowitz.c swar.c mst.c report.c dclist.c merge_mono.c sparse_mat.c sparse.c)
target_link_libraries(merge utils m)

add_executable(replay replay.c sparse.c)
target_link_libraries(replay utils)

# checkdep: checkdep.o ;
# linalg: linalg.o gauss.o readmat.o ; target_link_libraries( utils
# mkbitstrings: mkbitstrings.c ; target_link_libraries( utils
# hismerge: hismerge.c ; target_link_libraries( utils
# gauss.c: gauss.c ; $(CC) $(MY_CFLAGS) -DMULTI_ROW=3 -DNO_MAIN $< -c -o $@

add_executable(characters characters.c files.c gauss.c)
target_link_libraries(characters utils)

add_executable(allsqrt allsqrt.c files.c hashpair.c)
target_link_libraries(allsqrt utils)

add_executable(apply_perm apply_perm.c)
target_link_libraries(apply_perm utils)


# This library is used by bwc mostly, but since its files are relevant here as
# well...
add_library(mat_io readbuffer.c readmat.c readmat-easy.c)

add_executable(balance balance.c readbuffer.c rowset_heap.cpp)
target_link_libraries(balance utils m)

add_executable(transpose transpose.c readmat.c)
target_link_libraries(transpose utils m)

