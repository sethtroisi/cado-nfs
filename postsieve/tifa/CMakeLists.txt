SET (algo_sources 
    lib/algo/src/cfrac.c
    lib/algo/src/fermat.c
    lib/algo/src/qs.c
    lib/algo/src/siqs.c
    lib/algo/src/squfof.c
    lib/algo/src/tdiv.c
    lib/algo/src/tifa_factor.c)

SET (utils_sources
    lib/utils/src/array.c
    lib/utils/src/bernsteinisms.c
    lib/utils/src/factoring_machine.c
    lib/utils/src/funcs.c
    lib/utils/src/gauss_elim.c
    lib/utils/src/gmp_utils.c
    lib/utils/src/hashtable.c
    lib/utils/src/lindep.c
    lib/utils/src/linked_list.c
    lib/utils/src/matrix.c
    lib/utils/src/res_tdiv.c
    lib/utils/src/smooth_filter.c
    lib/utils/src/sqrt_cont_frac.c
    lib/utils/src/x_array_list.c
    lib/utils/src/x_tree.c)

SET (data_sources
    lib/data/src/first_primes.c)

SET (sources
    ${algo_sources} ${utils_sources} ${data_sources})

## Add -std=c99 and Revision stuff to compile flags for this target
SET (CFLAGS "-std=c99 -DREV=\\\"`svnversion ${CMAKE_CURRENT_SOURCE_DIR}`\\\"")
SET_SOURCE_FILES_PROPERTIES (${sources}
    COMPILE_FLAGS ${CFLAGS})

INCLUDE_DIRECTORIES (. lib lib/utils/include lib/algo/include lib/data/include)

ADD_LIBRARY (tifa ${sources})

