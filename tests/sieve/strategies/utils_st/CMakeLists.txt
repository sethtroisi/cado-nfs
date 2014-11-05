set (SUF_SRC ${CADO_NFS_SOURCE_DIR}/sieve/strategies/utils_st/)

include_directories(${CADO_NFS_SOURCE_DIR}/sieve/strategies/utils_st)
include_directories(${CADO_NFS_SOURCE_DIR})
include_directories(${CADO_NFS_SOURCE_DIR}/sieve/ecm)
include_directories(${CADO_NFS_SOURCE_DIR}/sieve)

add_executable(test_fm test_fm.c ${SUF_SRC}fm.c)
#build_test (test_fm test_fm.c ${CADO_NFS_SOURCE_DIR}/sieve/strategies/utils_st/fm.c "")
run_test (test_fm test_fm)

add_executable(test_strategy test_strategy.c ${SUF_SRC}tab_strategy.c ${SUF_SRC}fm.c ${SUF_SRC}strategy.c ${SUF_SRC}fm.c ${SUF_SRC}tab_fm.c)
run_test (test_strategy test_strategy)
