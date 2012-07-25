# test a bug appeared in Ubuntu/Linaro 4.6.3-1ubuntu5
message(STATUS "Testing known bugs for compiler")
try_run(gcc-ubuntu-bug_runs gcc-ubuntu-bug_compiles
            ${CADO_NFS_BINARY_DIR}/config
            ${CADO_NFS_SOURCE_DIR}/config/gcc-ubuntu-bug.c)
if(gcc-ubuntu-bug_compiles)
    if (gcc-ubuntu-bug_runs EQUAL 0)
        message(STATUS "Testing known bugs for compiler -- Not found")
        set(VOLATILE_IF_GCC_UBUNTU_BUG 0)
    else (gcc-ubuntu-bug_runs EQUAL 0)
        message(STATUS "Testing known bugs for compiler -- Found a bug")
        set(VOLATILE_IF_GCC_UBUNTU_BUG 1)
    endif (gcc-ubuntu-bug_runs EQUAL 0)
else(gcc-ubuntu-bug_compiles)
  message(STATUS "Testing known bugs for compiler -- Error (cannot compile)")
  set(VOLATILE_IF_GCC_UBUNTU_BUG 0)
endif(gcc-ubuntu-bug_compiles)
