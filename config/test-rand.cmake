# test if rand() gives true alea
message(STATUS "Testing if rand() gives true alea")
try_run(test-rand_runs test-rand_compiles
        ${CADO_NFS_BINARY_DIR}/config
        ${CADO_NFS_SOURCE_DIR}/config/test-rand.c)
if(test-rand_runs EQUAL 0)
   set(HAVE_RAND_BUG 0)
   message(STATUS "Testing if rand() gives true alea -- Probably")
else(test-rand_runs EQUAL 0)
   set(HAVE_RAND_BUG 1)
   message(STATUS "Testing if rand() gives true alea -- Buggy")
endif(test-rand_runs EQUAL 0)
