# test endianness
message(STATUS "Testing endianness")
try_run(endian_runs endian_compiles
        ${CADO_NFS_BINARY_DIR}/config
	${CADO_NFS_SOURCE_DIR}/config/endian.c)
if (endian_compiles)
   if (endian_runs EQUAL 0)
      message(STATUS "System is little endian")
      set(IS_BIG_ENDIAN 0)
   else (endian_runs EQUAL 0)
      message(STATUS "System is big endian")
      set(IS_BIG_ENDIAN 1)
   endif (endian_runs EQUAL 0)
else(endian_compiles)
endif(endian_compiles)
            
