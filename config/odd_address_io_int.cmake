# An int could be load & write from an odd adress ?
message(STATUS "Testing if an int could be load & write from an odd adress")
try_run(odd_address_io_int_runs odd_address_io_int_compiles
        ${CADO_NFS_BINARY_DIR}/config
        ${CADO_NFS_SOURCE_DIR}/config/odd_address_io_int.c)
if(odd_address_io_int_runs EQUAL 0)
   set(ODD_ADRESS_IO_INT 1)
   message(STATUS "Testing if an int could be load & write from an odd adress -- Yes")
else(odd_address_io_int_runs EQUAL 0)
   set(ODD_ADRESS_IO_INT 0)
   message(STATUS "Testing if an int could be load & write from an odd adress -- No")
endif(odd_address_io_int_runs EQUAL 0)
