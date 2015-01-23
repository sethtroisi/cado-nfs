
include(${PROJECT_SOURCE_DIR}/config/search_for_function.cmake)
set(CMAKE_REQUIRED_LIBRARIES ${gmp_libname})
search_for_function(mpn_xor_n HAVE_MPN_XOR_N)
if(NOT HAVE_MPN_XOR_N)
search_for_function(__gmpn_xor_n HAVE_MPN_XOR_N)
endif(NOT HAVE_MPN_XOR_N)
