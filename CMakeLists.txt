project(CADO_NFS)

# I have 2.6 -- yet I recall having used 2.4 before, and so far nothing that I
# now use appears to be so new that I expect it to cause trouble with 2.4. So
# my guess is that 2.4 is good enough. If it isn't, then so be it.
#
# XXX Maybe set_property ? (utils/)
cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
        cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)



## You can force a path to gmp.h using -DGMP_INCLUDE_PATH=/path/to/gmph
find_path(GMP_INCLUDE_PATH gmp.h
    /usr/include
    /usr/local/include)

## You can force the gmp lib using -DGMP_LIB=/opt/lib/libgmp.a
find_library(GMP_LIB gmp
    /usr/lib
    /usr/lib64
    /usr/local/lib
    /usr/local/lib64)

## You can override these defaults with -DCMAKE_C_FLAGS="-g -fomit-frame-pointer"

execute_process(COMMAND ./version.sh
        WORKING_DIRECTORY ${CADO_NFS_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE revision)
message("Configuring CADO version ${revision}")
if(CMAKE_COMPILER_IS_GNUCC)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g -W -Wall")
endif(CMAKE_COMPILER_IS_GNUCC)

# TODO ; set the COMPILE_DEFINITIONS property in the global scope ; see
# cmake-properties.txt
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DREV=\"\\\"${revision}\\\"\"")

include_directories (${CADO_NFS_SOURCE_DIR} ${GMP_INCLUDE_PATH}
    ${CADO_NFS_SOURCE_DIR}/utils)
link_directories (${CADO_NFS_BINARY_DIR}/utils)
link_directories (${CADO_NFS_BINARY_DIR}/gf2x)

add_subdirectory (utils)
add_subdirectory (polyselect)
add_subdirectory (sieve)
add_subdirectory (sqrt/naive)
add_subdirectory (cantor)
add_subdirectory (linalg)
add_subdirectory (linalg/bwc)
add_subdirectory (linalg/bwc/lingen)

add_custom_command(OUTPUT ${CADO_NFS_BINARY_DIR}/gf2x/Makefile
    COMMAND mkdir -p ${CADO_NFS_BINARY_DIR}/gf2x
    COMMAND cd gf2x && ${CADO_NFS_SOURCE_DIR}/gf2x/configure --disable-shared
    COMMENT "Configuring gf2x"
)

add_custom_target(gf2x ALL
    DEPENDS ${CADO_NFS_BINARY_DIR}/gf2x/Makefile
    WORKING_DIRECTORY gf2x
    COMMAND make
    COMMENT "Building gf2x"
)
