
project(CADO_NFS)

# May be overridden by environment.
set(CADO_VERSION "0.1-rc2")

#####################################################################
# cmake version -- I have 2.6 -- yet I recall having used 2.4 before, and so
# far nothing that I now use appears to be so new that I expect it to cause
# trouble with 2.4. So my guess is that 2.4 is good enough. If it isn't, then
# so be it.
#
# XXX Maybe set_property ? (utils/)
cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
        cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#####################################################################
# Get some overrides from the environment.

macro(get_override extname intname description)
string(COMPARE NOTEQUAL "$ENV{${extname}}" "" HAS_${intname}_OVERRIDE)
if(HAS_${intname}_OVERRIDE)
    set(${intname} "$ENV{${extname}}")
    message("${description} set to ${${intname}}")
endif(HAS_${intname}_OVERRIDE)
endmacro(get_override varname description)

# Syntax is (user-friendly name) (cmake name) (description)
# Make sure everything is exported from scripts/call_cmake.sh
get_override(PREFIX     CMAKE_INSTALL_PREFIX  "Install prefix")
get_override(CFLAGS     CMAKE_C_FLAGS         "CFLAGS")
get_override(CXXFLAGS   CMAKE_CXX_FLAGS       "CXXFLAGS")
get_override(CC         CMAKE_C_COMPILER      "CC")
get_override(CXX        CMAKE_CXX_COMPILER    "CXX")
get_override(CADO_VERSION    CADO_VERSION               "package version")

#####################################################################
# What are we building ?

execute_process(COMMAND scripts/version.sh
        WORKING_DIRECTORY ${CADO_NFS_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE CADO_REV)
message("Configuring CADO revision ${CADO_REV}")

#####################################################################
# Add the CFLAGS that are needed in all cases.

if(CMAKE_COMPILER_IS_GNUCC)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g -W -Wall")
endif(CMAKE_COMPILER_IS_GNUCC)

#####################################################################
# Set dependencies for external software.
# XXX FIXME !

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

#####################################################################
# Obtain information from the host.

include(CheckFunctionExists)

# posix_memalign
CHECK_FUNCTION_EXISTS(posix_memalign      HAVE_POSIX_MEMALIGN)

# SSE-2
message(STATUS "Testing whether sse-2 code can be used")
try_compile(sse2_compiles
            ${CADO_NFS_BINARY_DIR}/config
            ${CADO_NFS_SOURCE_DIR}/config/sse2.c)
if(sse2_compiles)
    message(STATUS "Testing whether sse-2 code can be used -- Yes")
else(sse2_compiles)
    try_compile(sse2_compiles
        ${CADO_NFS_BINARY_DIR}/config
        ${CADO_NFS_SOURCE_DIR}/config/sse2.c
        COMPILE_DEFINITIONS -msse2)
    if(sse2_compiles)
        message(STATUS "Trying whether sse-2 code can be used -- Yes, with -msse2")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    else(sse2_compiles)
        message(STATUS "Trying whether sse-2 code can be used -- No")
    endif(sse2_compiles)
endif(sse2_compiles)

# gcc-style x86_64 assembly.
message(STATUS "Trying to compile gcc-style x86_64 assembly")
try_compile(HAVE_GCC_STYLE_AMD64_ASM
    ${CADO_NFS_BINARY_DIR}/config
    ${CADO_NFS_SOURCE_DIR}/config/asm_sample.c
    COMPILE_DEFINITIONS "-x assembler")
if(HAVE_GCC_STYLE_AMD64_ASM)
    message(STATUS "Trying to compile gcc-style x86_64 assembly -- Success")
else(HAVE_GCC_STYLE_AMD64_ASM)
    message(STATUS "Trying to compile gcc-style x86_64 assembly -- Failed")
endif(HAVE_GCC_STYLE_AMD64_ASM)

#####################################################################
# now the real stuff.

# First, create the config file gathering everything we've got so far.
configure_file(cado_config_h.in cado_config.h ESCAPE_QUOTES @ONLY)

include_directories (${GMP_INCLUDE_PATH})
include_directories (${CADO_NFS_SOURCE_DIR})
include_directories (${CADO_NFS_SOURCE_DIR}/utils)
include_directories (${CADO_NFS_BINARY_DIR})

link_directories (${CADO_NFS_BINARY_DIR}/utils)
link_directories (${CADO_NFS_BINARY_DIR}/gf2x)

# Start with cantor since it so easily fails when -msse2 has accidentally been
# forgotten on a P4 build.
add_subdirectory (cantor)
add_subdirectory (utils)
add_subdirectory (polyselect)
add_subdirectory (sieve)
add_subdirectory (sqrt)
add_subdirectory (linalg)
add_subdirectory (merge)
add_subdirectory (linalg/bwc)

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

set(package_basename cado-nfs-${CADO_VERSION})

add_custom_target(dist
    COMMAND scripts/dist.sh ${package_basename}
    WORKING_DIRECTORY ${CADO_NFS_SOURCE_DIR}
)

install(PROGRAMS cadofactor.pl DESTINATION bin)

set(example_subdir share/${package_basename}/examples)

install(FILES params/params.c59 DESTINATION ${example_subdir})
install(FILES params/params.c79 DESTINATION ${example_subdir})

configure_file(
    ${CADO_NFS_SOURCE_DIR}/run_example.sh
    ${CADO_NFS_BINARY_DIR}/run_example.sh
    @ONLY)

install(PROGRAMS ${CADO_NFS_BINARY_DIR}/run_example.sh DESTINATION bin)
