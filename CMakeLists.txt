
project(CADO_NFS)

# May be overridden by environment.
set(CADO_VERSION "0.1-rc2")

#####################################################################
# cmake version -- set_property requires 2.6. Note that haters of progress
# report things (this include me, casually) have to wait for cmake 2.7, since
# it has only recently become possible to disable them:
# http://public.kitware.com/Bug/view.php?id=8726

cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
        cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#####################################################################
# What are we building ?

execute_process(COMMAND scripts/version.sh
        WORKING_DIRECTORY ${CADO_NFS_SOURCE_DIR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE CADO_REV)
message(STATUS "Configuring CADO revision ${CADO_REV}")

#####################################################################
# Get some overrides from the environment.

macro(get_override extname intname)
string(COMPARE NOTEQUAL "$ENV{${extname}}" "" HAS_${intname}_OVERRIDE)
if(HAS_${intname}_OVERRIDE)
    set(${intname} "$ENV{${extname}}")
    message(STATUS "${extname}=${${intname}}")
endif(HAS_${intname}_OVERRIDE)
endmacro(get_override varname)

# Syntax is (user-friendly name) (cmake name)
# Make sure everything is exported from scripts/call_cmake.sh
get_override(PREFIX               CMAKE_INSTALL_PREFIX)
get_override(CFLAGS               CMAKE_C_FLAGS)
get_override(CXXFLAGS             CMAKE_CXX_FLAGS)
get_override(CC                   CMAKE_C_COMPILER)
get_override(CXX                  CMAKE_CXX_COMPILER)
get_override(CADO_VERSION         CADO_VERSION)
get_override(MPI                  MPI)
get_override(PTHREADS             PTHREADS)
set(GF2X_CONFIGURE_FLAGS --silent --disable-shared --disable-dependency-tracking)
get_override(GF2X_CONFIGURE_FLAGS GF2X_CONFIGURE_FLAGS)

#####################################################################
# Add the CFLAGS that are needed in all cases.

if(CMAKE_COMPILER_IS_GNUCC)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g -W -Wall")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -W -Wall")
endif(CMAKE_COMPILER_IS_GNUCC)

#####################################################################
# Set dependencies for external software.

# You can force a path to gmp.h using the environment variables GMP, or
# GMP_INCDIR and GMP_LIBDIR
string(COMPARE NOTEQUAL "$ENV{GMP}" "" HAS_GMP_OVERRIDE)
if (HAS_GMP_OVERRIDE)
    message(STATUS "Adding $ENV{GMP} to the search path for Gnu MP")
    set(GMP_INCDIR_HINTS "$ENV{GMP}/include" ${GMP_INCDIR_HINTS})
    set(GMP_INCDIR_HINTS "$ENV{GMP}"         ${GMP_INCDIR_HINTS})
    set(GMP_LIBDIR_HINTS "$ENV{GMP}/lib"     ${GMP_LIBDIR_HINTS})
    set(GMP_LIBDIR_HINTS "$ENV{GMP}/.libs"   ${GMP_LIBDIR_HINTS})
endif(HAS_GMP_OVERRIDE)
string(COMPARE NOTEQUAL "$ENV{GMP_INCDIR}" "" HAS_GMP_INCDIR_OVERRIDE)
if (HAS_GMP_INCDIR_OVERRIDE)
    message(STATUS "Adding $ENV{GMP_INCDIR} to the search path for Gnu MP")
    set(GMP_INCDIR_HINTS "$ENV{GMP_INCDIR}" ${GMP_INCDIR_HINTS})
endif(HAS_GMP_INCDIR_OVERRIDE)
string(COMPARE NOTEQUAL "$ENV{GMP_LIBDIR}" "" HAS_GMP_LIBDIR_OVERRIDE)
if (HAS_GMP_LIBDIR_OVERRIDE)
    message(STATUS "Adding $ENV{GMP_LIBDIR} to the search path for Gnu MP")
    set(GMP_LIBDIR_HINTS "$ENV{GMP_LIBDIR}"     ${GMP_LIBDIR_HINTS})
endif(HAS_GMP_LIBDIR_OVERRIDE)
find_path   (GMP_INCDIR gmp.h HINTS ${GMP_INCDIR_HINTS} DOC "Gnu MP headers")
find_library(GMP_LIB    gmp   HINTS ${GMP_LIBDIR_HINTS} DOC "Gnu MP library")
# Yeah. CMake docs defines the ``PATH'' to a file as being its dirname. Very
# helpful documentation there :-((
get_filename_component(GMP_LIBDIR ${GMP_LIB} PATH)
message(STATUS "GMP_INCDIR=${GMP_INCDIR}")
message(STATUS "GMP_LIBDIR=${GMP_LIBDIR}")
string(COMPARE NOTEQUAL "${GMP_INCDIR}" GMP_INCDIR-NOTFOUND GMP_INCDIR_OK)
string(COMPARE NOTEQUAL "${GMP_LIBDIR}" GMP_LIBDIR-NOTFOUND GMP_LIBDIR_OK)
if(GMP_INCDIR_OK)
include_directories(${GMP_INCDIR})
endif(GMP_INCDIR_OK)
if(GMP_LIBDIR_OK)
link_directories(${GMP_LIBDIR})
endif(GMP_LIBDIR_OK)

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
    COMMAND cd gf2x && ${CADO_NFS_SOURCE_DIR}/gf2x/configure ${GF2X_CONFIGURE_FLAGS}
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
    ESCAPE_QUOTES
    @ONLY)

configure_file(
    ${CADO_NFS_SOURCE_DIR}/scripts/fixup_params_file.sh.in
    ${CADO_NFS_BINARY_DIR}/fixup_params_file.sh
    ESCAPE_QUOTES
    @ONLY)

install(PROGRAMS ${CADO_NFS_BINARY_DIR}/run_example.sh DESTINATION bin)
install(PROGRAMS ${CADO_NFS_BINARY_DIR}/fixup_params_file.sh DESTINATION bin)
