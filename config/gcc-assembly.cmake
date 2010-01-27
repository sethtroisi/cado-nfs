
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

