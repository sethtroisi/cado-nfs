
# gcc-style x86_64 inline assembly.
message(STATUS "Trying to compile gcc-style x86_64 inline assembly")
try_compile(HAVE_GCC_STYLE_AMD64_INLINE_ASM
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/inline-assembly.c)
if(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
    message(STATUS "Trying to compile gcc-style x86_64 inline assembly -- Success")
else(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
    message(STATUS "Trying to compile gcc-style x86_64 inline assembly -- Failed")
endif(HAVE_GCC_STYLE_AMD64_INLINE_ASM)

