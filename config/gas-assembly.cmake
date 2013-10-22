# gas-syntax assembly sources.
message(STATUS "Trying to compile gas-syntax assembly sources")
try_compile(HAVE_GAS_SYNTAX_ASSEMBLY_SOURCES
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/asm_sample.c
    COMPILE_DEFINITIONS "-x assembler")
if(HAVE_GAS_SYNTAX_ASSEMBLY_SOURCES)
    message(STATUS "Trying to compile gas-syntax assembly sources -- Success")
else(HAVE_GAS_SYNTAX_ASSEMBLY_SOURCES)
    message(STATUS "Trying to compile gas-syntax assembly sources -- Failed")
endif(HAVE_GAS_SYNTAX_ASSEMBLY_SOURCES)


