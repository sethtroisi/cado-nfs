
message(STATUS "Checking for alloca.h")
find_path   (ALLOCA_INCDIR alloca.h DOC "alloca.h")
string(COMPARE NOTEQUAL "${ALLOC_INCDIR}" ALLOCA_INCDIR-NOTFOUND HAVE_ALLOCA_H)
if (HAVE_ALLOCA_H)
message(STATUS "Checking for alloca.h -- found")
else (HAVE_ALLOCA_H)
message(STATUS "Checking for alloca.h -- not found")
endif (HAVE_ALLOCA_H)
