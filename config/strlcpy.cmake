# Test whether strlcpy() and strlcat() exist

INCLUDE (CheckIncludeFiles)
include(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(strlcpy	HAVE_STRLCPY)
CHECK_FUNCTION_EXISTS(strlcat	HAVE_STRLCAT)
