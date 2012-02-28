// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_IO_H__
#define __FPPOL_IO_H__

#include <stdio.h>



/* Fixed-size polynomials.
 *****************************************************************************/

// Representation of polynomials is in hexadecimal.

// Conversion to string.
// Generic prototype:
//   char *fppol<sz>_get_str(char *str, fppol<sz>_srcptr p);
#define __DECL_FPPOLxx_GET_STR(sz)                            \
  char *fppol##sz##_get_str(char *str, fppol##sz##_srcptr p);


// Conversion from string.
// Return 1 if successful.
// Generic prototype:
//   int fppol<sz>_set_str(fppol<sz>_ptr r, const char *str);
#define __DECL_FPPOLxx_SET_STR(sz)                             \
  int fppol##sz##_set_str(fppol##sz##_ptr r, const char *str);


// Output to stream.
// Generic prototype:
//   void fppol<sz>_out(FILE *f, fppol<sz>_srcptr p);
#define __DECL_FPPOLxx_OUT(sz)                         \
  void fppol##sz##_out(FILE *f, fppol##sz##_srcptr p);


// Input from stream.
// Return 1 if successful.
// Generic prototype:
//   int fppol<sz>_inp(fppol<sz>_ptr r, FILE *f);
#define __DECL_FPPOLxx_INP(sz)                     \
  int fppol##sz##_inp(fppol##sz##_ptr r, FILE *f);


// All declarations bundled up into a single macro.
#define __DECL_FPPOLxx_IO_ALL(sz)  \
        __DECL_FPPOLxx_GET_STR(sz) \
        __DECL_FPPOLxx_SET_STR(sz) \
        __DECL_FPPOLxx_OUT    (sz) \
        __DECL_FPPOLxx_INP    (sz)

__DECL_FPPOLxx_IO_ALL(16)
__DECL_FPPOLxx_IO_ALL(32)
__DECL_FPPOLxx_IO_ALL(64)

#undef __DECL_FPPOLxx_GET_STR
#undef __DECL_FPPOLxx_SET_STR
#undef __DECL_FPPOLxx_OUT
#undef __DECL_FPPOLxx_INP
#undef __DECL_FPPOLxx_IO_ALL



/* Multiprecision polynomials.
 *****************************************************************************/

// Conversion to string.
char *fppol_get_str(char *str, fppol_srcptr p);

// Conversion from string.
// Return 1 if successful.
int fppol_set_str(fppol_ptr r, const char *str);

// Output to stream.
void fppol_out(FILE *f, fppol_srcptr p);

// Input from stream.
// Return 1 if successful.
int fppol_inp(fppol_ptr r, FILE *f);

#endif  /* __FPPOL_IO_H__ */
