// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_IO_H__
#define __FPPOL_IO_H__

#include <stdio.h>



// /!\ Representation of polynomials is in hexadecimal.

// Length of string representation, not counting the null-terminator.
// Generic prototype:
//   size_t fppol<sz>_strlen(fppol<sz>_srcptr p);
#define __DECL_FPPOLxx_STRLEN(sz)                 \
  static inline                                   \
  size_t fppol##sz##_strlen(fppol##sz##_srcptr p) \
  { int d = fppol##sz##_deg(p);                   \
    return d < 0 ? 1 : (__FP_BITS*(d+1)+3)>>2; }


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
        __DECL_FPPOLxx_STRLEN (sz) \
        __DECL_FPPOLxx_GET_STR(sz) \
        __DECL_FPPOLxx_SET_STR(sz) \
        __DECL_FPPOLxx_OUT    (sz) \
        __DECL_FPPOLxx_INP    (sz)

__DECL_FPPOLxx_IO_ALL( 8)
__DECL_FPPOLxx_IO_ALL(16)
__DECL_FPPOLxx_IO_ALL(32)
__DECL_FPPOLxx_IO_ALL(64)
__DECL_FPPOLxx_IO_ALL()

#undef __DECL_FPPOLxx_STRLEN
#undef __DECL_FPPOLxx_GET_STR
#undef __DECL_FPPOLxx_SET_STR
#undef __DECL_FPPOLxx_OUT
#undef __DECL_FPPOLxx_INP
#undef __DECL_FPPOLxx_IO_ALL

#endif  /* __FPPOL_IO_H__ */
