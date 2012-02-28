// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_DIV_H__
#define __FPPOL_DIV_H__



// Division with remainder.
// Compute (q, r) such that a = b*q + r, with deg(r) < deg(b).
// Return 1 in the case of a successful computation (b != 0), 0 otherwise.
// Generic prototype:
//   int fppol<sz>_divrem(fppol<sz>_ptr    q, fppol<sz>_ptr    r,
//                        fppol<sz>_srcptr a, fppol<sz>_srcptr b);
#define __DECL_FPPOLxx_DIVREM(sz)                                     \
  int fppol##sz##_divrem(fppol##sz##_ptr    q, fppol##sz##_ptr    r,  \
                         fppol##sz##_srcptr a, fppol##sz##_srcptr b);


// Quotient only.
// Same as above, but computes only q.
// Generic prototype:
//   int fppol<sz>_div(fppol<sz>_ptr    q,
//                     fppol<sz>_srcptr a, fppol<sz>_srcptr b);
#define __DECL_FPPOLxx_DIV(sz)                                     \
  int fppol##sz##_div(fppol##sz##_ptr    q,                        \
                      fppol##sz##_srcptr a, fppol##sz##_srcptr b);


// Remainder only.
// Same as above, but computes only r.
// Generic prototype:
//   int fppol<sz>_rem(fppol<sz>_ptr    r,
//                     fppol<sz>_srcptr a, fppol<sz>_srcptr b);
#define __DECL_FPPOLxx_REM(sz)                                     \
  int fppol##sz##_rem(fppol##sz##_ptr    r,                        \
                      fppol##sz##_srcptr a, fppol##sz##_srcptr b);


// All declarations bundled up into a single macro.
#define __DECL_FPPOLxx_DIV_ALL(sz) \
  __DECL_FPPOLxx_DIVREM       (sz) \
  __DECL_FPPOLxx_DIV          (sz) \
  __DECL_FPPOLxx_REM          (sz)

__DECL_FPPOLxx_DIV_ALL(16)
__DECL_FPPOLxx_DIV_ALL(32)
__DECL_FPPOLxx_DIV_ALL(64)
__DECL_FPPOLxx_DIV_ALL()

#undef __DECL_FPPOLxx_DIVREM
#undef __DECL_FPPOLxx_DIV
#undef __DECL_FPPOLxx_REM
#undef __DECL_FPPOLxx_DIV_ALL

#endif   /* __FPPOL_DIV_H__ */
