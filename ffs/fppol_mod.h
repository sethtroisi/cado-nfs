// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_MOD_H__
#define __FPPOL_MOD_H__



// Modular left shift by one (i.e., multiplication by t).
// /!\ Assume that m is monic and that p is reduced modulo m.
// Generic prototype:
//   void fppol<sz>_shl1mod(fppol<sz>_ptr    r, fppol<sz>_srcptr p,
//                          fppol<sz>_srcptr m);
#define __DECL_FPPOLxx_SHL1MOD(sz)                                     \
  void fppol##sz##_shl1mod(fppol##sz##_ptr    r, fppol##sz##_srcptr p, \
                           fppol##sz##_srcptr m);


// Modular multiplication.
// /!\ Assume that m is monic and that p is reduced modulo m.
// Generic prototype:
//   void fppol<sz>_mulmod(fppol<sz>_ptr    r, fppol<sz>_srcptr p,
//                         fppol<sz>_srcptr q, fppol<sz>_srcptr m);
#define __DECL_FPPOLxx_MULMOD(sz)                                      \
  void fppol##sz##_mulmod(fppol##sz##_ptr    r, fppol##sz##_srcptr p,  \
                          fppol##sz##_srcptr q, fppol##sz##_srcptr m);


// Modular inverse.
// Return a boolean, telling whether the inverse exists or not.
// /!\ Assume that m is monic and that p is reduced modulo m.
// Generic prototype:
//   int fppol<sz>_invmod(fppol<sz>_ptr    r, fppol<sz>_srcptr p,
//                        fppol<sz>_srcptr m);
#define __DECL_FPPOLxx_INVMOD(sz)                                    \
  int fppol##sz##_invmod(fppol##sz##_ptr    r, fppol##sz##_srcptr p, \
                         fppol##sz##_srcptr m);


// All declarations bundled up into a single macro.
#define __DECL_FPPOLxx_MOD_ALL(sz) \
        __DECL_FPPOLxx_SHL1MOD(sz) \
        __DECL_FPPOLxx_MULMOD (sz) \
        __DECL_FPPOLxx_INVMOD (sz)

__DECL_FPPOLxx_MOD_ALL(16)
__DECL_FPPOLxx_MOD_ALL(32)
__DECL_FPPOLxx_MOD_ALL(64)
__DECL_FPPOLxx_MOD_ALL()

#undef __DECL_FPPOLxx_SHL1MOD
#undef __DECL_FPPOLxx_MULMOD
#undef __DECL_FPPOLxx_INVMOD
#undef __DECL_FPPOLxx_MOD_ALL

#endif   /* __FPPOL_MOD_H__ */
