// Always include "fppol.h" first so that everything gets loaded in the
// correct order.
#include "fppol.h"

#ifndef __FPPOL_GCD_H__
#define __FPPOL_GCD_H__

// GCD if p and q. 
// Warning: the result is not guaranteed to be monic.
#define __DECL_FPPOLxx_GCD(sz)                                      \
  void fppol##sz##_gcd(fppol##sz##_ptr    r, fppol##sz##_srcptr p,  \
                       fppol##sz##_srcptr q);

__DECL_FPPOLxx_GCD(16)
__DECL_FPPOLxx_GCD(32)
__DECL_FPPOLxx_GCD(64)
__DECL_FPPOLxx_GCD()

#undef __DECL_FPPOLxx_GCD

#endif   /* __FPPOL_GCD_H__ */
