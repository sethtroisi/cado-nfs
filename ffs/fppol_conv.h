// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"
//#include "fppol_set.h"  // for fppolxx_deg, fppolxx_is_valid

#ifndef __FPPOL_CONV_H__
#define __FPPOL_CONV_H__

#include <stdint.h>

#include "cppmeta.h"

// Conversion of an polynomial to an uint64_t (evaluation in 2^__FP_BITS)
// Generic prototype:
//   uint64_t fppol<sz>_get_ui_sparse(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_GET_UI_SPARSE(sz)                                      \
  static inline                                                              \
  uint64_t fppol##sz##_get_ui_sparse(fppol##sz##_srcptr p)                   \
  { uint64_t r; __FP_GET_UI_SPARSE(sz, r, p); return r; }

// Conversion of an uint64_t (evaluation in 2^__FP_BITS) in polynomial
// Generic prototype:
//   uint64_t fppol<sz>_set_ui_sparse(fppol<sz>_ptr r, uint64_t x);
#define __DEF_FPPOLxx_SET_UI_SPARSE(sz)                                      \
  static inline                                                              \
  uint64_t fppol##sz##_set_ui_sparse(fppol##sz##_ptr r, uint64_t x)          \
  { __FP_SET_UI_SPARSE(sz, r, x); return fppol##sz##_is_valid(r); }


#define __DEF_FPPOLxx_CONV_ALL(sz)      \
        __DEF_FPPOLxx_GET_UI_SPARSE(sz) \
        __DEF_FPPOLxx_SET_UI_SPARSE(sz)

__DEF_FPPOLxx_CONV_ALL( 8)
__DEF_FPPOLxx_CONV_ALL(16)
__DEF_FPPOLxx_CONV_ALL(32)
__DEF_FPPOLxx_CONV_ALL(64)

#undef __DEF_FPPOLxx_GET_UI_SPARSE
#undef __DEF_FPPOLxx_SET_UI_SPARSE
#undef __DEF_FPPOLxx_CONV_ALL

#endif  /* __FPPOL_CONV_H__ */
