// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_TEST_H__
#define __FPPOL_TEST_H__

#include "macros.h"



/* Fixed-size polynomials.
 *****************************************************************************/

// Bitwise OR of all bit vectors: coefficient-wise non-zero test.
// Generic prototype:
//   uint<sz>_t fppol<sz>_fold_or(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_FOLD_OR(sz)                        \
  static inline                                          \
  uint##sz##_t fppol##sz##_fold_or(fppol##sz##_srcptr p) \
  { uint##sz##_t t = p[0];                               \
    for (unsigned k = 1; k < __FP_BITS; ++k) t |= p[k];  \
    return t; }


// Degree.
// By convention, deg(0) = -1.
// Generic prototype:
//   int fppol<sz>_deg(fppol<sz>_srcptr p);
#ifndef __GNUC__
#define __DEF_FPPOLxx_DEG(sz)                \
  static inline                              \
  int fppol##sz##_deg(fppol##sz##_srcptr p)  \
  {                                          \
    uint##sz##_t t = fppol##sz##_fold_or(p); \
    if (!t) return -1;                       \
    int d = 0;                               \
    for (unsigned k = sz; k >>= 1; )         \
      if (t >> k) d += k, t >>= k;           \
    return d;                                \
  }
#else
// GCC provides useful builtins that translate into the BSF asm instructions
// on Intel / AMD cpus.
#define __DEF_FPPOLxx_DEG(sz)                               \
  static inline                                             \
  int fppol##sz##_deg(fppol##sz##_srcptr p)                 \
  {                                                         \
    uint##sz##_t t = fppol##sz##_fold_or(p);                \
    if (!t) return -1;                                      \
    if (sizeof(uint##sz##_t) <= sizeof(unsigned int))       \
      return sz-1-__builtin_clz(t);                         \
    if (sizeof(uint##sz##_t) == sizeof(unsigned long))      \
      return sz-1-__builtin_clzl(t);                        \
    if (sizeof(uint##sz##_t) == sizeof(unsigned long long)) \
      return sz-1-__builtin_clzll(t);                       \
    int d = 0;                                              \
    for (unsigned k = sz; k >>= 1; )                        \
      if (t >> k) d += k, t >>= k;                          \
    return d;                                               \
  }
#endif

// Test if zero.
// Generic prototype:
//   int fppol<sz>_is_zero(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_IS_ZERO(sz)               \
  static inline                                 \
  int fppol##sz##_is_zero(fppol##sz##_srcptr p) \
  { return !fppol##sz##_fold_or(p); }


// Test if in GF(p) (i.e., deg <= 0).
// Generic prototype:
//   int fppol<sz>_in_fp(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_IN_FP(sz)               \
  static inline                               \
  int fppol##sz##_in_fp(fppol##sz##_srcptr p) \
  { return !(fppol##sz##_fold_or(p) >> 1); }


// Test if equal.
// Generic prototype:
//   int fppol<sz>_eq(fppol<sz>_srcptr p, fppol<sz>_srcptr q);
#define __DEF_FPPOLxx_EQ(sz)                                     \
  static inline                                                  \
  int fppol##sz##_eq(fppol##sz##_srcptr p, fppol##sz##_srcptr q) \
  { int rc = p[0] == q[0];                                       \
    for (unsigned k = 1; k < __FP_BITS; ++k)                     \
      rc = rc && p[k] == q[k];                                   \
    return rc; }


// Test if monic.
// Generic prototype:
//   int fppol<sz>_is_monic(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_IS_MONIC(sz)                           \
  /* Forward declaration of fppol<sz>_get_coeff. */          \
  static inline                                              \
  void fppol##sz##_get_coeff(fp_ptr r, fppol##sz##_srcptr p, \
                             unsigned i);                    \
  static inline                                              \
  int fppol##sz##_is_monic(fppol##sz##_srcptr p)             \
  { int d = fppol##sz##_deg(p);                              \
    if (UNLIKELY(d == -1)) return 0;                         \
    fp_t lc;                                                 \
    fppol##sz##_get_coeff(lc, p, d);                         \
    return fp_is_one(lc); }


// Test if valid representation.
// Generic prototype:
//   int fppol<sz>_is_valid(fppol<sz>_srcptr p);
#define __DEF_FPPOLxx_IS_VALID(sz)                            \
  static inline                                               \
  int fppol##sz##_is_valid(MAYBE_UNUSED fppol##sz##_srcptr p) \
  { return __FP_IS_VALID(sz, p); }


// All definitions bundled up into a single macro.
#define __DEF_FPPOLxx_TEST_ALL(sz) \
        __DEF_FPPOLxx_FOLD_OR (sz) \
        __DEF_FPPOLxx_DEG     (sz) \
        __DEF_FPPOLxx_IS_ZERO (sz) \
        __DEF_FPPOLxx_IN_FP   (sz) \
        __DEF_FPPOLxx_EQ      (sz) \
        __DEF_FPPOLxx_IS_MONIC(sz) \
        __DEF_FPPOLxx_IS_VALID(sz)

__DEF_FPPOLxx_TEST_ALL(16)
__DEF_FPPOLxx_TEST_ALL(32)
__DEF_FPPOLxx_TEST_ALL(64)

#undef __DEF_FPPOLxx_FOLD_OR
#undef __DEF_FPPOLxx_DEG
#undef __DEF_FPPOLxx_IS_ZERO
#undef __DEF_FPPOLxx_IN_FP
#undef __DEF_FPPOLxx_EQ
#undef __DEF_FPPOLxx_IS_MONIC
#undef __DEF_FPPOLxx_IS_VALID
#undef __DEF_FPPOLxx_TEST_ALL



/* Multiprecision polynomials.
 *****************************************************************************/

// Degree.
static inline
int fppol_deg(fppol_srcptr p)
{ return p->deg; }

// Test if zero.
static inline
int fppol_is_zero(fppol_srcptr p)
{ return p->deg == -1; }

// Test if in GF(p) (i.e., deg <= 0).
static inline
int fppol_in_fp(fppol_srcptr p)
{ return p->deg <= 0; }

// Test if equal.
int fppol_eq(fppol_srcptr p, fppol_srcptr q);

// Test if monic.
int fppol_is_monic(fppol_srcptr p);

// Test if valid representation.
int fppol_is_valid(fppol_srcptr p);

#endif  /* __FPPOL_TEST_H__ */

