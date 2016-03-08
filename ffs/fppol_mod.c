#include "fppol_mod.h"
#include "cppmeta.h"
#include "macros.h"



// Modular left shift by one (i.e., multiplication by t).
// /!\ Assume that m is monic and that p is reduced modulo m.
#define __DEF_FPPOLxx_SHL1MOD(sz)                                      \
  void fppol##sz##_multmod(fppol##sz##_ptr    r, fppol##sz##_srcptr p, \
                           fppol##sz##_srcptr m)                       \
  {                                                                    \
    int degp = fppol##sz##_deg(p);                                     \
    ASSERT(fppol##sz##_deg(m) > degp);                                 \
    if (degp+1 < fppol##sz##_deg(m)) {                                 \
      fppol##sz##_mul_ti(r, p, 1);                                     \
      return;                                                          \
    }                                                                  \
                                                                       \
    fp_t          lc;                                                  \
    fppol##sz##_t t;                                                   \
    IF(sz, EMPTY, fppol_init(t);, )                                    \
    fppol##sz##_get_coeff(lc, p, degp);                                \
    fppol##sz##_mul_ti(t, p, 1);                                       \
    fppol##sz##_smul  (r, m, lc);                                      \
    fppol##sz##_sub   (r, t, r);                                       \
    IF(sz, EMPTY, fppol_clear(t);, )                                   \
  }


// Modular multiplication.
// /!\ Assume that m is monic and that p is reduced modulo m.
#define __DEF_FPPOLxx_MULMOD(sz)                                            \
  void fppol##sz##_mulmod(fppol##sz##_ptr    r, fppol##sz##_srcptr p,       \
                          fppol##sz##_srcptr q, fppol##sz##_srcptr m)       \
  {                                                                         \
    IF(sz, EMPTY,                                                           \
       fppol_mul(r, p, q);                                                  \
       fppol_rem(r, r, m);                                                  \
    ,                                                                       \
       CAT(CAT(fppol, __MUL_SIZE(sz, sz, )), _t)              rr;           \
       CAT(CAT(fppol, __MUL_SIZE(sz, sz, )), _t)              mm;           \
       IF(__MUL_SIZE(sz, sz, ), EMPTY, fppol_inits(rr, mm, NULL);, )        \
       CAT(CAT(fppol, __MUL_SIZE(sz, sz, )), _mul_##sz##x##sz)(rr, p, q);   \
       CAT(CAT(fppol, __MUL_SIZE(sz, sz, )), _set_##sz)       (mm, m);      \
       CAT(CAT(fppol, __MUL_SIZE(sz, sz, )), _rem)            (rr, rr, mm); \
       CAT(fppol##sz##_set_, __MUL_SIZE(sz, sz, mp))          (r, rr);      \
       IF(__MUL_SIZE(sz, sz, ), EMPTY, fppol_clears(rr, mm, NULL);, )       \
    )                                                                       \
  }


// Modular inverse.
// Return a boolean, telling whether the inverse exists or not.
// /!\ Assume that m is monic and that p is reduced modulo m.
// TODO: if need be, in the case where the inverse does not exist, r
// could be set to a strict factor of q (or 0 if p = 0).
#define __DEF_FPPOLxx_INVMOD(sz)                                     \
  int fppol##sz##_invmod(fppol##sz##_ptr    r, fppol##sz##_srcptr p, \
                         fppol##sz##_srcptr m)                       \
  {                                                                  \
    /* We will see later if these requirements are a burden to    */ \
    /* the caller or not.                                         */ \
    ASSERT_ALWAYS(fppol##sz##_is_monic(m));                          \
    ASSERT_ALWAYS(fppol##sz##_deg(p) < fppol##sz##_deg(m));          \
                                                                     \
    /* Let's be deterministic and fix 1/0 = 0. */                    \
    if (UNLIKELY(fppol##sz##_is_zero(p))) {                          \
      fppol##sz##_set_zero(r);                                       \
      return 0;                                                      \
    }                                                                \
                                                                     \
    if (UNLIKELY(fppol##sz##_in_fp(p))) {                            \
      fppol##sz##_sinv(r, p);                                        \
      return 1;                                                      \
    }                                                                \
                                                                     \
    fppol##sz##_t r0, r1, v0, v1, t;                                 \
    IF(sz, EMPTY, fppol_inits(r0, r1, v0, v1, t, NULL);, )           \
                                                                     \
    /* Euclid-reduce (r0, r1) and maintain the following:         */ \
    /*   u0*m + v0*p = r0   (u0, u1 not computed!)                */ \
    /*   u1*m + v1*p = r1                                         */ \
    fppol##sz##_set(r0, m);                                          \
    fppol##sz##_set(r1, p);                                          \
    fppol##sz##_set_zero(v0);                                        \
    fppol##sz##_set_one (v1);                                        \
    int d0 = fppol##sz##_deg(r0);                                    \
    int d1 = fppol##sz##_deg(r1);                                    \
    IF(CMP(FP_SIZE, 2), EQ, ,                                        \
       fp_t lc0; fp_t lc1; fp_t lc;                                  \
       fppol##sz##_get_coeff(lc0, r0, d0);                           \
       fppol##sz##_get_coeff(lc1, r1, d1);                           \
    )                                                                \
    for (int d = d0 - d1; ; ) {                                      \
      fppol##sz##_mul_ti(t, r1, d);                                  \
      IF(CMP(FP_SIZE, 2), EQ, ,                                      \
         fp_div(lc, lc0, lc1);                                       \
         fppol##sz##_smul(t, t, lc);                                 \
      )                                                              \
      fppol##sz##_sub(r0, r0, t);                                    \
      fppol##sz##_mul_ti(t, v1, d);                                  \
      IF(CMP(FP_SIZE, 2), EQ, ,                                      \
         fppol##sz##_smul(t, t, lc);                                 \
      )                                                              \
      fppol##sz##_sub(v0, v0, t);                                    \
      d0 = fppol##sz##_deg(r0);                                      \
      IF(CMP(FP_SIZE, 2), EQ, ,                                      \
         if (d0 >= 0) fppol##sz##_get_coeff(lc0, r0, d0);            \
      )                                                              \
                                                                     \
      if (d0 <= 0) {                                                 \
        IF(CMP(FP_SIZE, 2), EQ,                                      \
           fppol##sz##_set (r, v0);,                                 \
           fppol##sz##_sdiv(r, v0, lc0);                             \
        )                                                            \
        IF(sz, EMPTY, fppol_clears(r0, r1, v0, v1, t, NULL);, )      \
        return !d0;                                                  \
      }                                                              \
                                                                     \
      d = d0 - d1;                                                   \
      if (d <= 0) {                                                  \
        fppol##sz##_swap(r0, r1);                                    \
        fppol##sz##_swap(v0, v1);                                    \
        IF(CMP(FP_SIZE, 2), EQ, , fp_swap(lc0, lc1);)                \
        int dt = d0; d0 = d1; d1 = dt;                               \
        d = -d;                                                      \
      }                                                              \
    }                                                                \
    ASSERT_ALWAYS(0); /* Should never go there. */                   \
  }


// All declarations bundled up into a single macro.
#define __DEF_FPPOLxx_MOD_ALL(sz) \
        __DEF_FPPOLxx_SHL1MOD(sz) \
        __DEF_FPPOLxx_MULMOD (sz) \
        __DEF_FPPOLxx_INVMOD (sz)

__DEF_FPPOLxx_MOD_ALL( 8)
__DEF_FPPOLxx_MOD_ALL(16)
__DEF_FPPOLxx_MOD_ALL(32)
__DEF_FPPOLxx_MOD_ALL(64)
__DEF_FPPOLxx_MOD_ALL()

#undef __DEF_FPPOLxx_SHL1MOD
#undef __DEF_FPPOLxx_MULMOD
#undef __DEF_FPPOLxx_INVMOD
#undef __DEF_FPPOLxx_MOD_ALL
