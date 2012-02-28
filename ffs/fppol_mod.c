#include "fppol_mod.h"
#include "macros.h"



// Modular left shift by one (i.e., multiplication by t).
// /!\ Assume that m is monic and that p is reduced modulo m.
#define __DEF_FPPOLxx_SHL1MOD(sz)                                      \
  void fppol##sz##_shl1mod(fppol##sz##_ptr    r, fppol##sz##_srcptr p, \
                           fppol##sz##_srcptr m)                       \
  {                                                                    \
    int degp = fppol##sz##_deg(p);                                     \
    ASSERT(fppol##sz##_deg(m) > degp);                                 \
    if (degp+1 < fppol##sz##_deg(m)) {                                 \
      fppol##sz##_shl(r, p, 1);                                        \
      return;                                                          \
    }                                                                  \
                                                                       \
    fp_t          lc;                                                  \
    fppol##sz##_t t;                                                   \
    IF(sz, EMPTY, fppol_init(t);, )                                    \
    fppol##sz##_get_coeff(lc, p, degp);                                \
    fppol##sz##_shl (t, p, 1);                                         \
    fppol##sz##_smul(r, m, lc);                                        \
    fppol##sz##_sub (r, t, r);                                         \
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
       CAT(CAT(fppol, __MUL_SIZE(sz, sz, )), _mul_##sz##x##sz)(rr, p, q);   \
       CAT(CAT(fppol, __MUL_SIZE(sz, sz, )), _set_##sz)       (mm, m);      \
       CAT(CAT(fppol, __MUL_SIZE(sz, sz, )), _rem)            (rr, rr, mm); \
       CAT(fppol##sz##_set_, __MUL_SIZE(sz, sz, mp))          (r, rr);      \
    )                                                                       \
  }


// Modular inverse.
// Return a boolean, telling whether the inverse exists or not.
// /!\ Assume that m is monic and that p is reduced modulo m.
// TODO: if need be, in the case where the inverse does not exist, r
// could be set to a strict factor of q (or 0 if p = 0).
#define __DEF_FPPOLxx_INVMOD(sz)                                        \
  int fppol##sz##_invmod(fppol##sz##_ptr    r, fppol##sz##_srcptr p,    \
                         fppol##sz##_srcptr m)                          \
  {                                                                     \
    /* We will see later if these requirements are a burden to       */ \
    /* the caller or not.                                            */ \
    ASSERT_ALWAYS(fppol##sz##_is_monic(m));                             \
    ASSERT_ALWAYS(fppol##sz##_deg(p) < fppol##sz##_deg(m));             \
    /* Right now, we never bother with non-monic polynomials         */ \
    /* Let's assert that we are in characteristic 2.                 */ \
    /* TODO: implement for characteristic 3.                         */ \
    ASSERT_ALWAYS(FP_CHAR == 2);                                        \
                                                                        \
    /* Let's be deterministic and fix 1/0 = 0. */                       \
    if (UNLIKELY(fppol##sz##_is_zero(p))) {                             \
      fppol##sz##_set_zero(r);                                          \
      return 0;                                                         \
    }                                                                   \
                                                                        \
    fppol##sz##_t a, b, c, v, vp, nv;                                   \
    IF(sz, EMPTY, fppol_inits(a, b, c, v, vp, nv, NULL);, )             \
                                                                        \
    /* Euclide-reduce (a,b) and maintain the following:              */ \
    /*   u *m +  v*p = a    <- not computed!                         */ \
    /*   up*m + vp*p = b                                             */ \
    fppol##sz##_set(a, m);                                              \
    fppol##sz##_set(b, p);                                              \
    fppol##sz##_set_zero(v);                                            \
    fppol##sz##_set_one(vp);                                            \
    for (int dega = fppol##sz##_deg(a), degb = fppol##sz##_deg(b); ;) { \
      for (int d; (d = dega-degb) >= 0; ) {                             \
        fppol##sz##_shl(c, b, d); /* TODO: make monic. */               \
        fppol##sz##_sub(c, a, c);                                       \
        fppol##sz##_set(a, b);                                          \
        fppol##sz##_set(b, c);                                          \
        dega = degb;                                                    \
        degb = fppol##sz##_deg(b);                                      \
                                                                        \
        fppol##sz##_shl(nv, vp, d);                                     \
        fppol##sz##_sub(nv, v,  nv);                                    \
        fppol##sz##_set(v,  vp);                                        \
        fppol##sz##_set(vp, nv);                                        \
                                                                        \
        if (degb <= 0) {                                                \
          fppol##sz##_set(r, vp);                                       \
          IF(sz, EMPTY, fppol_clears(a, b, c, v, vp, nv, NULL);, )      \
          return !degb;                                                 \
        }                                                               \
      }                                                                 \
                                                                        \
      for (int d; (d = dega-degb) <= 0; ) {                             \
        fppol##sz##_shl(c, a, -d); /* TODO: make monic. */              \
        fppol##sz##_sub(b, b,  c);                                      \
        degb = fppol##sz##_deg(b);                                      \
                                                                        \
        fppol##sz##_shl(nv, v,  -d);                                    \
        fppol##sz##_sub(vp, vp, nv);                                    \
                                                                        \
        if (degb <= 0) {                                                \
          fppol##sz##_set(r, vp);                                       \
          IF(sz, EMPTY, fppol_clears(a, b, c, v, vp, nv, NULL);, )      \
          return !degb;                                                 \
        }                                                               \
      }                                                                 \
    }                                                                   \
    ASSERT_ALWAYS(0); /* Should never go there. */                      \
  }


// All declarations bundled up into a single macro.
#define __DEF_FPPOLxx_MOD_ALL(sz) \
        __DEF_FPPOLxx_SHL1MOD(sz) \
        __DEF_FPPOLxx_MULMOD (sz) \
        __DEF_FPPOLxx_INVMOD (sz)

__DEF_FPPOLxx_MOD_ALL(16)
__DEF_FPPOLxx_MOD_ALL(32)
__DEF_FPPOLxx_MOD_ALL(64)
__DEF_FPPOLxx_MOD_ALL()

#undef __DEF_FPPOLxx_SHL1MOD
#undef __DEF_FPPOLxx_MULMOD
#undef __DEF_FPPOLxx_INVMOD
#undef __DEF_FPPOLxx_MOD_ALL
