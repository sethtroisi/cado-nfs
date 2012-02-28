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
    /* ASSERT_ALWAYS(fppol##sz##_is_monic(m)); */                       \
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
        __DEF_FPPOLxx_INVMOD (sz)

__DEF_FPPOLxx_MOD_ALL(16)
__DEF_FPPOLxx_MOD_ALL(32)
__DEF_FPPOLxx_MOD_ALL(64)
__DEF_FPPOLxx_MOD_ALL()

#undef __DEF_FPPOLxx_SHL1MOD
#undef __DEF_FPPOLxx_MULMOD
#undef __DEF_FPPOLxx_INVMOD
#undef __DEF_FPPOLxx_MOD_ALL


// r = p mod q
// Specific to charac 2.
void fppol32_mod_64(fppol32_ptr r, fppol64_srcptr p, fppol32_srcptr q)
{
    ASSERT_ALWAYS(FP_CHAR == 2);
    int degq = fppol32_deg(q);
    uint64_t mask = ((uint64_t)1U)<<63;
    uint64_t pp = p[0];
    uint64_t qq = ((uint64_t)q[0]) << (63-degq);
    for (int i = 63; i >= degq; --i) {
        pp ^= (pp & mask) ? qq : 0;
        mask >>= 1;
        qq >>= 1;
    }
    r[0] = pp;
}

// r = p*q mod m
void fppol32_mulmod(fppol32_ptr r, fppol32_srcptr p, fppol32_srcptr q,
        fppol32_srcptr m)
{
    fppol64_t rr;
    fppol64_mul_32x32(rr, p, q);
    fppol32_mod_64(r, rr, m);
}


void fppol64_mulmod(fppol64_ptr r, fppol64_srcptr p, fppol64_srcptr q,
                fppol64_srcptr m)
{
    fppol_t pp, mm;
    fppol_init(pp);
    fppol_init(mm);

    fppol_set_64(pp, p);
    fppol_set_64(mm, m);
    fppol_mul_64(pp, pp, q);
    fppol_rem(pp, pp, mm);
    fppol64_set_mp(r, pp);
    
    fppol_clear(pp);
    fppol_clear(mm);
}
