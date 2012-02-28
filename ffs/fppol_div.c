#include "fppol_div.h"



// Division with remainder.
// Compute (q, r) such that a = b*q + r, with deg(r) < deg(b).
// Return 1 in the case of a successful computation (b != 0), 0 otherwise.
#define __DEF_FPPOLxx_DIVREM(sz)                                           \
  int fppol##sz##_divrem(fppol##sz##_ptr    q, fppol##sz##_ptr    r,       \
                         fppol##sz##_srcptr a, fppol##sz##_srcptr b)       \
  {                                                                        \
    fp_t ilcb; /* Inverse of the leading coeff of b. */                    \
    int  dega, degb;                                                       \
                                                                           \
    /* Trivial cases. */                                                   \
    degb = fppol##sz##_deg(b);                                             \
    if (UNLIKELY(degb == -1))                                              \
      return 0;                                                            \
    if (UNLIKELY(degb == 0)) {                                             \
      fppol##sz##_get_coeff(ilcb, b, 0);                                   \
      fppol##sz##_sdiv(q, a, ilcb);                                        \
      fppol##sz##_set_zero(r);                                             \
      return 1;                                                            \
    }                                                                      \
    dega = fppol##sz##_deg(a);                                             \
    if (dega < degb) {                                                     \
      fppol##sz##_set_zero(q);                                             \
      fppol##sz##_set(r, a);                                               \
      return 1;                                                            \
    }                                                                      \
                                                                           \
    /* Check if input alias, and if so, work on a tmp space. */            \
    fppol##sz##_t   tq, tr;                                                \
    fppol##sz##_ptr qq, rr;                                                \
    if (q == a || q == b) {                                                \
      IF(sz, EMPTY, fppol_init(tq);, )                                     \
      qq = &tq[0];                                                         \
    }                                                                      \
    else                                                                   \
      qq = &q[0];                                                          \
    if (r == a || r == b) {                                                \
      IF(sz, EMPTY, fppol_init(tr);, )                                     \
      rr = &tr[0];                                                         \
    }                                                                      \
    else                                                                   \
      rr = &r[0];                                                          \
                                                                           \
    /* General case. */                                                    \
    fp_t          lcr, cq;                                                 \
    fppol##sz##_t tmp;                                                     \
    IF(sz, EMPTY, fppol_init(tmp);, )                                      \
                                                                           \
    /* Invariant: a = b*qq + rr.                                        */ \
    /* TODO: it is not clear that computing the degree at each step is  */ \
    /* interesting. It might be that doing the reduction for each coeff */ \
    /* of rr, be it zero or not, is faster.                             */ \
    fppol##sz##_set(rr, a);                                                \
    fppol##sz##_set_zero(qq);                                              \
    fppol##sz##_get_coeff(ilcb, b, degb);                                  \
    fp_inv(ilcb, ilcb);                                                    \
    for (int degr = dega; degr >= degb; ) {                                \
      fppol##sz##_get_coeff(lcr, rr, degr);                                \
      fp_mul(cq, lcr, ilcb);                                               \
      fppol##sz##_set_coeff(qq, cq, degr-degb);                            \
      fppol##sz##_smul(tmp, b, cq);                                        \
      fppol##sz##_shl(tmp, tmp, degr-degb);                                \
      fppol##sz##_sub(rr, rr, tmp);                                        \
      degr = fppol##sz##_deg(rr);                                          \
    }                                                                      \
                                                                           \
    /* If input were aliases, do the final copy. */                        \
    if (q == a || q == b) {                                                \
      fppol##sz##_set(q, qq);                                              \
      IF(sz, EMPTY, fppol_clear(tq);, )                                    \
    }                                                                      \
    if (r == a || r == b) {                                                \
      fppol##sz##_set(r, rr);                                              \
      IF(sz, EMPTY, fppol_clear(tr);, )                                    \
    }                                                                      \
    return 1;                                                              \
  }


// Quotient only.
// Same as above, but computes only q.
#define __DEF_FPPOLxx_DIV(sz)                                     \
  int fppol##sz##_div(fppol##sz##_ptr    q,                       \
                      fppol##sz##_srcptr a, fppol##sz##_srcptr b) \
  {                                                               \
    fppol##sz##_t r;                                              \
    IF(sz, EMPTY, fppol_init(r);, )                               \
    int rc = fppol##sz##_divrem(q, r, a, b);                      \
    IF(sz, EMPTY, fppol_clear(r);, )                              \
    return rc;                                                    \
  }


// Remainder only.
// Same as above, but computes only r.
// TODO: calling divrem for that one is suboptimal.
// See fppol_mod.c for an example in characteristic 2.
#define __DEF_FPPOLxx_REM(sz)                                     \
  int fppol##sz##_rem(fppol##sz##_ptr    r,                       \
                      fppol##sz##_srcptr a, fppol##sz##_srcptr b) \
  {                                                               \
    fppol##sz##_t q;                                              \
    IF(sz, EMPTY, fppol_init(q);, )                               \
    int rc = fppol##sz##_divrem(q, r, a, b);                      \
    IF(sz, EMPTY, fppol_clear(q);, )                              \
    return rc;                                                    \
  }


// All declarations bundled up into a single macro.
#define __DEF_FPPOLxx_DIV_ALL(sz) \
        __DEF_FPPOLxx_DIVREM (sz) \
        __DEF_FPPOLxx_DIV    (sz) \
        __DEF_FPPOLxx_REM    (sz)

__DEF_FPPOLxx_DIV_ALL(16)
__DEF_FPPOLxx_DIV_ALL(32)
__DEF_FPPOLxx_DIV_ALL(64)
__DEF_FPPOLxx_DIV_ALL()

#undef __DEF_FPPOLxx_DIVREM
#undef __DEF_FPPOLxx_DIV
#undef __DEF_FPPOLxx_REM
#undef __DEF_FPPOLxx_DIV_ALL
