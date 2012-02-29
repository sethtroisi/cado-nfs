#include "fppol_internal.h"
#include "cppmeta.h"
#include "macros.h"
#include "fppol_gcd.h"

/*
 * Compute the GCD with an additive Euclidian algorithm.
 * The result is *not* made monic.
 */

#define __DEF_FPPOLxx_GCD(sz)                                        \
  void fppol##sz##_gcd(fppol##sz##_ptr    r, fppol##sz##_srcptr p,   \
                       fppol##sz##_srcptr q)                         \
  {                                                                  \
    if (UNLIKELY(fppol##sz##_is_zero(p))) {                          \
      fppol##sz##_set(r, q);                                         \
      return;                                                        \
    }                                                                \
    if (UNLIKELY(fppol##sz##_is_zero(q))) {                          \
      fppol##sz##_set(r, p);                                         \
      return;                                                        \
    }                                                                \
    int degp = fppol##sz##_deg(p);                                   \
    int degq = fppol##sz##_deg(q);                                   \
    if (UNLIKELY(degp == 0 || degq == 0)) {                          \
      fppol##sz##_set_ti(r, 1);                                      \
      return;                                                        \
    }                                                                \
    fppol##sz##_t r0, r1, t;                                         \
    IF(sz, EMPTY, fppol_inits(r0, r1, t, NULL);, )                   \
    int d0, d1;                                                      \
    if (degp >= degq) {                                              \
        fppol##sz##_set(r0, p);                                      \
        fppol##sz##_set(r1, q);                                      \
        d0 = degp;                                                   \
        d1 = degq;                                                   \
    } else {                                                         \
        fppol##sz##_set(r0, q);                                      \
        fppol##sz##_set(r1, p);                                      \
        d0 = degq;                                                   \
        d1 = degp;                                                   \
    }                                                                \
    IF(CMP(FP_SIZE, 2), EQ, ,                                        \
       fp_t lc0; fp_t lc1; fp_t lc;                                  \
       fppol##sz##_get_coeff(lc0, r0, d0);                           \
       fppol##sz##_get_coeff(lc1, r1, d1);                           \
    )                                                                \
    for (int d = d0 - d1; ; ) {                                      \
      fppol##sz##_shl(t, r1, d);                                     \
      IF(CMP(FP_SIZE, 2), EQ, ,                                      \
         fp_div(lc, lc0, lc1);                                       \
         fppol##sz##_smul(t, t, lc);                                 \
      )                                                              \
      fppol##sz##_sub(r0, r0, t);                                    \
      d0 = fppol##sz##_deg(r0);                                      \
      IF(CMP(FP_SIZE, 2), EQ, ,                                      \
         if (d0 >= 0) fppol##sz##_get_coeff(lc0, r0, d0);            \
      )                                                              \
                                                                     \
      if (d0 <= 0) {                                                 \
          if (d0 < 0)                                                \
            fppol##sz##_set (r, r1);                                 \
          else                                                       \
            fppol##sz##_set_ti (r, 0);                               \
          IF(sz, EMPTY, fppol_clears(r0, r1, t, NULL);, )            \
          return ;                                                   \
      }                                                              \
                                                                     \
      d = d0 - d1;                                                   \
      if (d <= 0) {                                                  \
        __fppol##sz##_swap(r0, r1);                                  \
        IF(CMP(FP_SIZE, 2), EQ, , __fp_swap(lc0, lc1);)              \
        int dt = d0; d0 = d1; d1 = dt;                               \
        d = -d;                                                      \
      }                                                              \
    }                                                                \
    ASSERT_ALWAYS(0); /* Should never go there. */                   \
  }


__DEF_FPPOLxx_GCD(16)
__DEF_FPPOLxx_GCD(32)
__DEF_FPPOLxx_GCD(64)
__DEF_FPPOLxx_GCD()

#undef __DEF_FPPOLxx_GCD
