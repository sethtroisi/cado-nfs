#include <string.h>

#include "fppol_mul.h"
#include "fppol_internal.h"

#if defined(USE_F2) && defined(HAVE_GF2X)
#include <gf2x.h>
#endif


/* Multiprecision polynomials.
 *****************************************************************************/

// Multiplication.
void fppol_mul(fppol_ptr r, fppol_srcptr p, fppol_srcptr q)
{
  // Multiplication by zero.
  if (UNLIKELY(fppol_is_zero(p) || fppol_is_zero(q))) {
    fppol_set_zero(r);
    return;
  }

  // Check if input alias, and if so, work on a tmp space.
  fppol_t   t;
  fppol_ptr rr;
  if (r == p || r == q) { 
    fppol_init(t);
    rr = &t[0];
  }
  else
    rr = &r[0];

  // Do the multiplication.
  int kp = p->deg>>6;
  int kq = q->deg>>6;
  __fppol_realloc_lazy(rr, (kp+kq+2)<<6);
  memset(rr->limbs, 0, (kp+kq+2) * sizeof(fppol64_t));
#if defined(USE_F2) && defined(HAVE_GF2X)
  int cc = (sizeof(unsigned long) == sizeof(uint64_t))?1:2;
  gf2x_mul((unsigned long *)rr->limbs, 
           (unsigned long *)p->limbs, cc*(kp+1),
           (unsigned long *)q->limbs, cc*(kq+1));
#else
  for (int i = 0; i <= kp; ++i)
    for (int j = 0; j <= kq; ++j)
      __FP_MUL_128_64x64(rr->limbs[i+j+1], rr->limbs[i+j],
                         p ->limbs[i],     q ->limbs[j], add);
#endif
  rr->deg = p->deg + q->deg;

  // If input were aliases, do the final copy.
  if (r == p || r == q) {
    fppol_set(r, rr);
    fppol_clear(t);
  }
}


// Multiplication by an <sq>-term polynomial.
#define __DEF_FPPOL_MUL_xx(sq)                                           \
  void fppol_mul_##sq(fppol_ptr r, fppol_srcptr p, fppol##sq##_srcptr q) \
  {                                                                      \
    /* Multiplication by zero. */                                        \
    if (UNLIKELY(fppol_is_zero(p) || fppol##sq##_is_zero(q))) {          \
      fppol_set_zero(r);                                                 \
      return;                                                            \
    }                                                                    \
                                                                         \
    /* Check if input alias, and if so, work on a tmp space. */          \
    fppol_t   t;                                                         \
    fppol_ptr rr;                                                        \
    if (r == p) {                                                        \
      fppol_init(t);                                                     \
      rr = &t[0];                                                        \
    }                                                                    \
    else                                                                 \
      rr = &r[0];                                                        \
                                                                         \
    /* Do the multiplication. */                                         \
    int kp = p->deg>>6;                                                  \
    __fppol_realloc_lazy(rr, (kp+2)<<6);                                 \
    memset(rr->limbs, 0, (kp+2) * sizeof(fppol64_t));                    \
    for (int i = 0; i <= kp; ++i)                                        \
      __FP_MUL_128_64x##sq(rr->limbs[i+1], rr->limbs[i],                 \
                           p ->limbs[i],   q, add);                      \
    rr->deg = p->deg + fppol##sq##_deg(q);                               \
                                                                         \
    /* If input were aliases, do the final copy. */                      \
    if (r == p) {                                                        \
      fppol_set(r, rr);                                                  \
      fppol_clear(t);                                                    \
    }                                                                    \
  }


// Multiplication of <sp> x <sq>-term polynomials.
#define __DEF_FPPOL_MUL_xxxyy(sp, sq)                                    \
  void fppol_mul_##sp##x##sq(fppol_ptr          r,                       \
                             fppol##sp##_srcptr p, fppol##sq##_srcptr q) \
  { CAT(CAT(fppol,  __MUL_SIZE(sp, sq, )), _t) t;                        \
    CAT(CAT(fppol,  __MUL_SIZE(sp, sq, )), _mul_##sp##x##sq)(t, p, q);   \
    CAT(fppol_set_, __MUL_SIZE(sp, sq, ))(r, t); }


// Multiplication of 64 x <sq>-term polynomials.
#define __DEF_FPPOL_MUL_64xxx(sq)                                             \
  void fppol_mul_64x##sq(fppol_ptr r, fppol64_srcptr p, fppol##sq##_srcptr q) \
  { __fppol_realloc_lazy(r, 128);                                             \
    __FP_MUL_128_64x##sq(r->limbs[1], r->limbs[0], p, q, );                   \
    r->deg = 63 + sq-1;                                                       \
    __fppol_update_degree(r); }


__DEF_FPPOL_MUL_xx   (16)
__DEF_FPPOL_MUL_xx   (32)
__DEF_FPPOL_MUL_xx   (64)
__DEF_FPPOL_MUL_xxxyy(16, 16)
__DEF_FPPOL_MUL_xxxyy(32, 16)
__DEF_FPPOL_MUL_xxxyy(32, 32)
__DEF_FPPOL_MUL_64xxx(16)
__DEF_FPPOL_MUL_64xxx(32)
__DEF_FPPOL_MUL_64xxx(64)

#undef __DEF_FPPOL_MUL_xx
#undef __DEF_FPPOL_MUL_xxxyy
#undef __DEF_FPPOL_MUL_64xxx
