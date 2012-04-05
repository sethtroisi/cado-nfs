#include <string.h>

#include "fppol_arith.h"
#include "fppol_internal.h"
#include "macros.h"



/* Multiprecision polynomials.
 *****************************************************************************/

// Opposite.
void fppol_opp(fppol_ptr r, fppol_srcptr p)
{
  __fppol_realloc_lazy(r, p->deg+1);
  for (int k = 0; k <= p->deg>>6; ++k)
    fppol64_opp(r->limbs[k], p->limbs[k]);
  r->deg = p->deg;
}


// Addition.
void fppol_add(fppol_ptr r, fppol_srcptr p, fppol_srcptr q)
{
  int need_update = (p->deg == q->deg), k;

  // Use aliases to have deg(qq) >= deg(pp).
  fppol_srcptr pp, qq;
  if (p->deg > q->deg) { pp = q; qq = p; }
  else                 { pp = p; qq = q; }

  __fppol_realloc_lazy(r, qq->deg+1);
  for (k = 0; k <= pp->deg>>6; ++k)
    fppol64_add(r->limbs[k], pp->limbs[k], qq->limbs[k]);
  for (; k <= qq->deg>>6; ++k)
    fppol64_set(r->limbs[k], qq->limbs[k]);
  r->deg = qq->deg;
  if (need_update)
    __fppol_update_degree(r);
}


// Subtraction.
void fppol_sub(fppol_ptr r, fppol_srcptr p, fppol_srcptr q)
{
  int need_update = (p->deg == q->deg), opp, k;

  // Use aliases to have deg(qq) >= deg(pp).
  fppol_srcptr pp, qq;
  if (p->deg > q->deg) { pp = q; qq = p; opp = 1; }
  else                 { pp = p; qq = q; opp = 0; }

  __fppol_realloc_lazy(r, qq->deg+1);
  for (k = 0; k <= pp->deg>>6; ++k)
    fppol64_sub(r->limbs[k], p->limbs[k], q->limbs[k]);
  for (; k <= qq->deg>>6; ++k) {
    if (opp) fppol64_set(r->limbs[k], qq->limbs[k]);
    else     fppol64_opp(r->limbs[k], qq->limbs[k]);
  }
  r->deg = qq->deg;
  if (need_update)
    __fppol_update_degree(r);
}


// Scalar multiplication.
void fppol_smul(fppol_ptr r, fppol_srcptr p, fp_srcptr x)
{
  if (fp_is_zero(x)) {
    fppol_set_zero(r);
    return;
  }
  __fppol_realloc_lazy(r, p->deg+1);
  for (int k = 0; k <= p->deg>>6; ++k)
    fppol64_smul(r->limbs[k], p->limbs[k], x);
  r->deg = p->deg;
}


// Scalar division.
void fppol_sdiv(fppol_ptr r, fppol_srcptr p, fp_srcptr x)
{
  if (fp_is_zero(x)) {
    fppol_set_zero(r);
    return;
  }
  __fppol_realloc_lazy(r, p->deg+1);
  for (int k = 0; k <= p->deg>>6; ++k)
    fppol64_sdiv(r->limbs[k], p->limbs[k], x);
  r->deg = p->deg;
}


// Scalar inversion.
void fppol_sinv(fppol_ptr r, fppol_srcptr p)
{
  __fppol_realloc_lazy(r, p->deg+1);
  for (int k = 0; k <= p->deg>>6; ++k)
    fppol64_sinv(r->limbs[k], p->limbs[k]);
  r->deg = p->deg;
}


// Addition without overlapping coefficients (i.e., coefficient-wise OR).
void fppol_add_disjoint(fppol_ptr r, fppol_srcptr p, fppol_srcptr q)
{
  // Use aliases to have deg(qq) >= deg(pp).
  fppol_srcptr pp, qq;
  if (p->deg > q->deg) { pp = q; qq = p; }
  else                 { pp = p; qq = q; }

  __fppol_realloc_lazy(r, qq->deg+1);
  int k;
  for (k = 0; k <= pp->deg>>6; ++k)
    fppol64_add_disjoint(r->limbs[k], pp->limbs[k], qq->limbs[k]);
  for (; k <= qq->deg>>6; ++k)
    fppol64_set(r->limbs[k], qq->limbs[k]);
  r->deg = qq->deg;
}


// Multiplication by t^i (i.e., left shift).
void fppol_mul_ti(fppol_ptr r, fppol_srcptr p, unsigned i)
{
  // If polynomial is zero, shifting's easy.
  if (UNLIKELY(fppol_is_zero(p))) {
    fppol_set_zero(r);
    return;
  }

  // Shift limb by limb, most significant limbs first.
  int kp = p->deg>>6, d = i>>6;
  fppol64_t t;
  __fppol_realloc_lazy(r, (kp+d+2)<<6);
  fppol64_set_zero(r->limbs[kp+d+1]);
  for (int k = kp; k >= 0; --k) {
    fppol64_div_ti      (t,               p->limbs[k], 64-(i&0x3f));
    fppol64_add_disjoint(r->limbs[k+d+1], r->limbs[k+d+1], t);
    fppol64_mul_ti      (r->limbs[k+d],   p->limbs[k],     i&0x3f);
  }
  memset(r->limbs, 0, d * sizeof(fppol64_t));
  r->deg = p->deg+i;
}


// Quotient of division by t^i (i.e., right shift).
void fppol_div_ti(fppol_ptr r, fppol_srcptr p, unsigned i)
{
  // If polynomial is too small, shifting's easy.
  if (UNLIKELY(p->deg < (signed)i)) {
    fppol_set_zero(r);
    return;
  }

  // Shift limb by limb, least significant limbs first.
  int kp = p->deg>>6, d = i>>6;
  fppol64_t t;
  __fppol_realloc_lazy(r, (kp-d+1)<<6);
  fppol64_div_ti(r->limbs[0], p->limbs[d], i&0x3f);
  for (int k = d+1; k <= kp; ++k) {
    fppol64_mul_ti      (t,               p->limbs[k], 64-(i&0x3f));
    fppol64_add_disjoint(r->limbs[k-d-1], r->limbs[k-d-1], t);
    fppol64_div_ti      (r->limbs[k-d],   p->limbs[k],     i&0x3f);
  }
  r->deg = p->deg-i;
}


// Remainder modulo t^i (i.e., i first coefficients).
void fppol_mod_ti(fppol_ptr r, fppol_srcptr p, unsigned i)
{
  // If polynomial is too small, masking's easy.
  if (UNLIKELY(p->deg < (signed)i)) {
    fppol_set(r, p);
    return;
  }

  __fppol_realloc_lazy(r, i);
  unsigned k;
  for (k = 0; k < i>>6; ++k)
    fppol64_set(r->limbs[k], p->limbs[k]);
  fppol64_mod_ti(r->limbs[k], p->limbs[k], i&0x3f);
}
