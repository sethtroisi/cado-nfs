#include <string.h>

#include "fppol_arith.h"
#include "fppol_internal.h"



/* Fixed-size polynomials.
 *****************************************************************************/

// Degree.
// By convention, deg(0) = -1.
#define __DEF_FPPOLxx_DEG(sz)                \
  int fppol##sz##_deg(fppol##sz##_srcptr p)  \
  {                                          \
    uint##sz##_t t = fppol##sz##_fold_or(p); \
    if (!t) return -1;                       \
    int d = 0;                               \
    for (unsigned k = sz; k >>= 1; )         \
      if (t >> k) d += k, t >>= k;           \
    return d;                                \
  }

__DEF_FPPOLxx_DEG(16)
__DEF_FPPOLxx_DEG(32)
__DEF_FPPOLxx_DEG(64)

#undef __DEF_FPPOLxx_DEG



/* Multiprecision polynomials.
 *****************************************************************************/

// Set to zero.
void fppol_set_zero(fppol_ptr r)
{
  __fppol_realloc_lazy(r, 0);
  r->deg = -1;
}


// Set to one.
void fppol_set_one(fppol_ptr r)
{
  __fppol_realloc_lazy(r, 1);
  fppol64_set_one(r->limb[0]);
  r->deg = 0;
}


// Set to t^i.
void fppol_set_ti(fppol_ptr r, unsigned i)
{
  __fppol_realloc_lazy(r, i+1);
  memset(r->limb, 0, (i>>6) * sizeof(fppol64_t));
  fppol64_set_ti(r->limb[i>>6], i&0x3f);
  r->deg = i;
}


// Set to another polynomial.
void fppol_set(fppol_ptr r, fppol_srcptr p)
{
  __fppol_realloc_lazy(r, p->deg+1);
  for (int k = 0; k <= p->deg>>6; ++k)
    fppol64_set(r->limb[k], p->limb[k]);
  r->deg = p->deg;
}


// Set to an <sz>-term polynomial.
#define __DEF_FPPOL_SET_xx(sz)                           \
  void fppol_set_##sz(fppol_ptr r, fppol##sz##_srcptr p) \
  { __fppol_realloc_lazy(r, sz);                         \
    fppol64_set_##sz(r->limb[0], p);                     \
    r->deg = fppol##sz##_deg(p); }

__DEF_FPPOL_SET_xx(16)
__DEF_FPPOL_SET_xx(32)
__DEF_FPPOL_SET_xx(64)

#undef __DEF_FPPOL_SET_xx


// Set to something coded in an uint64_t.
// Not sure this is meaningful. Will probably disappear soon.
void fppol_set_ui(fppol_ptr r, uint64_t x)
{
  __fppol_realloc_lazy(r, 1);
  // fppol64_set_ui(r->limb[0], ui);
  fppol64_set_zero(r->limb[0]);
  r->limb[0][0] = x;
  r->deg = fppol64_deg(r->limb[0]);
}


// Set degree-i coefficient.
void fppol_set_coeff(fppol_ptr r, fp_srcptr x, unsigned i)
{
  int ii = i;
  if (ii > r->deg && fp_is_zero(x)) return;
  if (ii>>6 > r->deg>>6) {
    __fppol_realloc_lazy(r, ii+1);
    memset(&r->limb[(r->deg>>6)+1], 0,
           ((ii>>6)-(r->deg>>6)) * sizeof(fppol64_t));
  }
  fppol64_set_coeff(r->limb[i>>6], x, i&0x3f);
  if (ii > r->deg)
    r->deg = ii;
  else if (fp_is_zero(x) && ii == r->deg)
    __fppol_update_degree(r);
}


// Test if equal.
int fppol_eq(fppol_srcptr p, fppol_srcptr q)
{
  if (p->deg != q->deg)
    return 0;
  int k;
  for (k = p->deg>>6; k >= 0 && fppol64_eq(p->limb[k], q->limb[k]); --k);
  return k < 0;
}


// Test if valid representation.
int fppol_is_valid(fppol_srcptr p)
{
  int k;
  for (k = p->deg>>6; k >= 0 && fppol64_is_valid(p->limb[k]); --k);
  return k < 0;
}


// Opposite.
void fppol_opp(fppol_ptr r, fppol_srcptr p)
{
  __fppol_realloc_lazy(r, p->deg+1);
  for (int k = 0; k <= p->deg>>6; ++k)
    fppol64_opp(r->limb[k], p->limb[k]);
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
    fppol64_add(r->limb[k], pp->limb[k], qq->limb[k]);
  for (; k <= qq->deg>>6; ++k)
    fppol64_set(r->limb[k], qq->limb[k]);
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
    fppol64_sub(r->limb[k], p->limb[k], q->limb[k]);
  for (; k <= qq->deg>>6; ++k) {
    if (opp) fppol64_set(r->limb[k], qq->limb[k]);
    else     fppol64_opp(r->limb[k], qq->limb[k]);
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
    fppol64_smul(r->limb[k], p->limb[k], x);
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
    fppol64_sdiv(r->limb[k], p->limb[k], x);
  r->deg = p->deg;
}


// Coefficient-wise OR.
static inline
void __fppol64_or(fppol64_ptr r, fppol64_srcptr p, fppol64_srcptr q)
{
  for (unsigned k = 0; k < __FP_BITS; ++k)
    r[k] = p[k] | q[k];
}


// Left shift.
void fppol_shl(fppol_ptr r, fppol_srcptr p, unsigned i)
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
  fppol64_set_zero(r->limb[kp+d+1]);
  for (int k = kp; k >= 0; --k) {
      fppol64_shr(t,              p->limb[k], 64-(i&0x3f));
    __fppol64_or (r->limb[k+d+1], r->limb[k+d+1], t);
      fppol64_shl(r->limb[k+d],   p->limb[k],     i&0x3f);
  }
  memset(r->limb, 0, d * sizeof(fppol64_t));
  r->deg = p->deg+i;
}


// Right shift.
void fppol_shr(fppol_ptr r, fppol_srcptr p, unsigned i)
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
  fppol64_shr(r->limb[0], p->limb[d], i&0x3f);
  for (int k = d+1; k <= kp; ++k) {
      fppol64_shl(t,              p->limb[k], 64-(i&0x3f));
    __fppol64_or (r->limb[k-d-1], r->limb[k-d-1], t);
      fppol64_shr(r->limb[k-d],   p->limb[k],     i&0x3f);
  }
  r->deg = p->deg-i;
}
