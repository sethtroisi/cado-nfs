#include <string.h>

#include "fppol_set.h"
#include "fppol_internal.h"



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
  fppol64_set_one(r->limbs[0]);
  r->deg = 0;
}


// Set to t^i.
void fppol_set_ti(fppol_ptr r, unsigned i)
{
  __fppol_realloc_lazy(r, i+1);
  memset(r->limbs, 0, (i>>6) * sizeof(fppol64_t));
  fppol64_set_ti(r->limbs[i>>6], i&0x3f);
  r->deg = i;
}


// Set to another polynomial.
void fppol_set(fppol_ptr r, fppol_srcptr p)
{
  __fppol_realloc_lazy(r, p->deg+1);
  for (int k = 0; k <= p->deg>>6; ++k)
    fppol64_set(r->limbs[k], p->limbs[k]);
  r->deg = p->deg;
}


// Set to an <sz>-term polynomial.
#define __DEF_FPPOL_SET_xx(sz)                           \
  void fppol_set_##sz(fppol_ptr r, fppol##sz##_srcptr p) \
  { __fppol_realloc_lazy(r, sz);                         \
    fppol64_set_##sz(r->limbs[0], p);                    \
    r->deg = fppol##sz##_deg(p); }

__DEF_FPPOL_SET_xx(16)
__DEF_FPPOL_SET_xx(32)
__DEF_FPPOL_SET_xx(64)

#undef __DEF_FPPOL_SET_xx


// Set degree-i coefficient.
void fppol_set_coeff(fppol_ptr r, fp_srcptr x, unsigned i)
{
  int ii = i;
  if (ii > r->deg && fp_is_zero(x)) return;
  if (ii>>6 > r->deg>>6) {
    __fppol_realloc_lazy(r, ii+1);
    memset(&r->limbs[(r->deg>>6)+1], 0,
           ((ii>>6)-(r->deg>>6)) * sizeof(fppol64_t));
  }
  fppol64_set_coeff(r->limbs[i>>6], x, i&0x3f);
  if (ii > r->deg)
    r->deg = ii;
  else if (fp_is_zero(x) && ii == r->deg)
    __fppol_update_degree(r);
}
