#include "fppol_test.h"
#include "fppol_set.h"



/* Multiprecision polynomials.
 *****************************************************************************/

// Test if equal.
int fppol_eq(fppol_srcptr p, fppol_srcptr q)
{
  if (p->deg != q->deg)
    return 0;
  int k;
  for (k = p->deg>>6; k >= 0 && fppol64_eq(p->limbs[k], q->limbs[k]); --k);
  return k < 0;
}


// Test if monic.
int fppol_is_monic(fppol_srcptr p)
{
  if (UNLIKELY(p->deg == -1)) return 0;
  fp_t lc;
  fppol_get_coeff(lc, p, p->deg);
  return fp_is_one(lc);
}


// Test if valid representation.
int fppol_is_valid(fppol_srcptr p)
{
  int k;
  for (k = p->deg>>6; k >= 0 && fppol64_is_valid(p->limbs[k]); --k);
  return k < 0;
}
