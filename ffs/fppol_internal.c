#include <stdlib.h>

#include "fppol_internal.h"



// Reallocate memory with space for n terms, ie degree n-1.
void __fppol_realloc(fppol_ptr r, unsigned n)
{
  r->alloc = (n+63)>>6;
  r->limb  = realloc(r->limb, r->alloc * sizeof(fppol64_t));
  ASSERT_ALWAYS(!n || r->limb != NULL);
}


// Update the degree: assume that the given degree is an upper bound, and 
// check the nullity of the coefficients.
// /!\ This function works by side-effects.
void __fppol_update_degree(fppol_ptr r)
{
  int k;
  for (k = r->deg>>6; k >= 0 && fppol64_is_zero(r->limb[k]); --k);
  r->deg = k >= 0 ? (k<<6) + fppol64_deg(r->limb[k]) : -1;
}
