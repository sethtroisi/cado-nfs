#include <stdlib.h>

#include "ffspol.h"



// Reallocate memory with space for n coefficients.
void __ffspol_realloc(ffspol_ptr r, unsigned n)
{
  for (unsigned i = n; i < r->alloc; ++i)
    fppol_clear(r->coeffs[i]);
  r->coeffs = realloc(r->coeffs, n * sizeof(fppol_t));
  ASSERT_ALWAYS(!n || r->coeffs != NULL);
  for (unsigned i = r->alloc; i < n; ++i)
    fppol_init(r->coeffs[i]);
  r->alloc = n;
}


// Reallocate memory only when expanding.
static inline
void __ffspol_realloc_lazy(ffspol_ptr r, unsigned n)
{ if (n > r->alloc) __ffspol_realloc(r, n); }


// Initialize polynomial.
void ffspol_init(ffspol_ptr r)
{
  r->alloc  = 0;
  r->coeffs = NULL;
}


// Initialize a NULL-terminated list of polynomials.
void ffspol_inits(ffspol_ptr r, ...)
{
  va_list ap;
  for (va_start(ap, r); r != NULL; r = va_arg(ap, ffspol_ptr))
    ffspol_init(r);
  va_end(ap);
}


// Initialize polynomial with space for n coefficients.
void ffspol_init2(ffspol_ptr r, unsigned n)
{
  r->alloc = n;
  r->limbs = malloc(r->alloc * sizeof(fppol_t));
  ASSERT_ALWAYS(!n || r->limbs != NULL);
}


// Free polynomial.
void ffspol_clear(ffspol_ptr r)
{
  for (unsigned i = 0; i < r->alloc; ++i)
    fppol_clear(r->coeffs[i]);
  free(r->coeffs);
}


// Free a NULL-terminated list of polynomials.
void ffspol_clears(ffspol_ptr r, ...)
{
  va_list ap;
  for (va_start(ap, r); r != NULL; r = va_arg(ap, ffspol_ptr))
    ffspol_clear(r);
  va_end(ap);
}
