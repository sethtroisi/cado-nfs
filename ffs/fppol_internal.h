// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_INTERNAL_H__
#define __FPPOL_INTERNAL_H__



// Reallocate memory with space for n terms, ie degree n-1.
void __fppol_realloc(fppol_ptr r, unsigned n);

// Reallocate memory only when expanding.
static inline
void __fppol_realloc_lazy(fppol_ptr r, unsigned n)
{ if (n > r->alloc<<6) __fppol_realloc(r, n); }

// Update the degree: assume that the given degree is an upper bound, and 
// check the nullity of the coefficients.
// /!\ This function works by side-effects.
void __fppol_update_degree(fppol_ptr r);

#endif  /* __FPPOL_INTERNAL_H__ */
