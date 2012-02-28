// Always include "fppol.h" first so that everything gets loaded in the correct
// order.
#include "fppol.h"

#ifndef __FPPOL_INTERNAL_H__
#define __FPPOL_INTERNAL_H__



// Swap two base field elements.
static inline
void __fp_swap(fp_ptr x, fp_ptr y)
{ fp_t t; fp_set(t, x); fp_set(x, y); fp_set(y, t); }

// Swap two <sz>-term polynomials.
// Generic prototype:
//   void __fppol<sz>_swap(fppol<sz>_ptr p, fppol<sz>_ptr q);
#define __DEF_FPPOLxx_SWAP(sz)                                  \
  static inline                                                 \
  void __fppol##sz##_swap(fppol##sz##_ptr p, fppol##sz##_ptr q) \
  { fppol##sz##_t t;       fppol##sz##_set(t, p);               \
    fppol##sz##_set(p, q); fppol##sz##_set(q, t); }

__DEF_FPPOLxx_SWAP(16)
__DEF_FPPOLxx_SWAP(32)
__DEF_FPPOLxx_SWAP(64)

#undef __DEF_FPPOLxx_SWAP

// Swap two multiprecision polynomials.
// /!\ Shallow copy only.
static inline
void __fppol_swap(fppol_ptr p, fppol_ptr q)
{ __fppol_struct t = *p; *p = *q; *q = t; }

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
