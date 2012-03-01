#ifndef __FFSPOL_H__
#define __FFSPOL_H__

#include "fppol.h"



/* Type definitions.
 *****************************************************************************/

// Polynomials over GF(p)[x,t], represented as polynomials in x whose
// coefficients lie in GF(p)[t] (using the multiprecision type fppol_t).
// By convention, deg(0) = -1.
typedef struct {
  int       deg;
  unsigned  alloc;
  fppol_t  *coeffs;
} __ffspol_struct;

// Type and pointer shorthands.
// - ffspol_t:      type of polynomials.
// - ffspol_ptr:    read/write pointer (internal type)
// - ffspol_srcptr: read-only  pointer (internal type)
typedef       __ffspol_struct  ffspol_t[1];
typedef       __ffspol_struct *ffspol_ptr;
typedef const __ffspol_struct *ffspol_srcptr;



/* Initialization/destruction.
 *****************************************************************************/

// Initialize polynomial.
void ffspol_init(ffspol_ptr r);

// Initialize a NULL-terminated list of polynomials.
void ffspol_inits(ffspol_ptr r, ...);

// Initialize polynomial with space for n coefficients.
void ffspol_init2(ffspol_ptr r, unsigned n);

// Free polynomial.
void ffspol_clear(ffspol_ptr r);

// Free a NULL-terminated list of polynomials.
void ffspol_clears(ffspol_ptr r, ...);

#endif   /* __FFSPOL_H__ */
