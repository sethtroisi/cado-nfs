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



/* Input / output.
 *****************************************************************************/

// Length of string representation, not counting the null-terminator.
static inline
size_t ffspol_strlen(ffspol_srcptr p)
{
  unsigned l = p->deg < 0 ? 1 : p->deg;
  for (int i = 0; i <= p->deg; ++i)
    l += fppol_strlen(p->coeffs[i]);
  return l;
}

// Conversion to string.
char *ffspol_get_str(char *str, ffspol_srcptr p);

// Conversion from string.
// Return 1 if successful.
int ffspol_set_str(ffspol_ptr r, const char *str);

// Output to stream.
void ffspol_out(FILE *f, ffspol_srcptr p);

// Input from stream.
// Return 1 if successful.
int ffspol_inp(ffspol_ptr r, FILE *f);

#endif   /* __FFSPOL_H__ */
