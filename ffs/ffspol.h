#ifndef __FFSPOL_H__
#define __FFSPOL_H__

#include "types.h"



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


void ffspol_eval(fppol_ptr z, ffspol_srcptr f, fppol_srcptr x);
void ffspol_eval_diff(fppol_ptr z, ffspol_srcptr f, fppol_srcptr x);
void ffspol_mul(ffspol_ptr z, ffspol_srcptr x, ffspol_srcptr y);
void ffspol_add(ffspol_ptr r, ffspol_srcptr p, ffspol_srcptr q);
void ffspol_smul(ffspol_ptr z, ffspol_srcptr x, fppol_srcptr p);


#endif   /* __FFSPOL_H__ */
