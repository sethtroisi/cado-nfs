#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "ffspol.h"



/* Miscellaneous functions for internal use only.
 *****************************************************************************/

// Reallocate memory with space for n coefficients.
static
void __ffspol_realloc(ffspol_ptr r, unsigned n)
{
  for (unsigned i = n; i < r->alloc; ++i)
    fppol_clear(r->coeffs[i]);
  r->coeffs = (fppol_t *)realloc(r->coeffs, n * sizeof(fppol_t));
  ASSERT_ALWAYS(!n || r->coeffs != NULL);
  for (unsigned i = r->alloc; i < n; ++i)
    fppol_init(r->coeffs[i]);
  r->alloc = n;
}


// Reallocate memory only when expanding.
static inline
void __ffspol_realloc_lazy(ffspol_ptr r, unsigned n)
{ if (n > r->alloc) __ffspol_realloc(r, n); }


// Update the degree: assume that the given degree is an upper bound, and 
// check the nullity of the coefficients.
// /!\ This function works by side-effects.
static inline
void __ffspol_update_degree(ffspol_ptr r)
{ for (; r->deg >= 0 && fppol_is_zero(r->coeffs[r->deg]); --r->deg); }



/* Initialization/destruction.
 *****************************************************************************/

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
  r->alloc  = n;
  r->coeffs = (fppol_t *)malloc(r->alloc * sizeof(fppol_t));
  ASSERT_ALWAYS(!n || r->coeffs != NULL);
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


/* Input / output.
 *****************************************************************************/

// Conversion to string.
char *ffspol_get_str(char *str, ffspol_srcptr p)
{
  if (str == NULL) {
    str = (char *)malloc((ffspol_strlen(p)+1) * sizeof(char));
    ASSERT_ALWAYS(str != NULL);
  }
  strcpy(str, "0");
  char *ptr = str;
  for (int i = 0; i <= p->deg; ++i) {
    if (i) *ptr++ = ',';
    fppol_get_str(ptr,  p->coeffs[i]);
    ptr += fppol_strlen(p->coeffs[i]);
  }
  return str;
}


// Conversion from string.
// Return 1 if successful.
int ffspol_set_str(ffspol_ptr r, const char *str)
{
  unsigned n = 0, lmax = 0;
  for (unsigned i = 0, j = 0; i <= strlen(str); ++i) {
    if (str[i] && str[i] != ',') continue;
    ++n;
    if (lmax < i-j) lmax = i-j;
    j = i+1;
  }
  __ffspol_realloc_lazy(r, n);
  char *buf = (char *)malloc((lmax+1) * sizeof(char));
  ASSERT_ALWAYS(buf != NULL);
  r->deg = -1;
  for (unsigned i = 0, j = 0; i <= strlen(str); ++i) {
    if (str[i] && str[i] != ',') continue;
    strncpy(buf, str+j, i-j); buf[i-j] = '\0';
    if (!fppol_set_str(r->coeffs[++r->deg], buf)) return 0;
    j = i+1;
  }
  __ffspol_update_degree(r);
  return 1;
}


// Output to stream.
void ffspol_out(FILE *f, ffspol_srcptr p)
{
  char *buf = ffspol_get_str(NULL, p);
  fprintf(f == NULL ? stdout : f, "%s", buf);
  free(buf);
}


// Input from stream.
// Return 1 if successful.
int ffspol_inp(ffspol_ptr r, FILE *f)
{
  if (f == NULL) f = stdin;
  int c = ',';
  for (r->deg = -1; c == ','; c = getc(f))
    if (!fppol_inp(r->coeffs[++r->deg], f)) return 0;
  ungetc(c, f);
  __ffspol_update_degree(r);
  return 1;
}
