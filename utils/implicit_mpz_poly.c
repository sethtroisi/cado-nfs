#include "cado.h"
#include <stdlib.h> /* for malloc and free */
#include "implicit_mpz_poly.h"
#include "gmp_aux.h"


/* put in c the content of f */
void
mp_poly_content (mpz_t c, mpz_t *f, const int d)
{
  int i;

  mpz_set (c, f[0]);
  for (i = 1; i <= d; i++)
    mpz_gcd (c, c, f[i]);
  mpz_abs (c, c);
}
