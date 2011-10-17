#include "cado.h"
#include <stdlib.h>
#include "macros.h"
#include "discriminant.h"

/* D <- |discriminant (f)| = |resultant(f,diff(f,x))/lc(f)| */
void
discriminant (mpz_t D, mpz_t *f0, const int d)
{
  mpz_t *f, *g, num, den;
  int df, dg, i, s;

  f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
  ASSERT_ALWAYS(f != NULL);
  for (i = 0; i < d + 1; i++)
    mpz_init (f[i]);

  g = (mpz_t*) malloc (d * sizeof (mpz_t));
  ASSERT_ALWAYS(g != NULL);
  for (i = 0; i < d; i++)
    mpz_init (g[i]);

  for (i = 0; i <= d; i++)
    mpz_set (f[i], f0[i]);
  df = d;
  for (i = 0; i < d; i++)
    mpz_mul_ui (g[i], f0[i + 1], i + 1);
  dg = d - 1;
  mpz_init_set_ui (num, 1UL);
  mpz_init_set_ui (den, 1UL);
  while (dg > 0)
    {
      s = df;
      while (df >= dg)
	{
	  /* f <- f * lc(g), except f[df] since we'll divide it afterwards */
	  for (i = 0; i < df; i++)
	    mpz_mul (f[i], f[i], g[dg]);
	  s -= dg;
	  for (i = 0; i < dg; i++)
	    mpz_submul (f[i + df - dg], g[i], f[df]);
	  df --;
	  /* normalize f */
	  while (df > 0 && mpz_cmp_ui (f[df], 0) == 0)
	    df --;
	}
      /* den <- den * lc(g)^deg(f) */
      s -= df;
      if (s > 0)
	for (i = 0; i < s; i++)
	  mpz_mul (num, num, g[dg]);
      else
	for (i = 0; i < -s; i++)
	  mpz_mul (den, den, g[dg]);
      /* swap f and g */
      for (i = 0; i <= dg; i++)
	mpz_swap (f[i], g[i]);
      i = df;
      df = dg;
      dg = i;
    }
  /* num/den*g^deg(f) */
  mpz_pow_ui (g[0], g[0], df);
  mpz_mul (g[0], g[0], num);
  ASSERT (mpz_divisible_p (g[0], den));
  mpz_divexact (g[0], g[0], den);
  /* divide by lc(f0) to get the discriminant */
  ASSERT (mpz_divisible_p (g[0], f0[d]));
  mpz_divexact (D, g[0], f0[d]);
  mpz_clear (num);
  mpz_clear (den);

  for (i = 0; i < d + 1; i++)
    mpz_clear (f[i]);
  free (f);

  for (i = 0; i < d; i++)
    mpz_clear (g[i]);
  free (g);
}
