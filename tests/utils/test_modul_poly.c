#include "cado.h"
#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "modul_poly.h"
#include "gmp_aux.h"
#include "tests_common.h"
#include "cado_poly.h"

void
test_modul_poly_is_irreducible (unsigned long iter)
{
  modul_poly_t f;
  modulusul_t p;
  int d, i, irred, n;
  unsigned long q;
  residueul_t r[MAXDEGREE];

  /* first try some hard-coded polynomials */
  modul_poly_init (f, MAXDEGREE);

  modul_initmod_ul (p, 3);
  modul_set_ul (f->coeff[0], 1, p);
  f->degree = 0;
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) != 0);

  modul_initmod_ul (p, 3);
  modul_set_ul (f->coeff[0], 0, p);
  f->degree = -1;
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) != 0);

  modul_initmod_ul (p, 3);
  modul_set_ul (f->coeff[2], 2, p);
  modul_set_ul (f->coeff[1], 1, p);
  modul_set_ul (f->coeff[0], 2, p);
  f->degree = 2;
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) == 0);

  modul_initmod_ul (p, 5);
  modul_set_ul (f->coeff[2], 2, p);
  modul_set_ul (f->coeff[1], 1, p);
  modul_set_ul (f->coeff[0], 2, p);
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) == 0);

  modul_initmod_ul (p, 7);
  modul_set_ul (f->coeff[2], 2, p);
  modul_set_ul (f->coeff[1], 1, p);
  modul_set_ul (f->coeff[0], 2, p);
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) != 0);

  modul_initmod_ul (p, 11);
  modul_set_ul (f->coeff[2], 2, p);
  modul_set_ul (f->coeff[1], 1, p);
  modul_set_ul (f->coeff[0], 2, p);
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) != 0);

  while (iter--)
    {
      d = 1 + lrand48 () % (MAXDEGREE - 1);
      q = lrand48 ();
      q = ulong_nextprime (q);
      /* modul_poly_cantor_zassenhaus only works for odd primes */
      q += (q == 2);
      modul_initmod_ul (p, q);
      for (i = 0; i <= d; i++)
        modul_set_ul (f->coeff[i], lrand48 (), p);
      while (modul_is0 (f->coeff[d], p))
        modul_set_ul (f->coeff[d], lrand48 (), p);
      f->degree = d;
      irred = modul_poly_is_irreducible (f, p);
      if (d == 1)
        ASSERT_ALWAYS(irred != 0);
      if (d == 1 || (d == 2 && irred == 0))
        {
          n = modul_poly_cantor_zassenhaus (r, f, p);
          ASSERT_ALWAYS(n == d);
        }
      modul_clearmod (p);
    }
  modul_poly_clear (f);
}

void
test_modul_poly_roots_ulong (unsigned long iter)
{
  unsigned long r[MAXDEGREE];
  mpz_t f[MAXDEGREE + 1];
  int d, i, n;
  modulusul_t p;
  residueul_t y, x;
  modul_poly_t fp;

  for (i = 0; i <= MAXDEGREE; i++)
    mpz_init (f[i]);

  /* hard-coded examples */
  mpz_set_ui (f[0], 0);
  mpz_set_ui (f[1], 0);
  mpz_set_ui (f[2], 1);
  d = 2;
  modul_initmod_ul (p, 113);
  n = modul_poly_roots_ulong (r, f, d, p);
  ASSERT_ALWAYS(n == 1 && r[0] == 0);
  modul_clearmod (p);

  while (iter--)
    {
      d = 1 + lrand48 () % (MAXDEGREE - 1);
      for (i = 0; i <= d; i++)
        mpz_urandomb (f[i], state, 64);
      while (mpz_cmp_ui (f[d], 0) == 0)
        mpz_urandomb (f[d], state, 64);
      modul_initmod_ul (p, ulong_nextprime (lrand48 ()));
      n = modul_poly_roots_ulong (r, f, d, p);
      ASSERT_ALWAYS(0 <= n && n <= d);
      modul_poly_init (fp, d);
      modul_poly_set_mod (fp, f, d, p);
      if (n > 0 && fp->degree > 1)
        ASSERT_ALWAYS(modul_poly_is_irreducible (fp, p) == 0);
      /* if n=0, f might be irreducible or not mod p (product of two
         degree-2 factors for example),
         if d=1, f is irreducible */
      if (d == 1 || (d <= 3 && n == 0))
        ASSERT_ALWAYS(modul_poly_is_irreducible(fp, p) != 0);
      for (i = 0; i < n; i++)
        {
          modul_set_ul (x, r[i], p);
          modul_poly_eval (y, fp, x, p);
          ASSERT_ALWAYS(modul_is0 (y, p));
        }
      modul_poly_clear (fp);
      modul_clearmod (p);
    }
  for (i = 0; i <= MAXDEGREE; i++)
    mpz_clear (f[i]);
}

int
main (int argc, const char *argv[])
{
  unsigned long iter = 1000;

  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter (&iter);
  test_modul_poly_is_irreducible (iter);
  test_modul_poly_roots_ulong (iter);
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}

