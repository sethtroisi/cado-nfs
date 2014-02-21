#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "rootfinder.h"
#include "portability.h"
#include "test_iter.h"
#include "tests_common.h"
#include "mpz_poly.h"
#include "macros.h"
#include "gmp_aux.h"

int cmp(mpz_t * a, mpz_t * b)
{
    return mpz_cmp(*a, *b);
}

/* Check roots of polynomial ff[d]*x^d + ... + ff[1]*x + ff[0] mod pp.
   If pp is the empty string, generate a random prime (and then generate
   random coefficients for ff[]).
   If nroots <> -1, it is the expected number of roots. */
void
test (int d, const char *pp, const char *ff[], int nroots)
{
  mpz_t p, *f, *r, v;
  int i, n, n0, n1;
  mpz_poly_t F;

  mpz_init (p);
  if (strlen (pp) > 0)
    mpz_set_str (p, pp, 0);
  else
    {
      mpz_urandomb (p, state, 65);
      mpz_nextprime (p, p);
    }
  ASSERT_ALWAYS(mpz_probab_prime_p (p, 5));
  mpz_init (v);
  f = (mpz_t *) malloc ((d + 1) * sizeof(mpz_t));
  r = (mpz_t *) malloc ((d + 1) * sizeof(mpz_t));
  F->coeff = f;
  F->deg = d;
  for (i = 0; i <= d; i++)
    {
      mpz_init (f[i]);
      if (strlen (pp) > 0)
        mpz_set_str (f[i], ff[d - i], 0);
      else
        {
          do mpz_urandomb (f[i], state, 65);
          while (i == d && mpz_cmp_ui (f[i], 0) == 0);
        }
      mpz_init (r[i]);
    }
  if (sizeof(long)==4)
    {
      gmp_printf ("Testing polynomial of degree %d modulo p=%Zd:", d, p);
      for (i = d; i >= 0; i--)
        gmp_printf (" %Zd*x^%d", f[i], i);
      printf ("\n");
      fflush (stdout);
    }
  n0 = poly_roots (NULL, f, d, p);
  if (sizeof(long)==4)
    {
      printf ("n0=%d\n", n0);
      fflush (stdout);
    }
  if (mpz_sizeinbase (p, 2) <= 64)
    {
      n1 = poly_roots_uint64 (NULL, f, d, mpz_get_uint64 (p));
      if (sizeof(long)==4)
        {
          printf ("n1=%d\n", n1);
          fflush (stdout);
        }
      ASSERT_ALWAYS(n1 == n0);
    }
  n = poly_roots (r, f, d, p);
  if (sizeof(long)==4)
    {
      printf ("n=%d\n", n);
      fflush (stdout);
    }
  ASSERT_ALWAYS(n == n0);
  if (nroots != -1)
    ASSERT_ALWAYS(n == nroots);
  for (i = 0; i < n; i++)
    {
      mpz_poly_eval (v, F, r[i]);
      ASSERT_ALWAYS(mpz_divisible_p (v, p));
    }
  for (i = 0; i <= d; i++)
    {
      mpz_clear (f[i]);
      mpz_clear (r[i]);
    }
  free (f);
  free (r);
  mpz_clear (p);
  mpz_clear (v);
}

int
main (int argc, const char *argv[])
{
    int d;
    const char* test0[] = {"4294967291", "1", "0", "-3"};
    const char* test1[] = {"18446744073709551557", "1", "2", "3", "5"};
    const char* test2[] = {"18446744073709551629", "1", "-1", "7", "-1"};
    unsigned long iter = 300;

    tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
    tests_common_get_iter (&iter);

    d = argc - 3;

    if (d < 1 && d != -2) {
	fprintf(stderr, "Usage: test-rootfinder p a_d [...] a_1 a_0\n");
	exit(1);
    }

    if (d >= 1)
      test (d, argv[1], argv + 2, -1);
    test (2, test0[0], test0 + 1, 2);
    test (3, test1[0], test1 + 1, 1);
    test (3, test2[0], test2 + 1, 0);

    while (iter--)
      {
        d = 1 + (lrand48 () % 7);
        test (d, "", test0 + 1, -1);
      }

    tests_common_clear ();
    exit (EXIT_SUCCESS);
}

#if 0
// magma code for producing test cases.
s:=1.2;            
p:=10;
while p lt 2^200 do
    for i in [1..100] do
        p:=NextPrime(p);
        d:=Random([2..7]);
        coeffs:=[Random(GF(p)):i in [0..d]];
        F:=PolynomialRing(GF(p))!coeffs;
        printf "in %o", p;
        for c in Reverse(coeffs) do printf " %o", c; end for;
        printf "\n";
        r:=Sort([x[1]: x in Roots(F) ]);
        printf "out";
        for c in r do printf " %o", c; end for;
        printf "\n";
    end for;
    p := Ceiling(p*s);
end while;


#endif

