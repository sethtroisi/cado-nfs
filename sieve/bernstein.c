/* Bernstein's smoothness test */

/* This is a first very naive implementation. To use it:
   1) create two files /tmp/cofac_r and /tmp/cofac_a of cofactors on the
      rational and algebraic sides respectively (should be of the same size)
   2) replace n = 479233 by the common number of lines of those files
   3) update the ranges for rational and algebraic primes on lines UPDATE
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "utils.h"

unsigned long
tree_height (unsigned long n)
{
  unsigned long h = 0;

  while (n > 1)
    {
      h ++;
      n = (n + 1) / 2;
    }
  return h;
}

/* Input:
   R[0], ..., R[n-1] are cofactors (will be destroyed)
   P is the product of primes
   Output:
   for all cofactors that have only factors in P, the new value of R[j] is 0
*/
void
smoothness_test (mpz_t *R, unsigned long n, mpz_t P)
{
  unsigned long h = tree_height (n), i, j, w[64];
  mpz_t **T;

  ASSERT_ALWAYS(h >= 2);

  T = (mpz_t**) malloc ((h + 1) * sizeof (mpz_t*));
  T[0] = R;

  w[0] = n;
  /* initialize tree */
  for (i = 1; i <= h; i++)
    {
      w[i] = 1 + ((n - 1) >> i);
      T[i] = (mpz_t*) malloc (w[i] * sizeof (mpz_t));
      for (j = 0; j < w[i]; j++)
        mpz_init (T[i][j]);
    }

  /* compute product tree */
  for (i = 1; i <= h; i++)
    {
      for (j = 0; j < w[i-1] / 2; j++)
        mpz_mul (T[i][j], T[i-1][2*j], T[i-1][2*j+1]);
      if (w[i-1] & 1)
        mpz_set (T[i][w[i]-1], T[i-1][w[i-1]-1]);
    }
  printf ("T[h][0] has %lu bits\n", mpz_sizeinbase (T[h][0], 2));

  /* compute remainder tree */
  mpz_mod (T[h][0], P, T[h][0]);
  for (i = h; i > 1; i--)
    {
      for (j = 0; j < w[i-1] / 2; j++)
        {
          mpz_mod (T[i-1][2*j], T[i][j], T[i-1][2*j]);
          mpz_mod (T[i-1][2*j+1], T[i][j], T[i-1][2*j+1]);
        }
      if (w[i-1] & 1)
        mpz_swap (T[i-1][w[i-1]-1], T[i][w[i]-1]);
    }

  /* special last loop for i=1 */
  for (j = 0; j < n / 2; j++)
    {
      mpz_mod (T[2][0], T[1][j], T[0][2*j]);
      if (mpz_cmp_ui (T[2][0], 0) != 0) /* put non-smooth residues to 0 */
        mpz_set_ui (T[0][2*j], 0);
      mpz_mod (T[2][0], T[1][j], T[0][2*j+1]);
      if (mpz_cmp_ui (T[2][0], 0) != 0) /* put non-smooth residues to 0 */
        mpz_set_ui (T[0][2*j+1], 0);
    }
  if (n & 1)
    {
      if (mpz_cmp_ui (T[1][w[1]-1], 0) != 0) /* put non-smooth residues to 0 */
        mpz_set_ui (T[0][n-1], 0);
    }

  for (i = 1; i <= h; i++)
    {
      for (j = 0; j < w[i]; j++)
        mpz_clear (T[i][j]);
      free (T[i]);
    }
  free (T);
}

void
prime_product (mpz_t P, unsigned long lim, unsigned long pmax)
{
  unsigned long p, n = 0, alloc = 0, newalloc, i;
  mpz_t *L = NULL;

  for (p = 2; p <= pmax; p = getprime (p))
    {
      if (p > lim)
        {
          if (n >= alloc)
            {
              newalloc = 2 * alloc + 1;
              L = realloc (L, newalloc * sizeof (mpz_t));
              while (alloc < newalloc)
                mpz_init (L[alloc++]);
            }
          mpz_set_ui (L[n++], p);
        }
    }
  /* FIXME: equilibrate the product */
  while (n > 1)
    {
      for (i = 0; i+1 < n; i+=2)
        mpz_mul (L[i/2], L[i], L[i+1]);
      if (n & 1)
        mpz_swap (L[n/2], L[n-1]);
      n = (n + 1) / 2;
    }
  mpz_swap (P, L[0]);
  for (i = 0; i < alloc; i++)
    mpz_clear (L[i]);
  free (L);
}

int
main ()
{
  FILE *cofac_r, *cofac_a;
  mpz_t *R, *A, P;
  unsigned long n = 479233; /* number of cofactors, UPDATE */
  unsigned long i;
  size_t ret;
  double s;

  cofac_r = fopen ("/tmp/cofac_r", "r");
  cofac_a = fopen ("/tmp/cofac_a", "r");
  R = malloc (n * sizeof (mpz_t));
  A = malloc (n * sizeof (mpz_t));
  for (i = 0; i < n; i++)
    {
      mpz_init (R[i]);
      mpz_init (A[i]);
      ret = mpz_inp_str (R[i], cofac_r, 10);
      ASSERT_ALWAYS (ret > 0);
      ret = mpz_inp_str (A[i], cofac_a, 10);
      ASSERT_ALWAYS (ret > 0);
    }
  fprintf (stderr, "Read %lu cofactors\n", n);
  mpz_init (P);

  /* rational side */
  s = seconds ();
  prime_product (P, 250000000, 280000000); /* UPDATE */
  fprintf (stderr, "Computed P of %lu bits took %.0f seconds\n",
           mpz_sizeinbase (P, 2), seconds () - s);
  s = seconds ();
  smoothness_test (R, n, P);
  printf ("smoothness_test took %.0f seconds\n", seconds () - s);

  for (i = 0; i < n; i++)
    if (mpz_cmp_ui (R[i], 0) == 0)
      mpz_set_ui (A[i], 1);

  /* algebraic side */
  s = seconds ();
  prime_product (P, 500000000, 530000000); /* UPDATE */
  fprintf (stderr, "Computed P of %lu bits took %.0f seconds\n",
           mpz_sizeinbase (P, 2), seconds () - s);
  s = seconds ();
  smoothness_test (A, n, P);
  printf ("smoothness_test took %.0f seconds\n", seconds () - s);

  mpz_clear (P);

  for (i = 0; i < n; i++)
    {
      if (mpz_cmp_ui (R[i], 0) != 0 && mpz_cmp_ui (A[i], 0) != 0)
        gmp_printf ("%Zd %Zd\n", R[i], A[i]);
      mpz_clear (R[i]);
      mpz_clear (A[i]);
    }
  free (R);
  free (A);
  fclose (cofac_r);
  fclose (cofac_a);
}
