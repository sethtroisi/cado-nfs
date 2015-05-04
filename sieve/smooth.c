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

#define NB_MILLER_RABIN  2

#define STATUS_SMOOTH  0
#define STATUS_UNKNOWN 1
#define STATUS_USELESS 2
#define STATUS_NEW 3

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
   R[0], ..., R[n-1] are cofactors
   P is the product of primes
   Output:
   Each R[j] has been divided by its P-smooth part
*/
void
smoothness_test (mpz_t *R, unsigned long n, mpz_t P, unsigned char *status)
{
  unsigned long h = tree_height (n), i, j, w[64];
  mpz_t **T;
  mpz_t w1, w2;

  mpz_init(w1);
  mpz_init(w2);

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
  fprintf (stderr, "T[h][0] has %lu bits\n", mpz_sizeinbase (T[h][0], 2));

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
      mpz_mod (w1, T[1][j], T[0][2*j]);
      mpz_gcd(w2, w1, T[0][2*j]);
      mpz_divexact(T[0][2*j], T[0][2*j], w2);
      if (mpz_cmp_ui (w2, 1) > 0)
        status[2*j] = STATUS_NEW;
      mpz_mod (w1, T[1][j], T[0][2*j+1]);
      mpz_gcd(w2, w1, T[0][2*j+1]);
      mpz_divexact(T[0][2*j+1], T[0][2*j+1], w2);
      if (mpz_cmp_ui (w2, 1) > 0)
        status[2*j+1] = STATUS_NEW;
    }
  if (n & 1)
    {
      mpz_gcd(w2, T[1][w[1]-1], T[0][n-1]);
      mpz_divexact(T[0][n-1], T[0][n-1], w2);
      if (mpz_cmp_ui (w2, 1) > 0)
        status[n-1] = STATUS_NEW;
    }

  mpz_clear(w1);
  mpz_clear(w2);
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
  prime_info pi;

  prime_info_init (pi);
  for (p = 2; p <= pmax; p = getprime_mt (pi))
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
  prime_info_clear (pi);

    /* FIXME: equilibrate the product */
  while (n > 1)
    {
      for (i = 0; i+1 < n; i+=2)
        mpz_mul (L[i/2], L[i], L[i+1]);
      if (n & 1)
        mpz_set(L[n/2], L[n-1]);
//        mpz_swap (L[n/2], L[n-1]);
      n = (n + 1) / 2;
    }
//  mpz_swap (P, L[0]);
  mpz_set (P, L[0]);

  for (i = 0; i < alloc; i++)
    mpz_clear (L[i]);
  free (L);
}

void update_status(mpz_t *R, mpz_t *A, unsigned char *b_status_r, unsigned char *b_status_a, unsigned long int n, unsigned int b_side_R,
                   unsigned long int rlim, unsigned long int lpbr, unsigned long int alim, unsigned long int lpba,
                   unsigned long int *nb_smooth_r, unsigned long int *nb_smooth_a, unsigned long int *nb_useless)
{
  mpz_t z_B3;
  mpz_t z_L2;

  unsigned long int i, L;


  mpz_init(z_B3);
  mpz_init(z_L2);

  if (b_side_R)
  {
    mpz_set_ui(z_B3, rlim * rlim);
    mpz_mul_ui(z_B3, z_B3, rlim);
    mpz_set_ui(z_L2, 0);
    mpz_setbit(z_L2, 2 * lpbr);
    L = 1UL << lpbr;

    for (i = 0; i < n; i++)
    {
      if (b_status_r[i] == STATUS_NEW)
      {
        if (mpz_cmp_ui (R[i], 1) == 0)
        {
          b_status_r[i] = STATUS_SMOOTH;
          (*nb_smooth_r)++;
        }
        else if (mpz_cmp_ui(R[i], L) <= 0)  // works only if rlim*rlim > 2^lpbr (which is assumed)
        {
          b_status_r[i] = STATUS_SMOOTH;
          mpz_set_ui(R[i], 1);
          (*nb_smooth_r)++;
        }
        else if ( (mpz_cmp_ui(R[i], L) > 0) && (mpz_cmp_ui(R[i], rlim * rlim) < 0) )
        {
          if (b_status_a[i] == STATUS_SMOOTH)
            (*nb_smooth_a)--;
          b_status_r[i] = STATUS_USELESS;
          mpz_set_ui(R[i], 1);
          b_status_a[i] = STATUS_USELESS;
          mpz_set_ui(A[i], 1);
          (*nb_useless)++;
        }
        else if ( (mpz_cmp(R[i], z_L2) > 0) && (mpz_cmp(R[i], z_B3) < 0) )
        {
          if (b_status_a[i] == STATUS_SMOOTH)
            (*nb_smooth_a)--;
          b_status_r[i] = STATUS_USELESS;
          mpz_set_ui(R[i], 1);
          b_status_a[i] = STATUS_USELESS;
          mpz_set_ui(A[i], 1);
          (*nb_useless)++;
        }
        else if ((mpz_cmp_ui(R[i], rlim * L) <= 0) && (mpz_probab_prime_p(R[i], NB_MILLER_RABIN) == 0) )
        {
          b_status_r[i] = STATUS_SMOOTH;
          mpz_set_ui(R[i], 1);
          (*nb_smooth_r)++;
        }
        else if ((mpz_cmp_ui(R[i], L) > 0) && (mpz_probab_prime_p(R[i], NB_MILLER_RABIN) != 0))
        {
          if (b_status_a[i] == STATUS_SMOOTH)
            (*nb_smooth_a)--;
          b_status_r[i] = STATUS_USELESS;
          mpz_set_ui(R[i], 1);
          b_status_a[i] = STATUS_USELESS;
          mpz_set_ui(A[i], 1);
          (*nb_useless)++;
        }
        else
          b_status_r[i] = STATUS_UNKNOWN;
      }
    }
  }
  else
  {
    mpz_set_ui(z_B3, alim * alim);
    mpz_mul_ui(z_B3, z_B3, alim);
    mpz_set_ui(z_L2, 0);
    mpz_setbit(z_L2, 2 * lpba);
    L = 1UL << lpba;

    for (i = 0; i < n; i++)
    {
      if (b_status_a[i] == STATUS_NEW)
      {
        if (mpz_cmp_ui (A[i], 1) == 0)
        {
          b_status_a[i] = STATUS_SMOOTH;
          (*nb_smooth_a)++;
        }
        else if (mpz_cmp_ui(A[i], L) <= 0)  // works only if alim*alim > 2^lpba (which is assumed)
        {
          b_status_a[i] = STATUS_SMOOTH;
          mpz_set_ui(A[i], 1);
          (*nb_smooth_a)++;
        }
        else if ( (mpz_cmp_ui(A[i], L) > 0) && (mpz_cmp_ui(A[i], alim * alim) < 0) )
        {
          if (b_status_r[i] == STATUS_SMOOTH)
            (*nb_smooth_r)--;
          b_status_a[i] = STATUS_USELESS;
          mpz_set_ui(A[i], 1);
          b_status_r[i] = STATUS_USELESS;
          mpz_set_ui(R[i], 1);
          (*nb_useless)++;
        }
        else if ( (mpz_cmp(A[i], z_L2) > 0) && (mpz_cmp(A[i], z_B3) < 0) )
        {
          if (b_status_r[i] == STATUS_SMOOTH)
            (*nb_smooth_r)--;
          b_status_a[i] = STATUS_USELESS;
          mpz_set_ui(A[i], 1);
          b_status_r[i] = STATUS_USELESS;
          mpz_set_ui(R[i], 1);
          (*nb_useless)++;
        }
        else if ((mpz_cmp_ui(A[i], alim * L) <= 0) && (mpz_probab_prime_p(A[i], NB_MILLER_RABIN) == 0) )
        {
          b_status_a[i] = STATUS_SMOOTH;
          mpz_set_ui(A[i], 1);
          (*nb_smooth_a)++;
        }
        else if ((mpz_cmp_ui(A[i], L) > 0) && (mpz_probab_prime_p(A[i], NB_MILLER_RABIN) != 0))
        {
          if (b_status_r[i] == STATUS_SMOOTH)
            (*nb_smooth_r)--;
          b_status_a[i] = STATUS_USELESS;
          mpz_set_ui(A[i], 1);
          b_status_r[i] = STATUS_USELESS;
          mpz_set_ui(R[i], 1);
          (*nb_useless)++;
        }
        else
          b_status_a[i] = STATUS_UNKNOWN;
      }
    }
  }

  mpz_clear(z_B3);
  mpz_clear(z_L2);
}

int
main (int argc, char *argv[])
{
  FILE *cofac;
  mpz_t *R, *A, P;
  unsigned long n = 5567917; /* number of cofactors, UPDATE */
  unsigned long n_init;
  unsigned long i;
  size_t ret;
  double s, t_smooth = 0, t_update = 0;
  unsigned int n0_pass;
  unsigned long int lpbr = 33;
  unsigned long int lpba = 33;
  unsigned long int rlim_low;
  unsigned long int alim_low;
  unsigned long int rlim_high;
  unsigned long int alim_high;
  unsigned long int rlim_step = 500000000;
  unsigned long int alim_step = 500000000;
  unsigned char *b_status_r;
  unsigned char *b_status_a;
  unsigned long int nb_smooth_r;
  unsigned long int nb_smooth_a;
  unsigned long int nb_useless;
  long *a;
  unsigned long *b;

  n_init = n;

  b_status_r = (unsigned char *) malloc(n * sizeof(unsigned char));
  b_status_a = (unsigned char *) malloc(n * sizeof(unsigned char));
  for (i = 0; i < n; i++)
  {
    b_status_r[i] = STATUS_NEW;
    b_status_a[i] = STATUS_NEW;
  }
  ASSERT_ALWAYS (argc == 2);
  cofac = fopen (argv[1], "r");
  a = malloc (n * sizeof (long));
  b = malloc (n * sizeof (unsigned long));
  R = malloc (n * sizeof (mpz_t));
  A = malloc (n * sizeof (mpz_t));
  for (i = 0; i < n; i++)
    {
      mpz_init (R[i]);
      mpz_init (A[i]);
      ret = gmp_fscanf (cofac, "%ld %lu %Zd %Zd\n", a+i, b+i, R[i], A[i]);
      ASSERT_ALWAYS (ret == 4);
    }
  fprintf (stderr, "Read %lu cofactors\n", n);
  mpz_init (P);


  /**** Initial pass ***/

  fprintf (stderr, "\n");
  fprintf (stderr, "Initial pass\n");

  rlim_high = 250000000;
  alim_high = 500000000;

  nb_smooth_r = 0;
  nb_smooth_a = 0;
  nb_useless = 0;

  t_update -= seconds ();
  update_status(R, A, b_status_r, b_status_a, n, 1, rlim_high, lpbr, alim_high, lpba, &nb_smooth_r, &nb_smooth_a, &nb_useless);
  update_status(R, A, b_status_r, b_status_a, n, 0, rlim_high, lpbr, alim_high, lpba, &nb_smooth_r, &nb_smooth_a, &nb_useless);
  t_update += seconds ();
  fprintf (stderr, "Update time %.0fs\n", t_update);

#define NPASS 20

  rlim_step = 1 + ((1UL << lpbr) - rlim_high - 1) / NPASS;
  alim_step = 1 + ((1UL << lpba) - alim_high - 1) / NPASS;

  n0_pass = 0;
  while (n0_pass < NPASS)
  {
    n0_pass++;

    fprintf (stderr, "\n");
    fprintf (stderr, "Pass %u\n", n0_pass);

    /* rational side */

    rlim_low = rlim_high;
    rlim_high += rlim_step;

    fprintf (stderr, "\n");
    fprintf (stderr, "rlim : %lu:%lu\n", rlim_low, rlim_high);
    fprintf (stderr, "nb_smooth_r = %lu ; nb_useless = %lu ; nb_unknown = %lu\n", nb_smooth_r, nb_useless, n_init - nb_smooth_r - nb_useless);

    s = seconds ();
    prime_product (P, rlim_low, rlim_high); /* UPDATE */
    fprintf (stderr, "Computed P of %lu bits took %.0f seconds\n",
             mpz_sizeinbase (P, 2), seconds () - s);
    s = seconds ();
    t_smooth -= seconds ();
    smoothness_test (R, n, P, b_status_r);
    t_smooth += seconds ();
    fprintf (stderr, "smoothness_test took %.0fs (total %.0fs so far)\n",
             seconds () - s, t_smooth);

    t_update -= seconds ();
    update_status(R, A, b_status_r, b_status_a, n, 1, rlim_high, lpbr, alim_high, lpba, &nb_smooth_r, &nb_smooth_a, &nb_useless);
    t_update += seconds ();
    fprintf (stderr, "nb_smooth_r = %lu ; nb_useless = %lu ; nb_unknown = %lu\n", nb_smooth_r, nb_useless, n_init - nb_smooth_r - nb_useless);
    fprintf (stderr, "Update time %.0fs so far\n", t_update);

    /* algebraic side */

    alim_low = alim_high;
    alim_high += alim_step;

    fprintf (stderr, "\n");
    fprintf (stderr, "alim : %lu:%lu\n", alim_low, alim_high);
    fprintf (stderr, "nb_smooth_a = %lu ; nb_useless = %lu ; nb_unknown = %lu\n", nb_smooth_a, nb_useless, n_init - nb_smooth_a - nb_useless);

    s = seconds ();
    prime_product (P, alim_low, alim_high); /* UPDATE */
    fprintf (stderr, "Computed P of %lu bits took %.0f seconds\n",
             mpz_sizeinbase (P, 2), seconds () - s);
    s = seconds ();
    t_smooth -= seconds ();
    smoothness_test (A, n, P, b_status_a);
    t_smooth += seconds ();
    fprintf (stderr, "smoothness_test took %.0fs (total %.0fs so far)\n",
             seconds () - s, t_smooth);

    t_update -= seconds ();
    update_status(R, A, b_status_r, b_status_a, n, 0, rlim_high, lpbr, alim_high, lpba, &nb_smooth_r, &nb_smooth_a, &nb_useless);
    t_update += seconds ();
    fprintf (stderr, "nb_smooth_a = %lu ; nb_useless = %lu ; nb_unknown = %lu\n", nb_smooth_a, nb_useless, n_init - nb_smooth_a - nb_useless);
    fprintf (stderr, "Update time %.0fs so far\n", t_update);
  }

  /* those relations that remain unknown are those for which we were not able
     to determine the status on each side, in any case they are non smooth */

  /* print smooth cofactors */
  unsigned long nb_smooth = 0;
  for (i = 0; i < n; i++)
    {
      if (b_status_r[i] == STATUS_SMOOTH && b_status_a[i] == STATUS_SMOOTH)
        {
          printf ("Smooth: a=%ld b=%lu\n", a[i], b[i]);
          nb_smooth ++;
        }
    }
  fprintf (stderr, "Found %lu relations\n", nb_smooth);

  mpz_clear (P);

  for (i = 0; i < n; i++)
    {
//      gmp_fprintf (stderr, "%Zd %Zd\n", R[i], A[i]);
      mpz_clear (R[i]);
      mpz_clear (A[i]);
    }

  free (b_status_r);
  free (b_status_a);
  free (R);
  free (A);
  free (a);
  free (b);
  fclose (cofac);
}

