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
smoothness_test (mpz_t *R, unsigned long n, mpz_t P)
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
      mpz_mod (w1, T[1][j], T[0][2*j]);
      mpz_gcd(w2, w1, T[0][2*j]);
      mpz_divexact(T[0][2*j], T[0][2*j], w2);
      mpz_mod (w1, T[1][j], T[0][2*j+1]);
      mpz_gcd(w2, w1, T[0][2*j+1]);
      mpz_divexact(T[0][2*j+1], T[0][2*j+1], w2);
    }
  if (n & 1)
    {
      mpz_gcd(w2, T[1][w[1]-1], T[0][n-1]);
      mpz_divexact(T[0][n-1], T[0][n-1], w2);
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

  unsigned long int i;


  mpz_init(z_B3);
  mpz_init(z_L2);

  if (b_side_R)
  {
    mpz_set_ui(z_B3, rlim * rlim);
    mpz_mul_ui(z_B3, z_B3, rlim);
    mpz_set_ui(z_L2, 0);
    mpz_setbit(z_L2, 2 * lpbr);

    for (i = 0; i < n; i++)
    {
      if (b_status_r[i] == STATUS_UNKNOWN)
      {
        if (mpz_cmp_ui (R[i], 1) == 0)
        {
          b_status_r[i] = STATUS_SMOOTH;
          (*nb_smooth_r)++;
        }
        else if (mpz_cmp_ui(R[i], (1UL << lpbr)) <= 0)  // works only if rlim*rlim > 2^lpbr (which is assumed)
        {
          b_status_r[i] = STATUS_SMOOTH;
          mpz_set_ui(R[i], 1);
          (*nb_smooth_r)++;
        }
        else if ( (mpz_cmp_ui(R[i], (1UL << lpbr)) > 0) && (mpz_cmp_ui(R[i], rlim * rlim) < 0) )
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
        else if ( (mpz_cmp_ui(R[i], rlim * (1UL << lpbr)) <= 0) && (mpz_probab_prime_p(R[i], NB_MILLER_RABIN) == 0) )
        {
          b_status_r[i] = STATUS_SMOOTH;
          mpz_set_ui(R[i], 1);
          (*nb_smooth_r)++;
        }
        else if ( (mpz_cmp_ui(R[i], (1UL << lpbr)) > 0) && (mpz_probab_prime_p(R[i], NB_MILLER_RABIN) != 0) )
        {
          if (b_status_a[i] == STATUS_SMOOTH)
            (*nb_smooth_a)--;
          b_status_r[i] = STATUS_USELESS;
          mpz_set_ui(R[i], 1);
          b_status_a[i] = STATUS_USELESS;
          mpz_set_ui(A[i], 1);
          (*nb_useless)++;
        }
      }
    }
  }
  else
  {
    mpz_set_ui(z_B3, alim * alim);
    mpz_mul_ui(z_B3, z_B3, alim);
    mpz_set_ui(z_L2, 0);
    mpz_setbit(z_L2, 2 * lpba);

    for (i = 0; i < n; i++)
    {
      if (b_status_a[i] == STATUS_UNKNOWN)
      {
        if (mpz_cmp_ui (A[i], 1) == 0)
        {
          b_status_a[i] = STATUS_SMOOTH;
          (*nb_smooth_a)++;
        }
        else if (mpz_cmp_ui(A[i], (1UL << lpba)) <= 0)  // works only if alim*alim > 2^lpba (which is assumed)
        {
          b_status_a[i] = STATUS_SMOOTH;
          mpz_set_ui(A[i], 1);
          (*nb_smooth_a)++;
        }
        else if ( (mpz_cmp_ui(A[i], (1UL << lpba)) > 0) && (mpz_cmp_ui(A[i], alim * alim) < 0) )
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
        else if ( (mpz_cmp_ui(A[i], alim * (1UL << lpba)) <= 0) && (mpz_probab_prime_p(A[i], NB_MILLER_RABIN) == 0) )
        {
          b_status_a[i] = STATUS_SMOOTH;
          mpz_set_ui(A[i], 1);
          (*nb_smooth_a)++;
        }
        else if ( (mpz_cmp_ui(A[i], (1UL << lpba)) > 0) && (mpz_probab_prime_p(A[i], NB_MILLER_RABIN) != 0) )
        {
          if (b_status_r[i] == STATUS_SMOOTH)
            (*nb_smooth_r)--;
          b_status_a[i] = STATUS_USELESS;
          mpz_set_ui(A[i], 1);
          b_status_r[i] = STATUS_USELESS;
          mpz_set_ui(R[i], 1);
          (*nb_useless)++;
        }
      }
    }
  }

  mpz_clear(z_B3);
  mpz_clear(z_L2);
}

int
main ()
{
  FILE *cofac_r, *cofac_a;
  mpz_t *R, *A, P;
  unsigned long n = 479233; /* number of cofactors, UPDATE */
  unsigned long n_init;
  unsigned long i;
  size_t ret;
  double s;
  unsigned int n0_pass;
  unsigned long int lpbr = 33;
  unsigned long int lpba = 33;
  unsigned long int rlim_low;
  unsigned long int alim_low;
  unsigned long int rlim_high;
  unsigned long int alim_high;
  unsigned long int rlim_step = 50000000;
  unsigned long int alim_step = 50000000;
  unsigned char *b_status_r;
  unsigned char *b_status_a;
  unsigned long int nb_smooth_r;
  unsigned long int nb_smooth_a;
  unsigned long int nb_useless;


  n_init = n;

  b_status_r = (unsigned char *) malloc(n * sizeof(unsigned char));
  b_status_a = (unsigned char *) malloc(n * sizeof(unsigned char));
  for (i = 0; i < n; i++)
  {
    b_status_r[i] = STATUS_UNKNOWN;
    b_status_a[i] = STATUS_UNKNOWN;
  }
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


  /**** Initial pass ***/

  printf("\n");
  printf("Initial pass\n");

  rlim_high = 250000000;
  alim_high = 500000000;

  nb_smooth_r = 0;
  nb_smooth_a = 0;
  nb_useless = 0;

  update_status(R, A, b_status_r, b_status_a, n, 1, rlim_high, lpbr, alim_high, lpba, &nb_smooth_r, &nb_smooth_a, &nb_useless);
  update_status(R, A, b_status_r, b_status_a, n, 0, rlim_high, lpbr, alim_high, lpba, &nb_smooth_r, &nb_smooth_a, &nb_useless);

  n0_pass = 0;
  while (n0_pass < 20)
  {
    n0_pass++;

    printf("\n");
    printf("Pass %u\n", n0_pass);

    /* rational side */

    rlim_low = rlim_high;
    rlim_high += rlim_step;

    printf("\n");
    printf("rlim : %lu:%lu\n", rlim_low, rlim_high);
    printf("nb_smooth_r = %lu ; nb_useless = %lu ; nb_unknown = %lu\n", nb_smooth_r, nb_useless, n_init - nb_smooth_r - nb_useless);

    s = seconds ();
    prime_product (P, rlim_low, rlim_high); /* UPDATE */
    fprintf (stderr, "Computed P of %lu bits took %.0f seconds\n",
             mpz_sizeinbase (P, 2), seconds () - s);
    s = seconds ();
    smoothness_test (R, n, P);
    printf ("smoothness_test took %.0f seconds\n", seconds () - s);

    update_status(R, A, b_status_r, b_status_a, n, 1, rlim_high, lpbr, alim_high, lpba, &nb_smooth_r, &nb_smooth_a, &nb_useless);
    printf("nb_smooth_r = %lu ; nb_useless = %lu ; nb_unknown = %lu\n", nb_smooth_r, nb_useless, n_init - nb_smooth_r - nb_useless);

    /* algebraic side */

    alim_low = alim_high;
    alim_high += alim_step;

    printf("\n");
    printf("alim : %lu:%lu\n", alim_low, alim_high);
    printf("nb_smooth_a = %lu ; nb_useless = %lu ; nb_unknown = %lu\n", nb_smooth_a, nb_useless, n_init - nb_smooth_a - nb_useless);

    s = seconds ();
    prime_product (P, alim_low, alim_high); /* UPDATE */
    fprintf (stderr, "Computed P of %lu bits took %.0f seconds\n",
             mpz_sizeinbase (P, 2), seconds () - s);
    s = seconds ();
    smoothness_test (A, n, P);
    printf ("smoothness_test took %.0f seconds\n", seconds () - s);

    update_status(R, A, b_status_r, b_status_a, n, 0, rlim_high, lpbr, alim_high, lpba, &nb_smooth_r, &nb_smooth_a, &nb_useless);
    printf("nb_smooth_a = %lu ; nb_useless = %lu ; nb_unknown = %lu\n", nb_smooth_a, nb_useless, n_init - nb_smooth_a - nb_useless);
  }

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
  fclose (cofac_r);
  fclose (cofac_a);
}

