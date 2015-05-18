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

void prime_product_init(prime_info pi, unsigned long p_max, unsigned long *p_last)
{
  unsigned long p;

  prime_info_init (pi);
  for (p = 2; p <= p_max; p = getprime_mt(pi));
  *p_last = p;
}

void prime_product(mpz_t P, prime_info pi, unsigned long p_max, unsigned long *p_last)
{
  unsigned long p, n = 0, alloc = 0, newalloc, i;
  mpz_t *L = NULL;

  p = *p_last;
  while (p <= p_max)
  {
    if (n >= alloc)
    {
      newalloc = 2 * alloc + 1;
      L = realloc (L, newalloc * sizeof (mpz_t));
      while (alloc < newalloc)
        mpz_init (L[alloc++]);
    }
    mpz_set_ui (L[n++], p);
    p = getprime_mt(pi);
  }
  *p_last = p;
  
  /* FIXME: equilibrate the product */
  while (n > 1)
  {
    for (i = 0; i+1 < n; i+=2)
      mpz_mul (L[i/2], L[i], L[i+1]);
    if (n & 1)
      mpz_swap (L[n/2], L[n-1]);
    n = (n + 1) / 2;
  }
  mpz_set (P, L[0]);
  
  for (i = 0; i < alloc; i++)
    mpz_clear (L[i]);
  free (L);
}

void update_status(mpz_t *R, mpz_t *A, unsigned char *b_status_r,
                   unsigned char *b_status_a, unsigned long int *n, unsigned long int *nb_rel_smooth,
                   unsigned long int rlim, unsigned long int lpbr,
                   unsigned long int *nb_smooth_r, unsigned long int *nb_smooth_a, unsigned long int *nb_useless,
                   unsigned int nb_type[5], unsigned long int *rel_index)
{
  mpz_t z_B3;
  mpz_t z_L2;
  
  unsigned long int tmp;
  unsigned long int i;
  unsigned long int B;  
  unsigned long int L;  


  mpz_init(z_B3);
  mpz_init(z_L2);
    
  memset(nb_type, 0, 5 * sizeof(unsigned int));
  
  mpz_set_ui(z_B3, rlim * rlim);
  mpz_mul_ui(z_B3, z_B3, rlim);
  mpz_set_ui(z_L2, 0);
  mpz_setbit(z_L2, 2 * lpbr);
  B = rlim;
  L = 1UL << lpbr;

  for (i = *nb_rel_smooth; i < *n; i++)
    {
      if (b_status_r[i] == STATUS_UNKNOWN)
      {
        if ( (mpz_cmp(R[i], z_L2) > 0) && (mpz_cmp(R[i], z_B3) < 0) )
        {
          if (b_status_a[i] == STATUS_SMOOTH)
            (*nb_smooth_a)--;
          b_status_r[i] = STATUS_USELESS;
          mpz_set_ui(R[i], 1);
          b_status_a[i] = STATUS_USELESS;
          mpz_set_ui(A[i], 1);
          (*nb_useless)++;
          mpz_swap(R[i], R[(*n)-1]);
          mpz_swap(A[i], A[(*n)-1]);
          tmp = b_status_r[i]; b_status_r[i] = b_status_r[(*n)-1] ; b_status_r[(*n)-1] = tmp;
          tmp = b_status_a[i]; b_status_a[i] = b_status_a[(*n)-1] ; b_status_a[(*n)-1] = tmp;
          tmp = rel_index[i]; rel_index[i] = rel_index[(*n)-1] ; rel_index[(*n)-1] = tmp;
          (*n)--; i--;
          nb_type[0]++;
        }
        else if ( (mpz_cmp_ui(R[i], L) > 0) && (mpz_cmp_ui(R[i], B * B) < 0) )
        {
          if (b_status_a[i] == STATUS_SMOOTH)
            (*nb_smooth_a)--;
          b_status_r[i] = STATUS_USELESS;
          mpz_set_ui(R[i], 1);
          b_status_a[i] = STATUS_USELESS;
          mpz_set_ui(A[i], 1);
          (*nb_useless)++;
          mpz_swap(R[i], R[(*n)-1]);
          mpz_swap(A[i], A[(*n)-1]);
          tmp = b_status_r[i]; b_status_r[i] = b_status_r[(*n)-1] ; b_status_r[(*n)-1] = tmp;
          tmp = b_status_a[i]; b_status_a[i] = b_status_a[(*n)-1] ; b_status_a[(*n)-1] = tmp;
          tmp = rel_index[i]; rel_index[i] = rel_index[(*n)-1] ; rel_index[(*n)-1] = tmp;
          (*n)--; i--;
          nb_type[1]++;
        }
        else if (mpz_cmp_ui(R[i], L) <= 0)  // assume L < B^2
        {
          b_status_r[i] = STATUS_SMOOTH;
          mpz_set_ui(R[i], 1);
          (*nb_smooth_r)++;
          if (b_status_a[i] == STATUS_SMOOTH)
          {
            mpz_swap(R[i], R[*nb_rel_smooth]);
            mpz_swap(A[i], A[*nb_rel_smooth]);
            tmp = b_status_r[i]; b_status_r[i] = b_status_r[*nb_rel_smooth] ; b_status_r[*nb_rel_smooth] = tmp;
            tmp = b_status_a[i]; b_status_a[i] = b_status_a[*nb_rel_smooth] ; b_status_a[*nb_rel_smooth] = tmp;
            tmp = rel_index[i]; rel_index[i] = rel_index[*nb_rel_smooth] ; rel_index[*nb_rel_smooth] = tmp;
            (*nb_rel_smooth)++;
          }
          nb_type[2]++;
        }
        /* now L < B^2 <= R[i] <= L^2 or B^3 <= R[i] */
        else if (mpz_probab_prime_p(R[i], NB_MILLER_RABIN) != 0)
        {
          if (b_status_a[i] == STATUS_SMOOTH)
            (*nb_smooth_a)--;
          b_status_r[i] = STATUS_USELESS;
          mpz_set_ui(R[i], 1);
          b_status_a[i] = STATUS_USELESS;
          mpz_set_ui(A[i], 1);
          (*nb_useless)++;
          mpz_swap(R[i], R[(*n)-1]);
          mpz_swap(A[i], A[(*n)-1]);
          tmp = b_status_r[i]; b_status_r[i] = b_status_r[(*n)-1] ; b_status_r[(*n)-1] = tmp;
          tmp = b_status_a[i]; b_status_a[i] = b_status_a[(*n)-1] ; b_status_a[(*n)-1] = tmp;
          tmp = rel_index[i]; rel_index[i] = rel_index[(*n)-1] ; rel_index[(*n)-1] = tmp;
          (*n)--; i--;
          nb_type[3]++;
        }
        /* now L < B^2 <= R[i] <= L^2 or B^3 <= R[i] and R[i] is composite */
        else if (mpz_cmp_ui(R[i], B * L) <= 0)
        {
          b_status_r[i] = STATUS_SMOOTH;
          mpz_set_ui(R[i], 1);
          (*nb_smooth_r)++;
          if (b_status_a[i] == STATUS_SMOOTH)
          {
            mpz_swap(R[i], R[*nb_rel_smooth]);
            mpz_swap(A[i], A[*nb_rel_smooth]);
            tmp = b_status_r[i]; b_status_r[i] = b_status_r[*nb_rel_smooth] ; b_status_r[*nb_rel_smooth] = tmp;
            tmp = b_status_a[i]; b_status_a[i] = b_status_a[*nb_rel_smooth] ; b_status_a[*nb_rel_smooth] = tmp;
            tmp = rel_index[i]; rel_index[i] = rel_index[*nb_rel_smooth] ; rel_index[*nb_rel_smooth] = tmp;
            (*nb_rel_smooth)++;
          }
          nb_type[4]++;
        }
      }
    }

  mpz_clear(z_B3);
  mpz_clear(z_L2);
}

int main(int argc, char* argv[])
{
  char line[200];
  unsigned long int nb_rel;
  unsigned long int nb_rel_all;
  unsigned long int nb_rel_new;
  unsigned long int nb_rel_read;
  unsigned long int nb_rel_step = 10000000;
  unsigned long int nb_rel_unknown;
  unsigned long int i;
  unsigned long int j;
  FILE *cofac;
  mpz_t *R, *A, P;
  size_t ret;
  double s;
  double t_smooth = 0;
  double t_update = 0;
  double start;
  unsigned int n0_pass;
  unsigned long int lpbr = 33;
  unsigned long int lpba = 33;
  unsigned long int rlim = 250000000; /* UPDATE */
  unsigned long int alim = 500000000; /* UPDATE */
  unsigned long int rlim_step = 750000000; /* UPDATE */
  unsigned long int alim_step = 750000000; /* UPDATE */
  unsigned long int rlim_new;
  unsigned long int alim_new;
  unsigned char *b_status_r;
  unsigned char *b_status_a;
  unsigned long int nb_rel_smooth;
  unsigned long int nb_smooth_r;
  unsigned long int nb_smooth_a;
  unsigned long int nb_useless;
  unsigned long int *rel_index;
  prime_info pi;
  unsigned long int prime;
  unsigned int nb_type[5];
  long int *a;
  unsigned long int *b;

  
  /* Init */
  
  start = seconds();

  fprintf (stderr, "\nInit : %.0f s\n\n", seconds() - start);

  ASSERT_ALWAYS (argc == 2);
  cofac = fopen (argv[1], "r");
  nb_rel_all = 0;
  while (fgets(line, 200, cofac) != NULL)
  {
    nb_rel_all++;
  }
  fseek(cofac, 0, SEEK_SET);

  a = malloc (nb_rel_all * sizeof (long));
  b = malloc (nb_rel_all * sizeof (unsigned long));
  R = malloc (nb_rel_all * sizeof (mpz_t));
  A = malloc (nb_rel_all * sizeof (mpz_t));
  for (i = j = 0; i < nb_rel_all; i++)
  {
    mpz_init (R[i]);
    mpz_init (A[i]);
    ret = gmp_fscanf (cofac, "%ld %lu %Zd %Zd\n", a+j, b+j, R[j], A[j]);
    ASSERT_ALWAYS (ret == 4);
#if 1 /* define to 1 to enable only two large primes on the rational side */
    if (mpz_sizeinbase (R[j], 2) > 66)
    {
    }
    else
#endif
    j++;
  }
  nb_rel_read = j;
  fprintf (stderr, "Read %lu cofactors\n", nb_rel_read);

  mpz_init (P);

  b_status_r = (unsigned char *) malloc(nb_rel_all * sizeof(unsigned char));
  b_status_a = (unsigned char *) malloc(nb_rel_all * sizeof(unsigned char));
  rel_index = (unsigned long int *) malloc(nb_rel_all * sizeof(unsigned long int));
  for (i = 0; i < nb_rel_all; i++)
  {
    b_status_r[i] = STATUS_UNKNOWN;
    b_status_a[i] = STATUS_UNKNOWN;
    rel_index[i] = i;
  }

  nb_rel_smooth = 0;
  nb_rel_unknown = nb_rel_read;

  nb_smooth_r = 0;
  nb_smooth_a = 0;
  nb_useless = 0;

  prime_info_init(pi);

  /* Initial one-side pass (to make rlim and alim be equal) */
  
  if (rlim < alim)
  {
    for (prime = 2; prime <= rlim; prime = getprime_mt(pi));

    while (rlim < alim)
    {
      rlim_new = (rlim + rlim_step < alim ? rlim + rlim_step : alim);

      fprintf (stderr, "\nInitial one-side pass : %0.f s\n\n", seconds() - start);

      fprintf (stderr, "rlim : %lu:%lu\n\n", rlim, rlim_new);

      s = seconds ();
      prime_product (P, pi, rlim_new, &prime);
      fprintf (stderr, "Computed P of %lu bits took %.0f seconds\n",
               mpz_sizeinbase (P, 2), seconds () - s);

      nb_rel = nb_rel_smooth;
      while (nb_rel < nb_rel_unknown)
      {
        nb_rel_new = (nb_rel + nb_rel_step < nb_rel_unknown ? nb_rel + nb_rel_step : nb_rel_unknown);
        s = seconds ();
        t_smooth -= seconds();
        smoothness_test (&(R[nb_rel]), nb_rel_new - nb_rel, P);
        t_smooth += seconds();
        fprintf (stderr, "smoothness_test (%lu cofactors) took %.0f seconds"
                 " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
        nb_rel = nb_rel_new;
      }
      
      rlim = rlim_new;
      t_update -= seconds();
      update_status(R, A, b_status_r, b_status_a, &nb_rel_unknown, &nb_rel_smooth, rlim, lpbr, &nb_smooth_r, &nb_smooth_a, &nb_useless, nb_type, rel_index);
      t_update += seconds();
      fprintf (stderr, "nb_smooth_r = %lu ; nb_useless = %lu ; nb_unknown = %lu ; nb_rel_smooth = %lu ;"
               " %u : %u : %u : %u : %u\n",
               nb_smooth_r, nb_useless, nb_rel_read - nb_smooth_r - nb_useless, nb_rel_smooth,
             nb_type[0], nb_type[1], nb_type[2], nb_type[3], nb_type[4]);
      fprintf (stderr, "t_update: %.0f seconds\n", t_update);
    }
  }
  else if (alim < rlim)
  {
    for (prime = 2; prime <= alim; prime = getprime_mt(pi));

    while (alim < rlim)
    {
      alim_new = (alim + alim_step < rlim ? alim + alim_step : rlim);

      fprintf (stderr, "\nInitial one-side pass : %0.f s\n\n", seconds() - start);

      fprintf (stderr, "\nalim : %lu:%lu\n", alim, alim_new);

      s = seconds ();
      prime_product (P, pi, alim_new, &prime);
      fprintf (stderr, "Computed P of %lu bits took %.0f seconds\n",
               mpz_sizeinbase (P, 2), seconds () - s);

      nb_rel = nb_rel_smooth;
      while (nb_rel < nb_rel_unknown)
      {
        nb_rel_new = (nb_rel + nb_rel_step < nb_rel_unknown ? nb_rel + nb_rel_step : nb_rel_unknown);
        s = seconds ();
        t_smooth -= seconds();
        smoothness_test (&(A[nb_rel]), nb_rel_new - nb_rel, P);
        t_smooth += seconds();
        fprintf (stderr, "smoothness_test (%lu cofactors) took %.0f seconds"
                 " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
        nb_rel = nb_rel_new;
      }
      
      alim = alim_new;
      t_update -= seconds();
      update_status(A, R, b_status_a, b_status_r, &nb_rel_unknown, &nb_rel_smooth, alim, lpba, &nb_smooth_a, &nb_smooth_r, &nb_useless, nb_type, rel_index);
      t_update += seconds();
      fprintf (stderr, "nb_smooth_a = %lu ; nb_useless = %lu ; nb_unknown = %lu ; nb_rel_smooth = %lu ;"
               " %u : %u : %u : %u : %u\n",
               nb_smooth_a, nb_useless, nb_rel_read - nb_smooth_a - nb_useless, nb_rel_smooth,
             nb_type[0], nb_type[1], nb_type[2], nb_type[3], nb_type[4]);
      fprintf (stderr, "t_update: %.0f seconds\n", t_update);
    }
  }
  else  // rlim = alim
  {
    for (prime = 2; prime <= rlim; prime = getprime_mt(pi));
  }

  
  /* Loop */

  n0_pass = 0;
  while ( (rlim < (1UL << lpbr)) || (alim < (1UL << lpba)) )
  {
    n0_pass++;
 
    fprintf (stderr, "\nPass %u : %.0f s\n", n0_pass, seconds() - start);
    
    /* rational side */

    rlim_new = (rlim + rlim_step < (1UL << lpbr) ? rlim + rlim_step : (1UL << lpbr));
    
    fprintf (stderr, "\nrlim : %lu:%lu\n", rlim, rlim_new);

    s = seconds ();
    prime_product (P, pi, rlim_new, &prime);
    fprintf (stderr, "Computed P of %lu bits took %.0f seconds\n",
             mpz_sizeinbase (P, 2), seconds () - s);

    nb_rel = nb_rel_smooth;
    while (nb_rel < nb_rel_unknown)
    {
      nb_rel_new = (nb_rel + nb_rel_step < nb_rel_unknown ? nb_rel + nb_rel_step : nb_rel_unknown);
      s = seconds ();
      t_smooth -= seconds();
      smoothness_test (&(R[nb_rel]), nb_rel_new - nb_rel, P);
      t_smooth += seconds();
      fprintf (stderr, "smoothness_test (%lu cofactors) took %.0f seconds"
               " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
      nb_rel = nb_rel_new;
    }

    rlim = rlim_new;
    t_update -= seconds();
    update_status(R, A, b_status_r, b_status_a, &nb_rel_unknown, &nb_rel_smooth, rlim, lpbr, &nb_smooth_r, &nb_smooth_a, &nb_useless, nb_type, rel_index);
    t_update += seconds();
    fprintf (stderr, "nb_smooth_r = %lu ; nb_useless = %lu ; nb_unknown = %lu ; nb_rel_smooth = %lu ;"
             " %u : %u : %u : %u : %u\n",
             nb_smooth_r, nb_useless, nb_rel_read - nb_smooth_r - nb_useless, nb_rel_smooth,
           nb_type[0], nb_type[1], nb_type[2], nb_type[3], nb_type[4]);
    fprintf (stderr, "t_update: %.0f seconds\n", t_update);

    /* algebraic side */

    alim_new = (alim + alim_step < (1UL << lpba) ? alim + alim_step : (1UL << lpba));

    fprintf (stderr, "\nalim : %lu:%lu\n", alim, alim_new);

    nb_rel = nb_rel_smooth;
    while (nb_rel < nb_rel_unknown)
    {
      nb_rel_new = (nb_rel + nb_rel_step < nb_rel_unknown ? nb_rel + nb_rel_step : nb_rel_unknown);
      s = seconds ();
      t_smooth -= seconds();
      smoothness_test (&(A[nb_rel]), nb_rel_new - nb_rel, P);
      t_smooth += seconds();
      fprintf (stderr, "smoothness_test (%lu cofactors) took %.0f seconds"
               " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
      nb_rel = nb_rel_new;
    }

    alim = alim_new;
    t_update -= seconds();
    update_status(A, R, b_status_a, b_status_r, &nb_rel_unknown, &nb_rel_smooth, alim, lpba, &nb_smooth_a, &nb_smooth_r, &nb_useless, nb_type, rel_index);
    t_update += seconds();
    fprintf (stderr, "nb_smooth_a = %lu ; nb_useless = %lu ; nb_unknown = %lu ; nb_rel_smooth = %lu ;"
             " %u : %u : %u : %u : %u\n",
             nb_smooth_a, nb_useless, nb_rel_read - nb_smooth_a - nb_useless, nb_rel_smooth,
             nb_type[0], nb_type[1], nb_type[2], nb_type[3], nb_type[4]);
    fprintf (stderr, "t_update: %.0f seconds\n", t_update);
  }

  fprintf (stderr, "\nCollect smooth relations: %.0f s\n\n", seconds() - start);

  for (i = 0; i < nb_rel_smooth; i++)
  {
    ASSERT_ALWAYS ( (b_status_r[i] == STATUS_SMOOTH) && (b_status_a[i] == STATUS_SMOOTH) );
    printf ("Smooth: index=%lu a=%ld b=%lu\n", rel_index[i], a[rel_index[i]], b[rel_index[i]]);
  }
  fprintf (stderr, "\nFound %u smooth relations\n", nb_rel_smooth);

  mpz_clear (P);

  for (i = 0; i < nb_rel_all; i++)
    {
      mpz_clear (R[i]);
      mpz_clear (A[i]);
    }
    
  free (b_status_r);
  free (b_status_a);
  free (R);
  free (A);
  fclose (cofac);

  fprintf (stderr, "\nThe end : %.0f s\n", seconds() - start);
}

