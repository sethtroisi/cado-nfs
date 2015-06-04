/* Bernstein's smoothness test */

/* This is a first very naive implementation. To use it:
   1) create a file (say cofac) containing lines of the following form:
      a b cofac_r cofac_a
      where cofac_r and cofac_a and cofactors on the rational and algebraic
      sides respectively
   2) run ./smooth cofac after updating the lines with "UPDATE".
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <gmp.h>
#include "smooth.h"
#include "utils.h"
#include "portability.h"

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

/* return the product tree formed from R[0..n-1].
   Put in w[i] the number of elements of level i:
   w[0] = n, w[1] = ceil(n/2), ... */
mpz_t**
product_tree (mpz_t *R, unsigned long n, unsigned long *w)
{
  unsigned long h = tree_height (n), i, j;
  mpz_t **T;

  T = (mpz_t**) malloc ((h + 1) * sizeof (mpz_t*));
  T[0] = R;

  /* initialize tree */
  w[0] = n;
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

  return T;
}

/* Clear the product tree (except T[0] which is assumed to have been
   allocated differently). */
void
clear_product_tree (mpz_t **T, unsigned long n, unsigned long *w)
{
  unsigned long i, j, h = tree_height (n);

  for (i = 1; i <= h; i++)
    {
      for (j = 0; j < w[i]; j++)
        mpz_clear (T[i][j]);
      free (T[i]);
    }
  free (T);
}

/* compute the remainder of P modulo the product tree T, up to T[1]
   (last step from T[1] to T[0] is special) */
void
remainder_tree (mpz_t **T, unsigned long n, unsigned long *w, mpz_t P)
{
  unsigned long h = tree_height (n), i, j;

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
}

#define MAX_DEPTH 32

/* Input:
   R[0], ..., R[n-1] are cofactors
   P is the product of primes
   Output:
   Each R[j] has been divided by its P-smooth part
*/
void
smoothness_test (mpz_t *R, unsigned long n, mpz_t P)
{
  unsigned long h = tree_height (n), j, w[MAX_DEPTH];
  mpz_t **T;
  mpz_t w1, w2;

  mpz_init (w1);
  mpz_init (w2);

  T = product_tree (R, n, w);

  fprintf (stderr, "T[h][0] has %lu bits\n", mpz_sizeinbase (T[h][0], 2));

  /* compute remainder tree */
  remainder_tree (T, n, w, P);

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

  mpz_clear (w1);
  mpz_clear (w2);
  clear_product_tree (T, n, w);
}

void
prime_product_init (prime_info pi, unsigned long p_max, unsigned long *p_last)
{
  unsigned long p;

  prime_info_init (pi);
  for (p = 2; p <= p_max; p = getprime_mt (pi));
  *p_last = p;
}

/* return a list L of all primes < pmax,
   and put in *n0 the number of such primes */
mpz_t*
prime_list (unsigned long pmax, unsigned long *n0)
{
  prime_info pi;
  mpz_t *L = NULL;
  unsigned long p, n, alloc, newalloc;

  prime_info_init (pi);
  for (p = 2, n = alloc = 0; p < pmax; p = getprime_mt (pi), n++)
    {
      if (n >= alloc)
        {
          newalloc = 3 * (alloc / 2) + 2;
          L = realloc (L, newalloc * sizeof (mpz_t));
          alloc = newalloc;
        }
      mpz_init (L[n]); /* initialize only when needed */
      mpz_set_ui (L[n], p);
    }
  L = realloc (L, n * sizeof (mpz_t));
  prime_info_clear (pi);
  *n0 = n;
  return L;
}

/* FIXME: since we don't need to keep the indidivual primes here,
   instead of allocating spaces for n mpz_t data structures, we need
   only to allocate O(log n) [see function compute_biproduct below] */
void
prime_product (mpz_t P, prime_info pi, unsigned long p_max,
               unsigned long *p_last)
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

/* invariant:
   relations 0 to *nb_rel_smooth-1 are smooth
   relations *nb_rel_smooth to *n-1 are unknown
   relations >= *n are non-smooth (useless)
*/
void
update_status (mpz_t *R, mpz_t *A,
               unsigned char *b_status_r, unsigned char *b_status_a,
               unsigned long *n, unsigned long *nb_rel_smooth,
               unsigned long rlim, unsigned long lpbr,
               unsigned long *nb_smooth_r, unsigned long *nb_smooth_a,
               unsigned long *nb_useless,
               unsigned int nb_type[5], int64_t *a, uint64_t *b)
{
  mpz_t z_B3;
  mpz_t z_L2;

  unsigned long tmp;
  unsigned long i;
  unsigned long B;
  unsigned long L;

  int64_t atmp;
  uint64_t btmp;

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
          /* relation i is useless, swap it with relation *n - 1 */
          mpz_swap(R[i], R[(*n)-1]);
          mpz_swap(A[i], A[(*n)-1]);
          tmp = b_status_r[i]; b_status_r[i] = b_status_r[(*n)-1] ; b_status_r[(*n)-1] = tmp;
          tmp = b_status_a[i]; b_status_a[i] = b_status_a[(*n)-1] ; b_status_a[(*n)-1] = tmp;
          atmp = a[i]; a[i] = a[(*n)-1]; a[(*n)-1] = atmp;
          btmp = b[i]; b[i] = b[(*n)-1]; b[(*n)-1] = btmp;
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
          /* relation i is useless, swap it with relation *n - 1 */
          mpz_swap(R[i], R[(*n)-1]);
          mpz_swap(A[i], A[(*n)-1]);
          tmp = b_status_r[i]; b_status_r[i] = b_status_r[(*n)-1] ; b_status_r[(*n)-1] = tmp;
          tmp = b_status_a[i]; b_status_a[i] = b_status_a[(*n)-1] ; b_status_a[(*n)-1] = tmp;
          atmp = a[i]; a[i] = a[(*n)-1]; a[(*n)-1] = atmp;
          btmp = b[i]; b[i] = b[(*n)-1]; b[(*n)-1] = btmp;
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
            /* relation i is smooth, swap it with relation *nb_rel_smooth-1 */
            mpz_swap(R[i], R[*nb_rel_smooth]);
            mpz_swap(A[i], A[*nb_rel_smooth]);
            tmp = b_status_r[i]; b_status_r[i] = b_status_r[*nb_rel_smooth] ; b_status_r[*nb_rel_smooth] = tmp;
            tmp = b_status_a[i]; b_status_a[i] = b_status_a[*nb_rel_smooth] ; b_status_a[*nb_rel_smooth] = tmp;
            atmp = a[i]; a[i] = a[*nb_rel_smooth]; a[*nb_rel_smooth] = atmp;
            btmp = b[i]; b[i] = b[*nb_rel_smooth]; b[*nb_rel_smooth] = btmp;
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
          /* relation i is useless, swap it with relation *n - 1 */
          tmp = b_status_r[i]; b_status_r[i] = b_status_r[(*n)-1] ; b_status_r[(*n)-1] = tmp;
          tmp = b_status_a[i]; b_status_a[i] = b_status_a[(*n)-1] ; b_status_a[(*n)-1] = tmp;
          atmp = a[i]; a[i] = a[(*n)-1]; a[(*n)-1] = atmp;
          btmp = b[i]; b[i] = b[(*n)-1]; b[(*n)-1] = btmp;
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
            /* relation i is smooth, swap it with relation *nb_rel_smooth-1 */
            tmp = b_status_r[i]; b_status_r[i] = b_status_r[*nb_rel_smooth] ; b_status_r[*nb_rel_smooth] = tmp;
            tmp = b_status_a[i]; b_status_a[i] = b_status_a[*nb_rel_smooth] ; b_status_a[*nb_rel_smooth] = tmp;
            atmp = a[i]; a[i] = a[*nb_rel_smooth]; a[*nb_rel_smooth] = atmp;
            btmp = b[i]; b[i] = b[*nb_rel_smooth]; b[*nb_rel_smooth] = btmp;
            (*nb_rel_smooth)++;
          }
          nb_type[4]++;
        }
      }
    }

  mpz_clear(z_B3);
  mpz_clear(z_L2);
}

void
cofac_list_init (cofac_list l)
{
  l->a = NULL;
  l->b = NULL;
  l->R = NULL;
  l->A = NULL;
  l->alloc = 0;
  l->size = 0;
}

void
cofac_list_realloc (cofac_list l, size_t newsize)
{
  unsigned long i;

  /* if we shrink the list, clear the mpz_t's */
  for (i = newsize; i < l->size; i++)
    {
      mpz_clear (l->R[i]);
      mpz_clear (l->A[i]);
    }
  l->a = realloc (l->a, newsize * sizeof (int64_t));
  l->b = realloc (l->b, newsize * sizeof (uint64_t));
  l->R = realloc (l->R, newsize * sizeof (mpz_t));
  l->A = realloc (l->A, newsize * sizeof (mpz_t));
  l->alloc = newsize;
  if (newsize < l->size)
    l->size = newsize;
}

void
cofac_list_add (cofac_list l, long a, unsigned long b, mpz_t R, mpz_t A)
{
  if (l->size == l->alloc)
    cofac_list_realloc (l, 2 * l->alloc + 1);
  l->a[l->size] = a;
  l->b[l->size] = b;
  mpz_init (l->R[l->size]);
  mpz_init (l->A[l->size]);
  mpz_swap (l->R[l->size], R);
  mpz_swap (l->A[l->size], A);
  (l->size)++;
}

void
cofac_list_clear (cofac_list l)
{
  unsigned long i;
  for (i = 0; i < l->size; i++)
    {
      mpz_clear (l->R[i]);
      mpz_clear (l->A[i]);
    }
  free (l->a);
  free (l->b);
  free (l->R);
  free (l->A);
}

/* estimate the bit-size of l[0]*...*l[n-1] */
unsigned long
estimate_bit_size (mpz_t *l, unsigned long n)
{
  double x = 1.0;
  unsigned long i, e = 0;

  for (i = 0; i < n; i++)
    {
      x *= mpz_get_d (l[i]);
      if (x > 1.34e154)
        {
          e += 512;
          x = ldexp (x, -512);
        }
    }
  return e + (unsigned long) ilogb (x) + 1;
}

/* return the number n of smooth relations in l,
   which should be at the end in locations 0, 1, ..., n-1 */
void
find_smooth (cofac_list l, int lpba, int lpbr,
             unsigned long alim, unsigned long rlim)
{
  unsigned long nb_rel;
  unsigned long nb_rel_new;
  unsigned long nb_rel_read = l->size;
  unsigned long nb_rel_step = 10000000;
  unsigned long nb_rel_unknown;
  unsigned long i;
  mpz_t P;
  double s;
  double t_smooth = 0;
  double t_update = 0;
  double start;
  unsigned int n0_pass;
  unsigned long rlim_step;
  unsigned long alim_step;
  unsigned long rlim_new;
  unsigned long alim_new;
  unsigned char *b_status_r;
  unsigned char *b_status_a;
  unsigned long nb_rel_smooth;
  unsigned long nb_smooth_r;
  unsigned long nb_smooth_a;
  unsigned long nb_useless;
  prime_info pi;
  unsigned long prime;
  unsigned int nb_type[5];

  start = seconds ();

  mpz_init (P);

  b_status_r = (unsigned char *) malloc(nb_rel_read * sizeof(unsigned char));
  b_status_a = (unsigned char *) malloc(nb_rel_read * sizeof(unsigned char));
  for (i = 0; i < nb_rel_read; i++)
  {
    b_status_r[i] = STATUS_UNKNOWN;
    b_status_a[i] = STATUS_UNKNOWN;
  }

  nb_rel_smooth = 0;
  nb_rel_unknown = nb_rel_read;

  nb_smooth_r = 0;
  nb_smooth_a = 0;
  nb_useless = 0;

  prime_info_init (pi);

  /* Initial one-side pass (to make rlim and alim be equal) */

  unsigned long rsize, asize;
  /* the product of primes in [rlim, rlim + step] has a number of bits about:
     sum(log(p)/log(2)/log(p), p = rlim..rlim + step) ~ step/log(2)
     thus we need step ~ rsize*log(2) */
  rsize = estimate_bit_size (l->R, l->size);
  rlim_step = (unsigned long) ((double) rsize * log (2.0));
  asize = estimate_bit_size (l->A, l->size);
  alim_step = (unsigned long) ((double) asize * log (2.0));
  fprintf (stderr, "rlim_step = %lu alim_step = %lu\n", rlim_step, alim_step);

  /* the code below assumes max(rlim,alim) <= min(2^lpbr,2^lpba) */
  ASSERT_ALWAYS(rlim <= (1UL << lpbr));
  ASSERT_ALWAYS(rlim <= (1UL << lpba));
  ASSERT_ALWAYS(alim <= (1UL << lpbr));
  ASSERT_ALWAYS(alim <= (1UL << lpba));

  if (rlim < alim)
  {
    for (prime = 2; prime <= rlim; prime = getprime_mt (pi));

    while (rlim < alim)
    {
      rlim_new = (rlim + rlim_step < alim ? rlim + rlim_step : alim);

      fprintf (stderr, "\nInitial one-side pass: %0.fs\n\n", seconds() - start);

      fprintf (stderr, "rlim: %lu:%lu\n\n", rlim, rlim_new);

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
        smoothness_test (&(l->R[nb_rel]), nb_rel_new - nb_rel, P);
        t_smooth += seconds();
        fprintf (stderr, "smoothness_test (%lu cofactors) took %.0f seconds"
                 " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
        nb_rel = nb_rel_new;
      }

      rlim = rlim_new;
      t_update -= seconds();
      update_status (l->R, l->A, b_status_r, b_status_a, &nb_rel_unknown,
                     &nb_rel_smooth, rlim, lpbr, &nb_smooth_r, &nb_smooth_a,
                     &nb_useless, nb_type, l->a, l->b);
      t_update += seconds();
      fprintf (stderr, "nb_smooth_r = %lu ; nb_useless = %lu ; nb_unknown = %lu ; nb_rel_smooth = %lu ;"
               " %u: %u : %u : %u : %u\n",
               nb_smooth_r, nb_useless, nb_rel_read - nb_smooth_r - nb_useless, nb_rel_smooth,
             nb_type[0], nb_type[1], nb_type[2], nb_type[3], nb_type[4]);
      fprintf (stderr, "t_update: %.0f seconds\n", t_update);
    }
  }
  else if (alim < rlim)
  {
    for (prime = 2; prime <= alim; prime = getprime_mt (pi));

    while (alim < rlim)
    {
      alim_new = (alim + alim_step < rlim ? alim + alim_step : rlim);

      fprintf (stderr, "\nInitial one-side pass: %0.fs\n\n", seconds() - start);

      fprintf (stderr, "\nalim: %lu:%lu\n", alim, alim_new);

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
        smoothness_test (&(l->A[nb_rel]), nb_rel_new - nb_rel, P);
        t_smooth += seconds();
        fprintf (stderr, "smoothness_test (%lu cofactors) took %.0f seconds"
                 " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
        nb_rel = nb_rel_new;
      }

      alim = alim_new;
      t_update -= seconds();
      update_status (l->A, l->R, b_status_a, b_status_r, &nb_rel_unknown,
                     &nb_rel_smooth, alim, lpba, &nb_smooth_a, &nb_smooth_r,
                     &nb_useless, nb_type, l->a, l->b);
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
    for (prime = 2; prime <= rlim; prime = getprime_mt (pi));
  }

  ASSERT_ALWAYS (rlim == alim);

  /* the loop below assumes rlim_step = alim_step */
  rlim_step = alim_step = (rlim_step + alim_step) / 2;

  /* Loop */

  n0_pass = 0;
  while ( (rlim < (1UL << lpbr)) || (alim < (1UL << lpba)) )
  {
    n0_pass++;

    fprintf (stderr, "\nPass %u: %.0f s\n", n0_pass, seconds() - start);

    /* rational side */

    rlim_new = (rlim + rlim_step < (1UL << lpbr) ? rlim + rlim_step : (1UL << lpbr));

    fprintf (stderr, "\nrlim: %lu:%lu\n", rlim, rlim_new);

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
      smoothness_test (&(l->R[nb_rel]), nb_rel_new - nb_rel, P);
      t_smooth += seconds();
      fprintf (stderr, "smoothness_test (%lu cofactors) took %.0f seconds"
               " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
      nb_rel = nb_rel_new;
    }

    rlim = rlim_new;
    t_update -= seconds();
    update_status (l->R, l->A, b_status_r, b_status_a, &nb_rel_unknown,
                   &nb_rel_smooth, rlim, lpbr, &nb_smooth_r, &nb_smooth_a,
                   &nb_useless, nb_type, l->a, l->b);
    t_update += seconds();
    fprintf (stderr, "nb_smooth_r = %lu ; nb_useless = %lu ; nb_unknown = %lu ; nb_rel_smooth = %lu ;"
             " %u : %u : %u : %u : %u\n",
             nb_smooth_r, nb_useless, nb_rel_read - nb_smooth_r - nb_useless, nb_rel_smooth,
           nb_type[0], nb_type[1], nb_type[2], nb_type[3], nb_type[4]);
    fprintf (stderr, "t_update: %.0f seconds\n", t_update);

    /* algebraic side */

    alim_new = (alim + alim_step < (1UL << lpba) ? alim + alim_step : (1UL << lpba));

    fprintf (stderr, "\nalim: %lu:%lu\n", alim, alim_new);

    nb_rel = nb_rel_smooth;
    while (nb_rel < nb_rel_unknown)
    {
      nb_rel_new = (nb_rel + nb_rel_step < nb_rel_unknown ? nb_rel + nb_rel_step : nb_rel_unknown);
      s = seconds ();
      t_smooth -= seconds();
      smoothness_test (&(l->A[nb_rel]), nb_rel_new - nb_rel, P);
      t_smooth += seconds();
      fprintf (stderr, "smoothness_test (%lu cofactors) took %.0f seconds"
               " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
      nb_rel = nb_rel_new;
    }

    alim = alim_new;
    t_update -= seconds();
    update_status (l->A, l->R, b_status_a, b_status_r, &nb_rel_unknown,
                   &nb_rel_smooth, alim, lpba, &nb_smooth_a, &nb_smooth_r,
                   &nb_useless, nb_type, l->a, l->b);
    t_update += seconds();
    fprintf (stderr, "nb_smooth_a = %lu ; nb_useless = %lu ; nb_unknown = %lu ; nb_rel_smooth = %lu ;"
             " %u : %u : %u : %u : %u\n",
             nb_smooth_a, nb_useless, nb_rel_read - nb_smooth_a - nb_useless, nb_rel_smooth,
             nb_type[0], nb_type[1], nb_type[2], nb_type[3], nb_type[4]);
    fprintf (stderr, "t_update: %.0f seconds\n", t_update);
  }
  prime_info_clear (pi);

  fprintf (stderr, "\nCollect smooth relations: %.0f s\n\n", seconds() - start);

#if 0
  for (i = 0; i < nb_rel_smooth; i++)
  {
    ASSERT_ALWAYS ( (b_status_r[i] == STATUS_SMOOTH) && (b_status_a[i] == STATUS_SMOOTH) );
    printf ("Smooth: a=%" PRId64 " b=%" PRIu64 "\n", l->a[i], l->b[i]);
  }
#endif
  fprintf (stderr, "\nFound %lu smooth relations\n", nb_rel_smooth);

  cofac_list_realloc (l, nb_rel_smooth);

  mpz_clear (P);
  free (b_status_r);
  free (b_status_a);
}

void declare_usage (param_list pl)
{
  param_list_decl_usage (pl, "poly", "Polynomial file");
  param_list_decl_usage (pl, "lpbr", "rational large prime bound is 2^lpbr");
  param_list_decl_usage (pl, "lpba", "algebraic large prime bound is 2^lpba");
  param_list_decl_usage (pl, "rlim", "rational factor base bound");
  param_list_decl_usage (pl, "alim", "algebraic factor base bound");
}

void
compute_biproduct (mpz_t Q, cofac_list L)
{
  unsigned long n = L->size, i, j;
  mpz_t T[MAX_DEPTH];
  unsigned long acc[MAX_DEPTH] = {0,};

  for (i = 0; i < MAX_DEPTH; i++)
    mpz_init_set_ui (T[i], 1);
  acc[0] = 0;

  /* Each T[i] accumulates so far acc[i] pairs of norms,
     with acc[i] <= 2^(i+1).
     When acc[i] reaches 2^(i+1) values, we move them to T[i+1] */
  for (i = 0; i < n; i++)
    {
      mpz_mul (Q, L->R[i], L->A[i]);
      mpz_mul (T[0], T[0], Q);
      acc[0] ++;
      j = 0;
      for (j = 0; acc[j] == (2UL << j); j++)
        {
          ASSERT_ALWAYS (j+1 < MAX_DEPTH);
          mpz_mul (T[j+1], T[j+1], T[j]);
          acc[j+1] += acc[j];
          mpz_set_ui (T[j], 1);
          acc[j] = 0;
        }
    }
  /* final accumulation */
  mpz_set_ui (Q, 1);
  for (j = 0; j < MAX_DEPTH; j++)
    mpz_mul (Q, Q, T[j]);

  for (i = 0; i < MAX_DEPTH; i++)
    mpz_clear (T[i]);
}

/* return the number of characters read, -1 if error */
int
print_norm (char *s, mpz_t N, mpz_t *P, unsigned long n, unsigned int lpb)
{
  unsigned long j;
  int len = 0; /* number of characters printed */

  mpz_abs (N, N);
  if (mpz_cmp_ui (N, 1) == 0)
    return 0;

  if (mpz_sizeinbase (N, 2) <= lpb)
    return gmp_sprintf (s, "%Zx", N);

  for (j = 0; j < n; j++)
    if (mpz_divisible_p (N, P[j]))
      {
        do {
          if (len > 0)
            len += sprintf (s + len, ",");
          len += gmp_sprintf (s + len, "%Zx", P[j]);
          mpz_divexact (N, N, P[j]);
        }
        while (mpz_divisible_p (N, P[j]));
        /* if N < 2^lpb, it is necessarily prime */
        if (mpz_sizeinbase (N, 2) <= lpb)
          {
            len += gmp_sprintf (s + len, ",%Zx", N);
            return len;
          }
      }
  if (mpz_cmp_ui (N, 1) != 0)
    {
      gmp_fprintf (stderr, "Error, remaining unfactored part %Zd\n", N);
      return -1;
    }
  return len;
}

void
print_relations (cofac_list L, mpz_t *P, unsigned long n, int lpbr, int lpba)
{
  unsigned long i;
  int ret, len;
  double s = seconds ();
  char line[1024];
  unsigned long printed = 0;

  for (i = 0; i < L->size; i++)
    {
      len = sprintf (line, "%ld:%lu:", L->a[i], L->b[i]);
      ret = print_norm (line + len, L->R[i], P, n, lpbr);
      if (ret >= 0)
        {
          len += ret;
          len += sprintf (line + len, ":");
          ret = print_norm (line + len, L->A[i], P, n, lpba);
        }
      if (ret >= 0)
        {
          printf ("%s\n", line);
          printed ++;
        }
      else
        fprintf (stderr, "Error for %ld:%lu\n", L->a[i], L->b[i]);
    }
  fprintf (stderr, "Printed %lu relations in %.0f s\n",
           printed, seconds () - s);
}

/* Given a list L of bi-smooth cofactors, print the corresponding relations
   on stdout. The algorithm is as follows:
   (1) compute the product Q of all norms F(a,b) and G(a,b) corresponding
       to pairs (a,b) in L
   (2) construct a product tree of all primes p up to max(2^lpba, 2^lpbr)
   (3) going down that product tree, compute the remainder of Q modulo all
       primes p
   (4) store in a list T the primes p for which Q mod p = 0, these are the only
       primes that can appear in any F(a,b) or G(a,b)
   (5) use trial division (with primes in T) to factor the norms
*/
void
factor (cofac_list L, const char *poly_file, int lpba, int lpbr)
{
  cado_poly pol;
  unsigned long n = L->size, i, j, nprimes, w[MAX_DEPTH];
  mpz_t Q, *LP, **T;
  double s;

  cado_poly_init (pol);
  if (cado_poly_read (pol, poly_file) == 0)
    {
      fprintf (stderr, "Could not read polynomial file\n");
      exit (1);
    }

  fprintf (stderr, "factor: computing full norms...");
  fflush (stderr);
  s = seconds ();
  /* compute all norms F(a,b) and G(a,b) */
  for (i = 0; i < n; i++)
    {
      mpz_poly_homogeneous_eval_siui (L->A[i], pol->alg, L->a[i], L->b[i]);
      mpz_poly_homogeneous_eval_siui (L->R[i], pol->rat, L->a[i], L->b[i]);
    }
  fprintf (stderr, "done in %.0f s\n", seconds () - s);
  fflush (stderr);

  /* compute the product of all norms */
  mpz_init (Q);
  fprintf (stderr, "factor: computing product of all norms...");
  fflush (stderr);
  s = seconds ();
  compute_biproduct (Q, L);
  fprintf (stderr, "done in %.0f s\n", seconds () - s);
  fflush (stderr);

  /* compute all primes up to max(2^lpba, 2^lpbr) */
  lpba = (lpba > lpbr) ? lpba : lpbr;
  fprintf (stderr, "factor: computing primes up to %lu...", 1UL << lpba);
  fflush (stderr);
  s = seconds ();
  LP = prime_list (1UL << lpba, &nprimes);
  fprintf (stderr, "done in %.0f s\n", seconds () - s);
  fflush (stderr);

  /* form a product tree from LP */
  fprintf (stderr, "factor: computing product tree of primes...");
  fflush (stderr);
  s = seconds ();
  T = product_tree (LP, nprimes, w);
  fprintf (stderr, "done in %.0f s\n", seconds () - s);
  fflush (stderr);

  /* compute the remainder tree Q mod T */
  fprintf (stderr, "factor: computing remainder tree...");
  fflush (stderr);
  s = seconds ();
  remainder_tree (T, nprimes, w, Q);
  fprintf (stderr, "done in %.0f s\n", seconds () - s);
  fflush (stderr);

  /* now scan all primes appearing in norms */
  fprintf (stderr, "factor: scan smooth relations...");
  fflush (stderr);
  s = seconds ();
  for (i = j = 0; j < nprimes / 2; j++)
    {
      mpz_mod (Q, T[1][j], LP[2*j]);
      if (mpz_cmp_ui (Q, 0) == 0)
        {
          mpz_swap (LP[i], LP[2*j]);
          i++;
        }
      mpz_mod (Q, T[1][j], LP[2*j+1]);
      if (mpz_cmp_ui (Q, 0) == 0)
        {
          mpz_swap (LP[i], LP[2*j+1]);
          i++;
        }
    }
  if (nprimes & 1)
    {
      mpz_mod (Q, T[1][j], LP[2*j]);
      if (mpz_cmp_ui (Q, 0) == 0)
        {
          mpz_swap (LP[i], LP[2*j]);
          i++;
        }
    }
  fprintf (stderr, "done: found %lu primes in %.0f s\n", i, seconds () - s);
  fflush (stderr);
  /* the primes appearing in relations are LP[0..i-1] */

  print_relations (L, LP, i, lpbr, lpba);

  clear_product_tree (T, nprimes, w);
  for (i = 0; i < nprimes; i++)
    mpz_clear (LP[i]);
  free (LP);
  mpz_clear (Q);
  cado_poly_clear (pol);
}

int
main (int argc, char* argv[])
{
  cofac_list L;
  double start;
  FILE *cofac;
  int64_t a;
  uint64_t b;
  mpz_t R, A;
  param_list pl;
  const char *poly_file;
  int lpba, lpbr;
  unsigned long alim, rlim;

  start = seconds ();

  param_list_init (pl);
  declare_usage (pl);

  argc--, argv++;
  for( ; argc ; )
    {
      if (param_list_update_cmdline (pl, &argc, &argv))
        continue;
      else
        break;
    }

  poly_file = param_list_lookup_string (pl, "poly");
  ASSERT_ALWAYS(param_list_parse_int (pl, "lpba", &lpba));
  ASSERT_ALWAYS(param_list_parse_int (pl, "lpbr", &lpbr));
  ASSERT_ALWAYS(param_list_parse_ulong (pl, "alim", &alim));
  ASSERT_ALWAYS(param_list_parse_ulong (pl, "rlim", &rlim));

  if (poly_file == NULL)
    {
      fprintf (stderr, "Error, missing -poly <file>\n");
      exit (1);
    }

  /* Initialization */
  cofac_list_init (L);

  ASSERT_ALWAYS (argc == 1);
  cofac = fopen (argv[0], "r");

  mpz_init (R);
  mpz_init (A);
  while (1)
  {
    int ret;
    ret = gmp_fscanf (cofac, "%ld %lu %Zd %Zd\n", &a, &b, R, A);
    if (ret != 4)
      break;
    cofac_list_add (L, a, b, R, A);
  }
  fclose (cofac);
  cofac_list_realloc (L, L->size);
  fprintf (stderr, "Read %lu cofactors in %.0fs\n", L->size,
           seconds () - start);
  fflush (stderr);

  start = seconds ();
  find_smooth (L, lpba, lpbr, alim, rlim);
  fprintf (stderr, "Detecting smooth cofactors took %.0f s\n",
           seconds() - start);

  start = seconds ();
  factor (L, poly_file, lpba, lpbr);
  fprintf (stderr, "Factoring smooth cofactors took %.0f s\n",
           seconds() - start);

  mpz_clear (R);
  mpz_clear (A);
  cofac_list_clear (L);

  param_list_clear (pl);

  return 0;
}
