/* Bernstein's smoothness test */

/* This is a work-in-progress implementation. To use it:
   1) create a file (say cofac) containing lines of the following form:
      a b cofac_0 cofac_1
      where cofac_0 and cofac_1 are cofactors on sides 0 and 1 respectively
   2) run ./smooth -poly ... -cofac cofac -lpb0 ... -lpb1 ... -lim0 ...
      -lim1 ... [-v] [-split ...]

   The algorithm is the following:
   (a) compute the product P of all primes in [B, L]
   (b) compute the product tree T of all cofactors
   (c) compute the remainder tree of P mod T: at the leaves we have 0 for
       smooth cofactors (assuming all primes in [B,L] have multiplicity 1)
   In practice we split P into smaller chunks P_j of bit-size near that of all
   cofactors, and compute a gcd at leaves to reveal primes in P_j.
   If the number of cofactors doubles, the cost goes from M(n)*log(n) to
   M(2n)*log(2n), thus essentially doubles.
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <gmp.h>
#include "facul.h"
#include "facul_doit.h"
#include "smooth.h"
#include "utils.h"
#include "portability.h"

#define NB_MILLER_RABIN 1

#define STATUS_SMOOTH  0
#define STATUS_UNKNOWN 1
#define STATUS_USELESS 2

static void
mpz_list_init (mpz_list L)
{
  L->l = NULL;
  L->alloc = L->size = 0;
}

/* at any point, all values up to L->alloc should be mpz_init()'ed,
   even beyond L->size */
static void
mpz_list_add (mpz_list L, unsigned long n)
{
  if (L->size == L->alloc)
    {
      unsigned long old = L->alloc, i;
      L->alloc = 3 * (L->alloc / 2) + 2;
      L->l = realloc (L->l, L->alloc * sizeof (mpz_t));
      for (i = old; i < L->alloc; i++)
        mpz_init (L->l[i]);
    }
  ASSERT_ALWAYS(L->size < L->alloc);
  mpz_set_ui (L->l[L->size], n);
  L->size ++;
}

static void
mpz_list_clear (mpz_list L)
{
  size_t i;

  for (i = 0; i < L->alloc; i++)
    mpz_clear (L->l[i]);
  free (L->l);
  L->alloc = L->size = 0;
}

static unsigned long
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
static mpz_t**
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
static void
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
static void
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
   Each R[j] has been divided by its P-smooth part.
*/
static void
smoothness_test (mpz_t *R, unsigned long n, mpz_t P, int verbose)
{
  unsigned long h = tree_height (n), j, w[MAX_DEPTH];
  mpz_t **T;
  mpz_t w1, w2;

  mpz_init (w1);
  mpz_init (w2);

  T = product_tree (R, n, w);

  if (verbose)
    fprintf (stderr, "cofactor product has %lu bits\n",
             mpz_sizeinbase (T[h][0], 2));

  /* compute remainder tree */
  remainder_tree (T, n, w, P);

  /* special last loop for i=1, with T[0] = R */
  for (j = 0; j < n / 2; j++)
    {
      mpz_mod (w1, T[1][j], R[2*j]);
      mpz_gcd (w2, w1, R[2*j]);
      mpz_divexact (R[2*j], R[2*j], w2);
      mpz_mod (w1, T[1][j], R[2*j+1]);
      mpz_gcd (w2, w1, R[2*j+1]);
      mpz_divexact (R[2*j+1], R[2*j+1], w2);
    }
  /* special case if n is odd */
  if (n & 1)
    {
      mpz_gcd (w2, T[1][w[1]-1], R[n-1]);
      mpz_divexact (R[n-1], R[n-1], w2);
    }

  mpz_clear (w1);
  mpz_clear (w2);
  clear_product_tree (T, n, w);
}

/* add in the list L all primes pmin <= p < pmax.
   Assume pmin is the current prime in 'pi'
   (pmin=2 when 'pi' was just initialized).
   Return the current prime in 'pi' at the end, i.e.,
   the smallest prime >= pmax. */
static unsigned long
prime_list (mpz_list L, prime_info pi, unsigned long pmin,
            unsigned long pmax)
{
  unsigned long p;

  for (p = pmin; p < pmax; p = getprime_mt (pi))
    mpz_list_add (L, p);
  return p;
}

/* FIXME: since we don't need to keep the indidivual primes here,
   instead of allocating spaces for n mpz_t data structures, we need
   only to allocate O(log n) [see function compute_biproduct below] */
static unsigned long
prime_product (mpz_t P, prime_info pi, unsigned long p_max,
               unsigned long p_last)
{
  unsigned long i;
  mpz_list L;

  mpz_list_init (L);
  p_last = prime_list (L, pi, p_last, p_max);

  /* FIXME: equilibrate the product */
  mpz_t *l = L->l;
  unsigned long n = L->size;
  while (n > 1)
  {
    for (i = 0; i+1 < n; i+=2)
      mpz_mul (l[i/2], l[i], l[i+1]);
    if (n & 1)
      mpz_swap (l[n/2], l[n-1]);
    n = (n + 1) / 2;
  }
  mpz_set (P, l[0]);

  mpz_list_clear (L);
  return p_last;
}

/* invariant:
   relations 0 to *nb_rel_smooth-1 are smooth
   relations *nb_rel_smooth to *n-1 are unknown
   relations >= *n are non-smooth (useless)
*/
static void
update_status (mpz_t *R, mpz_t *A,
               unsigned char *b_status_r, unsigned char *b_status_a,
               unsigned long *n, unsigned long *nb_rel_smooth,
               unsigned long lim0, unsigned long lpb0,
               unsigned long *nb_smooth_r, unsigned long *nb_smooth_a,
               unsigned long *nb_useless, int64_t *a, uint64_t *b)
{
  mpz_t z_B2, z_BL, z_B3, z_L2;

  unsigned long tmp;
  unsigned long i;
  unsigned long B;
  unsigned long L;

  int64_t atmp;
  uint64_t btmp;

  mpz_init(z_B2);
  mpz_init(z_BL);
  mpz_init(z_B3);
  mpz_init(z_L2); /* set to 0 */

  mpz_set_ui (z_B2, lim0);
  mpz_mul_ui (z_B2, z_B2, lim0);
  mpz_mul_ui (z_B3, z_B2, lim0);
  mpz_setbit (z_L2, 2 * lpb0);
  B = lim0;
  L = 1UL << lpb0;
  mpz_set_ui (z_BL, B);
  mpz_mul_ui (z_BL, z_BL, L);
  ASSERT_ALWAYS(mpz_cmp_ui (z_B2, L) >= 0);

  for (i = *nb_rel_smooth; i < *n; i++)
    {
      if (b_status_r[i] == STATUS_UNKNOWN)
      {
        /* if L^2 < R[i] < B^3 or L < R[i] < B^2, then R[i] cannot be smooth */
        if ((mpz_cmp (R[i], z_L2) > 0 && mpz_cmp (R[i], z_B3) < 0) ||
            (mpz_cmp_ui (R[i], L) > 0 && mpz_cmp (R[i], z_B2) < 0))
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
        }
        /* if R[i] < L, then R[i] is smooth (we assume L <= B^2) */
        else if (mpz_cmp_ui (R[i], L) <= 0)
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
        }
        /* now L <= B^2 <= R[i] <= L^2 or B^3 <= R[i] */
        else if (mpz_probab_prime_p (R[i], NB_MILLER_RABIN) != 0)
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
        }
        /* now L <= B^2 <= R[i] <= L^2 or B^3 <= R[i] and R[i] is composite */
        else if (mpz_cmp (R[i], z_BL) <= 0)
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
        }
      }
    }

  mpz_clear(z_B2);
  mpz_clear(z_BL);
  mpz_clear(z_B3);
  mpz_clear(z_L2);
}

static void
cofac_list_init (cofac_list l)
{
  l->a = NULL;
  l->b = NULL;
  l->R = NULL;
  l->A = NULL;
  l->alloc = 0;
  l->size = 0;
}

static void
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

static void
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

static void
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
static unsigned long
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
static void
find_smooth (cofac_list l, int lpb1, int lpb0,
             unsigned long lim1, unsigned long lim0, int split, int verbose)
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
  double t_prime = 0;
  double start;
  unsigned int n0_pass;
  unsigned long lim0_step;
  unsigned long lim1_step;
  unsigned long lim0_new;
  unsigned long lim1_new;
  unsigned char *b_status_r;
  unsigned char *b_status_a;
  unsigned long nb_rel_smooth;
  unsigned long nb_smooth_r;
  unsigned long nb_smooth_a;
  unsigned long nb_useless;
  prime_info pi;
  unsigned long prime;

  start = seconds ();

  mpz_init (P);

  b_status_r = (unsigned char *) malloc(nb_rel_read * sizeof(unsigned char));
  b_status_a = (unsigned char *) malloc(nb_rel_read * sizeof(unsigned char));
  for (i = 0; i < nb_rel_read; i++)
  {
    b_status_r[i] = STATUS_UNKNOWN;
    b_status_a[i] = STATUS_UNKNOWN;
    ASSERT(mpz_cmp_ui (l->R[i], 0) > 0);
    ASSERT(mpz_cmp_ui (l->A[i], 0) > 0);
  }

  nb_rel_smooth = 0;
  nb_rel_unknown = nb_rel_read;

  nb_smooth_r = 0;
  nb_smooth_a = 0;
  nb_useless = 0;

  prime_info_init (pi);

  /* Initial one-side pass (to make lim0 and lim1 be equal) */

  unsigned long rsize, asize;
  /* the product of primes in [lim0, lim0 + step] has a number of bits about:
     sum(log(p)/log(2)/log(p), p = lim0..lim0 + step) ~ step/log(2)
     thus we need step ~ rsize*log(2) */
  if (split == 0)
    {
      rsize = estimate_bit_size (l->R, l->size);
      lim0_step = (unsigned long) ((double) rsize * log (2.0));
      asize = estimate_bit_size (l->A, l->size);
      lim1_step = (unsigned long) ((double) asize * log (2.0));
    }
  else
    {
       lim0_step = 1 + ((uint64_t) 1UL << lpb0) / split;
       lim1_step = 1 + ((uint64_t) 1UL << lpb1) / split;
    }

  /* the code below assumes max(lim0,lim1) <= min(2^lpb0,2^lpb1) */
  ASSERT_ALWAYS(lim0 <= (1UL << lpb0));
  ASSERT_ALWAYS(lim0 <= (1UL << lpb1));
  ASSERT_ALWAYS(lim1 <= (1UL << lpb0));
  ASSERT_ALWAYS(lim1 <= (1UL << lpb1));

  if (lim0 < lim1)
  {
    for (prime = 2; prime <= lim0; prime = getprime_mt (pi));

    while (lim0 < lim1)
    {
      lim0_new = (lim0 + lim0_step < lim1 ? lim0 + lim0_step : lim1);

      if (verbose)
        {
          fprintf (stderr, "\nInitial one-side pass: %0.fs\n\n", seconds() - start);
          fprintf (stderr, "lim0: %lu:%lu\n\n", lim0, lim0_new);
        }

      s = seconds ();
      prime = prime_product (P, pi, lim0_new, prime);
      s = seconds () - s;
      t_prime += s;
      if (verbose)
        fprintf (stderr, "Computing prime product of %lu bits took %.0fs (total %.0fs so far)\n",
                 mpz_sizeinbase (P, 2), s, t_prime);

      nb_rel = nb_rel_smooth;
      while (nb_rel < nb_rel_unknown)
      {
        nb_rel_new = (nb_rel + nb_rel_step < nb_rel_unknown ? nb_rel + nb_rel_step : nb_rel_unknown);
        s = seconds ();
        t_smooth -= seconds();
        smoothness_test (&(l->R[nb_rel]), nb_rel_new - nb_rel, P, verbose);
        t_smooth += seconds();
        if (verbose)
          fprintf (stderr, "smoothness test (%lu cofactors) took %.0f seconds"
                   " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
        nb_rel = nb_rel_new;
      }

      lim0 = lim0_new;
      t_update -= seconds();
      update_status (l->R, l->A, b_status_r, b_status_a, &nb_rel_unknown,
                     &nb_rel_smooth, lim0, lpb0, &nb_smooth_r, &nb_smooth_a,
                     &nb_useless, l->a, l->b);
      t_update += seconds();
      if (verbose)
        {
          fprintf (stderr, "smooth_r:%lu useless:%lu unknown:%lu rel_smooth:%lu\n",
                   nb_smooth_r, nb_useless, nb_rel_read - nb_smooth_r - nb_useless, nb_rel_smooth);
          fprintf (stderr, "t_update: %.0f seconds\n", t_update);
        }
    }
  }
  else if (lim1 < lim0)
  {
    for (prime = 2; prime <= lim1; prime = getprime_mt (pi));

    while (lim1 < lim0)
    {
      lim1_new = (lim1 + lim1_step < lim0 ? lim1 + lim1_step : lim0);

      if (verbose)
        {
          fprintf (stderr, "\nInitial one-side pass: %0.fs\n\n", seconds() - start);
          fprintf (stderr, "\nlim1: %lu:%lu\n", lim1, lim1_new);
        }

      s = seconds ();
      prime = prime_product (P, pi, lim1_new, prime);
      s = seconds () - s;
      t_prime += s;
      if (verbose)
        fprintf (stderr, "Computing prime product of %lu bits took %.0fs (total %.0fs so far)\n",
                 mpz_sizeinbase (P, 2), s, t_prime);

      nb_rel = nb_rel_smooth;
      while (nb_rel < nb_rel_unknown)
      {
        nb_rel_new = (nb_rel + nb_rel_step < nb_rel_unknown ? nb_rel + nb_rel_step : nb_rel_unknown);
        s = seconds ();
        t_smooth -= seconds();
        smoothness_test (&(l->A[nb_rel]), nb_rel_new - nb_rel, P, verbose);
        t_smooth += seconds();
        if (verbose)
          fprintf (stderr, "smoothness test (%lu cofactors) took %.0f seconds"
                   " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
        nb_rel = nb_rel_new;
      }

      lim1 = lim1_new;
      t_update -= seconds();
      update_status (l->A, l->R, b_status_a, b_status_r, &nb_rel_unknown,
                     &nb_rel_smooth, lim1, lpb1, &nb_smooth_a, &nb_smooth_r,
                     &nb_useless, l->a, l->b);
      t_update += seconds();
      if (verbose)
        {
          fprintf (stderr, "smooth_a:%lu useless:%lu unknown:%lu rel_smooth:%lu\n",
                   nb_smooth_a, nb_useless, nb_rel_read - nb_smooth_a - nb_useless, nb_rel_smooth);
          fprintf (stderr, "t_update: %.0f seconds\n", t_update);
        }
    }
  }
  else  // lim0 = lim1
  {
    for (prime = 2; prime <= lim0; prime = getprime_mt (pi));
  }

  ASSERT_ALWAYS (lim0 == lim1);

  /* the loop below assumes lim0_step = lim1_step */
  if (lim0_step < lim1_step)
    lim0_step = lim1_step;
  else
    lim1_step = lim0_step;

  /* Loop */

  n0_pass = 0;
  while ( (lim0 < (1UL << lpb0)) || (lim1 < (1UL << lpb1)) )
  {
    n0_pass++;

    if (verbose)
      fprintf (stderr, "\nPass %u: %.0f s\n", n0_pass, seconds() - start);

    /* side 0 */

    lim0_new = (lim0 + lim0_step < (1UL << lpb0) ? lim0 + lim0_step : (1UL << lpb0));

    if (verbose)
      fprintf (stderr, "\nlim0: %lu:%lu\n", lim0, lim0_new);

    s = seconds ();
    prime = prime_product (P, pi, lim0_new, prime);
    s = seconds () - s;
    t_prime += s;
    if (verbose)
      fprintf (stderr, "Computing prime product of %lu bits took %.0fs (total %.0fs so far)\n",
               mpz_sizeinbase (P, 2), s, t_prime);

    nb_rel = nb_rel_smooth;
    while (nb_rel < nb_rel_unknown)
    {
      nb_rel_new = (nb_rel + nb_rel_step < nb_rel_unknown ? nb_rel + nb_rel_step : nb_rel_unknown);
      s = seconds ();
      t_smooth -= seconds();
      smoothness_test (&(l->R[nb_rel]), nb_rel_new - nb_rel, P, verbose);
      t_smooth += seconds();
      if (verbose)
        fprintf (stderr, "smoothness_test (%lu cofactors) took %.0f seconds"
                 " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
      nb_rel = nb_rel_new;
    }

    lim0 = lim0_new;
    t_update -= seconds();
    update_status (l->R, l->A, b_status_r, b_status_a, &nb_rel_unknown,
                   &nb_rel_smooth, lim0, lpb0, &nb_smooth_r, &nb_smooth_a,
                   &nb_useless, l->a, l->b);
    t_update += seconds();
    if (verbose)
      {
        fprintf (stderr, "smooth_r:%lu useless:%lu unknown:%lu rel_smooth:%lu\n",
                 nb_smooth_r, nb_useless, nb_rel_read - nb_smooth_r - nb_useless, nb_rel_smooth);
        fprintf (stderr, "t_update: %.0f seconds\n", t_update);
      }

    /* side 1 */

    lim1_new = (lim1 + lim1_step < (1UL << lpb1) ? lim1 + lim1_step : (1UL << lpb1));

    if (verbose)
      fprintf (stderr, "\nlim1: %lu:%lu\n", lim1, lim1_new);

    nb_rel = nb_rel_smooth;
    while (nb_rel < nb_rel_unknown)
    {
      nb_rel_new = (nb_rel + nb_rel_step < nb_rel_unknown ? nb_rel + nb_rel_step : nb_rel_unknown);
      s = seconds ();
      t_smooth -= seconds();
      smoothness_test (&(l->A[nb_rel]), nb_rel_new - nb_rel, P, verbose);
      t_smooth += seconds();
      if (verbose)
        fprintf (stderr, "smoothness_test (%lu cofactors) took %.0f seconds"
                 " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
      nb_rel = nb_rel_new;
    }

    lim1 = lim1_new;
    t_update -= seconds();
    update_status (l->A, l->R, b_status_a, b_status_r, &nb_rel_unknown,
                   &nb_rel_smooth, lim1, lpb1, &nb_smooth_a, &nb_smooth_r,
                   &nb_useless, l->a, l->b);
    t_update += seconds();
    if (verbose)
      {
        fprintf (stderr, "smooth_a:%lu useless:%lu unknown:%lu rel_smooth:%lu\n",
                 nb_smooth_a, nb_useless, nb_rel_read - nb_smooth_a - nb_useless, nb_rel_smooth);
        fprintf (stderr, "t_update: %.0f seconds\n", t_update);
      }
  }
  prime_info_clear (pi);

  cofac_list_realloc (l, nb_rel_smooth);

  mpz_clear (P);
  free (b_status_r);
  free (b_status_a);

  fprintf (stderr, "t_prime %.0fs, t_smooth %.0fs, t_update %.0fs\n",
           t_prime, t_smooth, t_update);
}

static void
declare_usage (param_list pl)
{
  param_list_decl_usage (pl, "poly", "Polynomial file");
  param_list_decl_usage (pl, "cofac", "Cofactor file");
  param_list_decl_usage (pl, "lpb0", "side-0 large prime bound is 2^lpb0");
  param_list_decl_usage (pl, "lpb1", "side-1 large prime bound is 2^lpb1");
  param_list_decl_usage (pl, "lim0", "side-0 factor base bound");
  param_list_decl_usage (pl, "lim1", "side-1 factor base bound");
  param_list_decl_usage (pl, "v",    "(switch) verbose mode");
  param_list_decl_usage (pl, "split", "number of splits");
}

/* the list SP (small primes) contains all primes < B */
static void
print_smooth (mpz_t *factors, mpz_t n, facul_method_t *methods,
              unsigned long nb_methods,
              struct modset_t *fm, struct modset_t *cfm,
              int lpb, double BB, double BBB, mpz_list SP)
{
  int res_fac, j;
  unsigned long i;
  int m = 0; /* number of printed factors */

  /* remove small primes < B */
  for (i = 0; i < SP->size; i++)
    {
      while (mpz_divisible_p (n, SP->l[i]))
        {
          mpz_divexact (n, n, SP->l[i]);
          if (m++ > 0)
            printf (",");
          gmp_printf ("%Zx", SP->l[i]);
        }
    }

  /* any factor < B^2 is necessarily prime */

  for (i = 0; i < nb_methods && mpz_cmp_d (n, BB) >= 0; i++)
    {
      res_fac = facul_doit_onefm_mpz (factors, n, methods[i], fm, cfm,
                                      lpb, BB, BBB);

      if (fm->arith != CHOOSE_NONE)
        modset_clear (fm);
      if (cfm->arith != CHOOSE_NONE)
        modset_clear (cfm);

      if (res_fac == FACUL_NOT_SMOOTH) /* should not happen */
        {
          gmp_fprintf (stderr, "Error, non-smooth cofactor %Zd\n", n);
          exit (1);
        }

      for (j = 0; j < res_fac; j++)
        {
          if (m++ > 0)
            printf (",");
          gmp_printf ("%Zx", factors[j]);
          mpz_divexact (n, n, factors[j]);
        }
    }

  if (mpz_cmp_ui (n, 1) > 0)
    {
      if (mpz_cmp_d (n, BB) >= 0) /* by construction BB >= L */
        {
          gmp_fprintf (stderr, "Error, unfactored %Zd, please increase NB_MAX_METHODS\n", n);
          exit (1);
        }
      if (m++ > 0)
        printf (",");
      gmp_printf ("%Zx", n);
    }

}

/* Given a list L of bi-smooth cofactors, print the corresponding relations
   on stdout. */
static void
factor (cofac_list L, const char *poly_file, int lpb1, int lpb0)
{
  cado_poly pol;
  unsigned long n = L->size, i, pmax;
  int nb_methods;
  facul_method_t *methods;
#define MAX_FACTORS 16
  mpz_t factors[MAX_FACTORS];
  struct modset_t fm, cfm;
  double BB, BBB;
  mpz_list SP;
  prime_info pi;

  /* we trial divide by all primes < L^(1/2), so that any factor < L
     is necessarily prime */
  pmax = (lpb1 > lpb0) ? lpb1 : lpb0;
  pmax = (unsigned long) ceil (pow (2.0, (double) pmax / 2.0));
  mpz_list_init (SP);
  prime_info_init (pi);
  prime_list (SP, pi, 2, pmax);

  nb_methods = 50;
  if (nb_methods >= NB_MAX_METHODS)
    nb_methods = NB_MAX_METHODS - 1;
  methods = facul_make_default_strategy (nb_methods - 3, 0);
  cado_poly_init (pol);
  if (cado_poly_read (pol, poly_file) == 0)
    {
      fprintf (stderr, "Could not read polynomial file\n");
      exit (1);
    }

  for (i = 0; i < MAX_FACTORS; i++)
    mpz_init (factors[i]);
  BB = (double) pmax * (double) pmax;
  BBB = BB * (double) pmax;

  /* compute all norms F(a,b) and G(a,b) */
  for (i = 0; i < n; i++)
    {
      printf ("%" PRId64 ",%" PRIu64 ":", L->a[i], L->b[i]);

      /* at this point L->R[i] and L->A[i] contain the product of all prime
         factors > lim0 and lim1 respectively (apart from the special-q) */

      mpz_poly_homogeneous_eval_siui (L->R[i], pol->pols[0], L->a[i], L->b[i]);
      print_smooth (factors, L->R[i], methods, nb_methods,
                    &fm, &cfm, lpb0, BB, BBB, SP);
      printf (":");

      mpz_poly_homogeneous_eval_siui (L->A[i], pol->pols[1], L->a[i], L->b[i]);
      print_smooth (factors, L->A[i], methods, nb_methods,
                    &fm, &cfm, lpb1, BB, BBB, SP);
      printf ("\n");
      fflush (stdout);
    }

  for (i = 0; i < MAX_FACTORS; i++)
    mpz_clear (factors[i]);

  prime_info_clear (pi);
  mpz_list_clear (SP);
  cado_poly_clear (pol);
  facul_clear_aux_methods (methods);
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
  const char *poly_file, *cofac_file;
  int lpb1, lpb0, verbose, split = 0;
  unsigned long lim1, lim0;

  start = seconds ();

  param_list_init (pl);
  declare_usage (pl);

  param_list_configure_switch (pl, "-v", NULL);

  argc--, argv++;
  for( ; argc ; )
    {
      if (param_list_update_cmdline (pl, &argc, &argv))
        continue;
      else
        break;
    }

  poly_file = param_list_lookup_string (pl, "poly");
  cofac_file = param_list_lookup_string (pl, "cofac");
  ASSERT_ALWAYS(param_list_parse_int (pl, "lpb1", &lpb1));
  ASSERT_ALWAYS(param_list_parse_int (pl, "lpb0", &lpb0));
  ASSERT_ALWAYS(param_list_parse_ulong (pl, "lim1", &lim1));
  ASSERT_ALWAYS(param_list_parse_ulong (pl, "lim0", &lim0));
  param_list_parse_int (pl, "split", &split);
  verbose = param_list_parse_switch (pl, "-v");

  if (poly_file == NULL)
    {
      fprintf (stderr, "Error, missing -poly <file>\n");
      exit (1);
    }

  if (cofac_file == NULL)
    {
      fprintf (stderr, "Error, missing -cofac <file>\n");
      exit (1);
    }

  /* Initialization */
  cofac_list_init (L);

  cofac = fopen (cofac_file, "r");

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
  find_smooth (L, lpb1, lpb0, lim1, lim0, split, verbose);
  fprintf (stderr, "Detecting %lu smooth cofactors took %.0f s\n", L->size,
           seconds() - start);

  start = seconds ();
  factor (L, poly_file, lpb1, lpb0);
  fprintf (stderr, "Factoring %lu smooth cofactors took %.0f s\n", L->size,
           seconds() - start);

  mpz_clear (R);
  mpz_clear (A);
  cofac_list_clear (L);

  /* Clear the ecm addition chains that are stored as global variables
     to avoid re-computations. */
  free_saved_chains ();

  param_list_clear (pl);

  return 0;
}
