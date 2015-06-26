/* Bernstein's smoothness test */

/* This is a first very naive implementation. To use it:
   1) create a file (say cofac) containing lines of the following form:
      a b cofac_r cofac_a
      where cofac_r and cofac_a and cofactors on the rational and algebraic
      sides respectively
   2) run ./smooth cofac after updating the lines with "UPDATE".

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
#include "smooth.h"
#include "utils.h"
#include "portability.h"

#define NB_MILLER_RABIN 1

#define STATUS_SMOOTH  0
#define STATUS_UNKNOWN 1
#define STATUS_USELESS 2

void
mpz_list_init (mpz_list L)
{
  L->l = NULL;
  L->alloc = L->size = 0;
}

/* at any point, all values up to L->alloc should be mpz_init()'ed,
   even beyond L->size */
void
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

void
mpz_list_clear (mpz_list L)
{
  size_t i;

  for (i = 0; i < L->alloc; i++)
    mpz_clear (L->l[i]);
  free (L->l);
  L->alloc = L->size = 0;
}

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
   Each R[j] has been divided by its P-smooth part.
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

void
prime_product_init (prime_info pi, unsigned long p_max, unsigned long *p_last)
{
  unsigned long p;

  prime_info_init (pi);
  for (p = 2; p <= p_max; p = getprime_mt (pi));
  *p_last = p;
}

/* add in the list L all primes pmin <= p < pmax.
   Assume pmin is the current prime in 'pi'
   (pmin=2 when 'pi' was just initialized).
   Return the current prime in 'pi' at the end, i.e.,
   the smallest prime >= pmax. */
unsigned long
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
unsigned long
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
void
update_status (mpz_t *R, mpz_t *A,
               unsigned char *b_status_r, unsigned char *b_status_a,
               unsigned long *n, unsigned long *nb_rel_smooth,
               unsigned long rlim, unsigned long lpbr,
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

  mpz_set_ui (z_B2, rlim);
  mpz_mul_ui (z_B2, z_B2, rlim);
  mpz_mul_ui (z_B3, z_B2, rlim);
  mpz_setbit (z_L2, 2 * lpbr);
  B = rlim;
  L = 1UL << lpbr;
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

  /* Initial one-side pass (to make rlim and alim be equal) */

  unsigned long rsize, asize;
  /* the product of primes in [rlim, rlim + step] has a number of bits about:
     sum(log(p)/log(2)/log(p), p = rlim..rlim + step) ~ step/log(2)
     thus we need step ~ rsize*log(2) */
  rsize = estimate_bit_size (l->R, l->size);
  rlim_step = (unsigned long) ((double) rsize * log (2.0));
  asize = estimate_bit_size (l->A, l->size);
  alim_step = (unsigned long) ((double) asize * log (2.0));

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
      prime = prime_product (P, pi, rlim_new, prime);
      fprintf (stderr, "Computing prime product of %lu bits took %.0f seconds\n",
               mpz_sizeinbase (P, 2), seconds () - s);

      nb_rel = nb_rel_smooth;
      while (nb_rel < nb_rel_unknown)
      {
        nb_rel_new = (nb_rel + nb_rel_step < nb_rel_unknown ? nb_rel + nb_rel_step : nb_rel_unknown);
        s = seconds ();
        t_smooth -= seconds();
        smoothness_test (&(l->R[nb_rel]), nb_rel_new - nb_rel, P);
        t_smooth += seconds();
        fprintf (stderr, "smoothness test (%lu cofactors) took %.0f seconds"
                 " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
        nb_rel = nb_rel_new;
      }

      rlim = rlim_new;
      t_update -= seconds();
      update_status (l->R, l->A, b_status_r, b_status_a, &nb_rel_unknown,
                     &nb_rel_smooth, rlim, lpbr, &nb_smooth_r, &nb_smooth_a,
                     &nb_useless, l->a, l->b);
      t_update += seconds();
      fprintf (stderr, "smooth_r:%lu useless:%lu unknown:%lu rel_smooth:%lu\n",
               nb_smooth_r, nb_useless, nb_rel_read - nb_smooth_r - nb_useless, nb_rel_smooth);
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
      prime = prime_product (P, pi, alim_new, prime);
      fprintf (stderr, "Computing prime product of %lu bits took %.0f seconds\n",
               mpz_sizeinbase (P, 2), seconds () - s);

      nb_rel = nb_rel_smooth;
      while (nb_rel < nb_rel_unknown)
      {
        nb_rel_new = (nb_rel + nb_rel_step < nb_rel_unknown ? nb_rel + nb_rel_step : nb_rel_unknown);
        s = seconds ();
        t_smooth -= seconds();
        smoothness_test (&(l->A[nb_rel]), nb_rel_new - nb_rel, P);
        t_smooth += seconds();
        fprintf (stderr, "smoothness test (%lu cofactors) took %.0f seconds"
                 " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
        nb_rel = nb_rel_new;
      }

      alim = alim_new;
      t_update -= seconds();
      update_status (l->A, l->R, b_status_a, b_status_r, &nb_rel_unknown,
                     &nb_rel_smooth, alim, lpba, &nb_smooth_a, &nb_smooth_r,
                     &nb_useless, l->a, l->b);
      t_update += seconds();
      fprintf (stderr, "smooth_a:%lu useless:%lu unknown:%lu rel_smooth:%lu\n",
               nb_smooth_a, nb_useless, nb_rel_read - nb_smooth_a - nb_useless, nb_rel_smooth);
      fprintf (stderr, "t_update: %.0f seconds\n", t_update);
    }
  }
  else  // rlim = alim
  {
    for (prime = 2; prime <= rlim; prime = getprime_mt (pi));
  }

  ASSERT_ALWAYS (rlim == alim);

  /* the loop below assumes rlim_step = alim_step */
  if (rlim_step < alim_step)
    rlim_step = alim_step;
  else
    alim_step = rlim_step;

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
    prime = prime_product (P, pi, rlim_new, prime);
    fprintf (stderr, "Computing prime product of %lu bits took %.0f seconds\n",
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
                   &nb_useless, l->a, l->b);
    t_update += seconds();
    fprintf (stderr, "smooth_r:%lu useless:%lu unknown:%lu rel_smooth:%lu\n",
             nb_smooth_r, nb_useless, nb_rel_read - nb_smooth_r - nb_useless, nb_rel_smooth);
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
                   &nb_useless, l->a, l->b);
    t_update += seconds();
    fprintf (stderr, "smooth_a:%lu useless:%lu unknown:%lu rel_smooth:%lu\n",
             nb_smooth_a, nb_useless, nb_rel_read - nb_smooth_a - nb_useless, nb_rel_smooth);
    fprintf (stderr, "t_update: %.0f seconds\n", t_update);
  }
  prime_info_clear (pi);

  fprintf (stderr, "\nCollect smooth relations: %.0f s\n\n", seconds() - start);

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
compute_biproduct (mpz_t Q, cofac_list L, int alg)
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
#if 0
      mpz_mul (Q, L->R[i], L->A[i]);
      mpz_mul (T[0], T[0], Q);
#else
      if (alg)
        mpz_mul (T[0], T[0], L->A[i]);
      else
        mpz_mul (T[0], T[0], L->R[i]);
#endif
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
      static int count = 0;
      if (count++ < 10)
        gmp_fprintf (stderr, "Error, remaining unfactored part %Zd\n", N);
      return -1;
    }
  return len;
}

void
print_relations (cofac_list L, mpz_t *PR, unsigned long nr,
                 mpz_t *PA, unsigned long na,
                 int lpbr, int lpba)
{
  unsigned long i;
  int ret, len;
  double s = seconds ();
  char line[1024];
  unsigned long printed = 0;

  for (i = 0; i < L->size; i++)
    {
      len = sprintf (line, "%" PRId64 ":%" PRIu64 ":", L->a[i], L->b[i]);
      ret = print_norm (line + len, L->R[i], PR, nr, lpbr);
      if (ret >= 0)
        {
          len += ret;
          len += sprintf (line + len, ":");
          ret = print_norm (line + len, L->A[i], PA, na, lpba);
        }
      if (ret >= 0)
        {
          printf ("%s\n", line);
          printed ++;
        }
      else
        {
          static int count = 0;
          if (count++ < 10)
            fprintf (stderr, "Error for %" PRId64 ":%" PRIu64 "\n",
                     L->a[i], L->b[i]);
        }
    }
  fprintf (stderr, "Printed %lu relations in %.0f s\n",
           printed, seconds () - s);
}

/* put in LP[0..i-1] the good primes, and return i */
unsigned long
scan_primes (mpz_t *LP, mpz_t *T, unsigned long nprimes)
{
  mpz_t t;
  unsigned long i, j;

  mpz_init (t);
  for (i = j = 0; j < nprimes / 2; j++)
    {
      mpz_mod (t, T[j], LP[2*j]);
      if (mpz_cmp_ui (t, 0) == 0)
        mpz_swap (LP[i++], LP[2*j]);
      mpz_mod (t, T[j], LP[2*j+1]);
      if (mpz_cmp_ui (t, 0) == 0)
        mpz_swap (LP[i++], LP[2*j+1]);
    }
  if (nprimes & 1)
    {
      mpz_mod (t, T[j], LP[2*j]);
      if (mpz_cmp_ui (t, 0) == 0)
        mpz_swap (LP[i++], LP[2*j]);
    }
  mpz_clear (t);
  return i;
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

Note: we should use instead the following algorithm from
"How to find smooth parts of integers", Dan Bernstein, draft, 2004,
http://cr.yp.to/factorization/smoothparts-20040510.pdf.
Given cofactors r_1, ..., r_n and primes p_1, ..., p_k:
(1) compute R = r_1 * ... * r_n
(2) compute R mod p_1, ..., R mod p_k and keep only those p_j for which
    R mod p_j = 0
(3) if n >= 2, cut r_1, ..., r_n in two parts, and do the same recursively for
    the first and second parts
*/
void
factor (cofac_list L, const char *poly_file, int lpba, int lpbr)
{
  cado_poly pol;
  unsigned long n = L->size, i, nprimes, w[MAX_DEPTH];
  mpz_t QR, QA, t;
  double s, s_product = 0, s_remainder = 0;

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
  mpz_init (QR);
  mpz_init (QA);
  fprintf (stderr, "factor: computing product of all norms...");
  fflush (stderr);
  s = seconds ();
  compute_biproduct (QR, L, 0);
  compute_biproduct (QA, L, 1);
  fprintf (stderr, "done in %.0f s\n", seconds () - s);
  fflush (stderr);

  /* compute all primes up to max(2^lpba, 2^lpbr) */
  lpba = (lpba > lpbr) ? lpba : lpbr;
  unsigned long Qbits;
  prime_info pi;
  unsigned long pmin = 2, pmax = 1UL << lpba;
  mpz_list PR, PA;
  mpz_t *LPR, *LPA;

  Qbits = mpz_cmp (QA, QR) > 0 ? mpz_sizeinbase (QA, 2)
    : mpz_sizeinbase (QR, 2);

  mpz_list_init (PR);
  mpz_list_init (PA);
  prime_info_init (pi);
  mpz_init (t);
  while (pmin < pmax)
    {
      /* we consider primes in [pmin, pmin + step], whose product has about
         step/log(2) bits, thus we want step ~ Qbits * log(2) */
      unsigned long phigh = pmin + (unsigned long) ((double) Qbits * log (2.0));
      if (phigh > pmax)
        phigh = pmax;
      nprimes = PR->size;
      pmin = prime_list (PR, pi, pmin, phigh);
      nprimes = PR->size - nprimes; /* number of new primes */
      LPR = PR->l + (PR->size - nprimes); /* start location of new primes */

      /* copy new primes in PA */
      for (i = 0; i < nprimes; i++)
        mpz_list_add (PA, mpz_get_ui (PR->l[PR->size - nprimes + i]));
      LPA = PA->l + (PA->size - nprimes);

      /* form a product tree from LPR */
      mpz_t **T;
      s_product -= seconds ();
      T = product_tree (LPR, nprimes, w);
      s_product += seconds ();

      /* compute the remainder tree Q mod T */
      s_remainder -= seconds ();
      remainder_tree (T, nprimes, w, QR);
      /* now scan all primes appearing in norms */
      i = scan_primes (LPR, T[1], nprimes);
      s_remainder += seconds ();

      clear_product_tree (T, nprimes, w);
      /* the primes appearing in relations are LP[0..i-1] */
      PR->size = (LPR + i) - PR->l; /* reajust the 'appearing' primes */

      /* same thing for algebraic side */
      s_product -= seconds ();
      T = product_tree (LPA, nprimes, w);
      s_product += seconds ();

      /* compute the remainder tree Q mod T */
      s_remainder -= seconds ();
      remainder_tree (T, nprimes, w, QA);
      /* now scan all primes appearing in norms */
      i = scan_primes (LPA, T[1], nprimes);
      s_remainder += seconds ();

      clear_product_tree (T, nprimes, w);
      PA->size = (LPA + i) - PA->l;
    }
  fprintf (stderr, "factor: computing product tree of primes took %1.0f s\n",
           s_product);
  fprintf (stderr, "factor: computing remainder tree took %1.0f s\n",
           s_remainder);
  fflush (stderr);
  mpz_clear (t);
  prime_info_clear (pi);
  mpz_clear (QR);
  mpz_clear (QA);
  cado_poly_clear (pol);

  print_relations (L, PR->l, PR->size, PA->l, PA->size, lpbr, lpba);

  mpz_list_clear (PR);
  mpz_list_clear (PA);
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
