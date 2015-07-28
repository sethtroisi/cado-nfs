#include "cado.h"
#include <stdio.h>
#include "batch.h"
#include "utils.h"

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

/* same as prime_list, but only adds primes p for which f has at least one
   root modulo p */
static unsigned long
prime_list_poly (mpz_list L, prime_info pi, unsigned long pmin,
                 unsigned long pmax, mpz_poly_t f)
{
  unsigned long p;

  if (f->deg == 1)
    return prime_list (L, pi, pmin, pmax);

  for (p = pmin; p < pmax; p = getprime_mt (pi))
    if (mpz_poly_roots_ulong (NULL, f, p) > 0)
      mpz_list_add (L, p);
  return p;
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
  mpz_set (l->R[l->size], R);
  mpz_set (l->A[l->size], A);
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

/* same as prime_product, but keeps only primes for which the given polynomial
   has factors modulo p */
unsigned long
prime_product_poly (mpz_t P, prime_info pi, unsigned long p_max,
                    unsigned long p_last, mpz_poly_t f)
{
  unsigned long i;
  mpz_list L;

  mpz_list_init (L);
  p_last = prime_list_poly (L, pi, p_last, p_max, f);

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

#if 0
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
#else
/* Compute the remainder_tree using the "scaled" variant
   (http://cr.yp.to/arith/scaledmod-20040820.pdf).
   At the root, we compute a floating-point
   approximation of P/T[h][0] with m+guard bits, where m = nbits(T[h][0]). */
static void
remainder_tree (mpz_t **T, unsigned long n, unsigned long *w, mpz_t P)
{
  unsigned long h = tree_height (n), i, j, guard;
  unsigned long **nbits;
  mpz_t Q;

  guard = h;
  nbits = malloc ((h + 1) * sizeof (unsigned long*));
  for (i = 0; i <= h; i++)
    {
      nbits[i] = malloc (w[i] * sizeof (unsigned long));
      for (j = 0; j < w[i]; j++)
        nbits[i][j] = mpz_sizeinbase (T[i][j], 2);
    }

  mpz_init (Q);
  mpz_mod (Q, P, T[h][0]); /* first reduce modulo T[h][0] in case P is huge */
  mpz_mul_2exp (Q, Q, nbits[h][0] + guard);
  mpz_tdiv_q (T[h][0], Q, T[h][0]);
  /* P/T[h][0] ~ Q/2^(m+guard) */
  for (i = h; i > 1; i--)
    {
      for (j = 0; j < w[i-1] / 2; j++)
        {
          /* T[i][j]/2^(nbits[i][j] + guard) ~ P/T[i][j] */
          mpz_mul (T[i-1][2*j], T[i][j], T[i-1][2*j]);
          /* same for the right part */
          mpz_mul (T[i-1][2*j+1], T[i][j], T[i-1][2*j+1]);
          /* swap */
          mpz_swap (T[i-1][2*j], T[i-1][2*j+1]);

          /* get the fractional part, i.e., the low nbits[i][j] + guard bits */
          mpz_tdiv_r_2exp (T[i-1][2*j], T[i-1][2*j], nbits[i][j] + guard);
          /* now keep only nbits[i-1][2*j] + guard significant bits */
          mpz_div_2exp (T[i-1][2*j], T[i-1][2*j], nbits[i][j] - nbits[i-1][2*j]);

          mpz_tdiv_r_2exp (T[i-1][2*j+1], T[i-1][2*j+1], nbits[i][j] + guard);
          mpz_div_2exp (T[i-1][2*j+1], T[i-1][2*j+1], nbits[i][j] - nbits[i-1][2*j+1]);

        }
      if (w[i-1] & 1)
        mpz_swap (T[i-1][w[i-1]-1], T[i][w[i]-1]);
    }

  /* from X[1][j] ~ P/T[1][j]*2^(nbits[1][j] + guard) mod 2^(nbits[1][j] + guard),
     get X[1][j]*T[1][j]/2^(nbits[1][j] + guard) ~ P mod T[1][j] */
  for (j = 0; j < w[1]; j++)
    {
      mpz_mul (T[1][j], T[1][j], T[0][2*j]);
      if (2*j+1 < w[0])
        mpz_mul (T[1][j], T[1][j], T[0][2*j+1]);
      /* now X[1][j] ~ P*2^(nbits[1][j] + guard) */
      mpz_div_2exp (T[1][j], T[1][j], nbits[1][j]);
      /* now X[1][j] ~ P*2^guard */
      mpz_add_ui (T[1][j], T[1][j], (1UL << guard) - 1UL);
      mpz_div_2exp (T[1][j], T[1][j], guard);
    }

  mpz_clear (Q);
  for (i = 0; i <= h; i++)
    free (nbits[i]);
  free (nbits);
}
#endif

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
  static double t_rem = 0;

  mpz_init (w1);
  mpz_init (w2);

  /* special case for n=1 */
  if (n == 1)
    {
      mpz_gcd (w1, R[0], P);
      mpz_divexact (R[0], R[0], w1);
      goto clear_and_exit;
    }

  T = product_tree (R, n, w);

  if (verbose > 1)
    fprintf (stderr, "# cofactor product has %zu bits\n",
             mpz_sizeinbase (T[h][0], 2));

  /* compute remainder tree */
  t_rem -= seconds ();
  remainder_tree (T, n, w, P);
  t_rem += seconds ();
  if (verbose > 1)
    fprintf (stderr, "# remainder_tree: %.1fs\n", t_rem);

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

  clear_product_tree (T, n, w);
 clear_and_exit:
  mpz_clear (w1);
  mpz_clear (w2);
}

#define NB_MILLER_RABIN 1

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

/* return the number n of smooth relations in l,
   which should be at the end in locations 0, 1, ..., n-1 */
void
find_smooth (cofac_list l, int lpb0, int lpb1,
             unsigned long lim0, unsigned long lim1,
             FILE *batch0, FILE *batch1, int verbose)
{
  unsigned long nb_rel;
  unsigned long nb_rel_new;
  unsigned long nb_rel_read = l->size;
  unsigned long nb_rel_unknown;
  unsigned long i;
  mpz_t P;
  double s;
  double t_smooth = 0;
  double t_update = 0;
  double t_prime = 0;
  double start;
  unsigned int n0_pass;
  unsigned long lim0_new;
  unsigned long lim1_new;
  unsigned char *b_status_r;
  unsigned char *b_status_a;
  unsigned long nb_rel_smooth;
  unsigned long nb_smooth_r;
  unsigned long nb_smooth_a;
  unsigned long nb_useless;
  int ret;

  ASSERT_ALWAYS(batch0 != NULL);
  ASSERT_ALWAYS(batch1 != NULL);

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

  /* the code below assumes max(lim0,lim1) <= min(2^lpb0,2^lpb1) */
  ASSERT_ALWAYS(lim0 <= (1UL << lpb0));
  ASSERT_ALWAYS(lim0 <= (1UL << lpb1));
  ASSERT_ALWAYS(lim1 <= (1UL << lpb0));
  ASSERT_ALWAYS(lim1 <= (1UL << lpb1));

  /* Loop */

  n0_pass = 0;
  while ( (lim0 < (1UL << lpb0)) || (lim1 < (1UL << lpb1)) )
  {
    n0_pass++;

    if (verbose > 1)
      fprintf (stderr, "#\n# Pass %u: %.0f s\n", n0_pass, seconds() - start);

    /* side 0 */

    s = seconds ();
    ret = gmp_fscanf (batch0, "%lu %Zx\n", &lim0_new, P);
    ASSERT_ALWAYS(ret == 2);
    s = seconds () - s;
    t_prime += s;
    if (verbose > 1)
      fprintf (stderr, "# Reading prime product of %zu bits took %.0fs (total %.0fs so far)\n",
               mpz_sizeinbase (P, 2), s, t_prime);

    nb_rel = nb_rel_smooth;
    nb_rel_new = nb_rel_unknown;
    s = seconds ();
    t_smooth -= seconds();
    smoothness_test (&(l->R[nb_rel]), nb_rel_new - nb_rel, P, verbose);
    t_smooth += seconds();
    if (verbose > 1)
      fprintf (stderr, "# smoothness_test (%lu cofactors) took %.0f seconds"
               " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
    nb_rel = nb_rel_new;

    lim0 = lim0_new;
    t_update -= seconds();
    update_status (l->R, l->A, b_status_r, b_status_a, &nb_rel_unknown,
                   &nb_rel_smooth, lim0, lpb0, &nb_smooth_r, &nb_smooth_a,
                   &nb_useless, l->a, l->b);
    t_update += seconds();
    if (verbose > 1)
      {
        fprintf (stderr, "# smooth_r:%lu useless:%lu unknown:%lu rel_smooth:%lu\n",
                 nb_smooth_r, nb_useless, nb_rel_read - nb_smooth_r - nb_useless, nb_rel_smooth);
        fprintf (stderr, "# t_update: %.0f seconds\n", t_update);
      }

    /* side 1 */

    s = seconds ();
    ret = gmp_fscanf (batch1, "%lu %Zx\n", &lim1_new, P);
    ASSERT_ALWAYS(ret == 2);
    s = seconds () - s;
    t_prime += s;
    if (verbose > 1)
      fprintf (stderr, "# Reading prime product of %zu bits took %.0fs (total %.0fs so far)\n",
               mpz_sizeinbase (P, 2), s, t_prime);

    nb_rel = nb_rel_smooth;
    nb_rel_new = nb_rel_unknown;
    s = seconds ();
    t_smooth -= seconds();
    smoothness_test (&(l->A[nb_rel]), nb_rel_new - nb_rel, P, verbose);
    t_smooth += seconds();
    if (verbose > 1)
      fprintf (stderr, "# smoothness_test (%lu cofactors) took %.0f seconds"
               " (total %.0f so far)\n", nb_rel_new - nb_rel, seconds () - s, t_smooth);
    nb_rel = nb_rel_new;

    lim1 = lim1_new;
    t_update -= seconds();
    update_status (l->A, l->R, b_status_a, b_status_r, &nb_rel_unknown,
                   &nb_rel_smooth, lim1, lpb1, &nb_smooth_a, &nb_smooth_r,
                   &nb_useless, l->a, l->b);
    t_update += seconds();
    if (verbose > 1)
      {
        fprintf (stderr, "# smooth_a:%lu useless:%lu unknown:%lu rel_smooth:%lu\n",
                 nb_smooth_a, nb_useless, nb_rel_read - nb_smooth_a - nb_useless, nb_rel_smooth);
        fprintf (stderr, "# t_update: %.0f seconds\n", t_update);
      }
  }
  cofac_list_realloc (l, nb_rel_smooth);

  mpz_clear (P);
  free (b_status_r);
  free (b_status_a);

  if (verbose)
    {
      fprintf (stderr, "# find_smooth %.1fs (t_prime %.1fs, t_smooth %.1fs, t_update %.1fs)\n",
               seconds () - start, t_prime, t_smooth, t_update);
      fprintf (stderr, "# out of %lu rels, found %lu smooth\n",
               nb_rel_read, nb_rel_smooth);
    }
}

/* return the new value of m */
static int
print_smooth_aux (mpz_t *factors, mpz_t n, facul_method_t *methods,
                  struct modset_t *fm, struct modset_t *cfm,
                  int lpb, double BB, double BBB, int m)
{
  unsigned long i;
  int j, res_fac;

  for (i = 0; methods[i].method != 0 && mpz_cmp_d (n, BB) >= 0; i++)
    {
      res_fac = facul_doit_onefm_mpz (factors, n, methods[i], fm, cfm,
                                      lpb, BB, BBB);

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

      if (fm->arith != CHOOSE_NONE) /* fm is a non-trivial composite factor */
        {
          mpz_t t;
          struct modset_t cfm2;
          mpz_init (t);
          modset_get_z (t, fm);
          mpz_divexact (n, n, t);
          m = print_smooth_aux (factors, t, methods + i + 1, fm, &cfm2,
                                lpb, BB, BBB, m);
          modset_clear (fm);
          fm->arith = CHOOSE_NONE;
          mpz_clear (t);
        }

      if (cfm->arith != CHOOSE_NONE)
        {
          mpz_t t;
          struct modset_t fm2;
          mpz_init (t);
          modset_get_z (t, cfm);
          mpz_divexact (n, n, t);
          m = print_smooth_aux (factors, t, methods + i + 1, &fm2, cfm,
                                lpb, BB, BBB, m);
          modset_clear (cfm);
          cfm->arith = CHOOSE_NONE;
          mpz_clear (t);
        }
    }

  return m;
}

/* Print the prime factor of the input 'n' separated by spaces.
   The list SP (small primes) contains all primes < B.
   m is the number of already printed factors.
   Return 0 in case of error.
*/
static int
print_smooth (mpz_t *factors, mpz_t n, facul_method_t *methods,
              struct modset_t *fm, struct modset_t *cfm,
              int lpb, double BB, double BBB, mpz_list SP, int m)
{
  int res_fac, j;
  unsigned long i;

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

  for (i = 0; methods[i].method != 0 && mpz_cmp_d (n, BB) >= 0; i++)
    {
      res_fac = facul_doit_onefm_mpz (factors, n, methods[i], fm, cfm,
                                      lpb, BB, BBB);

      if (res_fac == FACUL_NOT_SMOOTH) /* should not happen */
        {
          gmp_fprintf (stderr, "Error, non-smooth cofactor %Zd\n", n);
          return 0;
        }

      for (j = 0; j < res_fac; j++)
        {
          if (m++ > 0)
            printf (",");
          gmp_printf ("%Zx", factors[j]);
          mpz_divexact (n, n, factors[j]);
        }

      if (fm->arith != CHOOSE_NONE) /* fm is a non-trivial composite factor */
        {
          mpz_t t;
          struct modset_t cfm2;
          mpz_init (t);
          modset_get_z (t, fm);
          mpz_divexact (n, n, t);
          m = print_smooth_aux (factors, t, methods + i + 1, fm, &cfm2,
                                lpb, BB, BBB, m);
          modset_clear (fm);
          fm->arith = CHOOSE_NONE;
          mpz_clear (t);
        }

      if (cfm->arith != CHOOSE_NONE)
        {
          mpz_t t;
          struct modset_t fm2;
          mpz_init (t);
          modset_get_z (t, cfm);
          mpz_divexact (n, n, t);
          m = print_smooth_aux (factors, t, methods + i + 1, &fm2, cfm,
                                lpb, BB, BBB, m);
          modset_clear (cfm);
          cfm->arith = CHOOSE_NONE;
          mpz_clear (t);
        }

    }

  if (mpz_cmp_ui (n, 1) > 0)
    {
      if (mpz_cmp_d (n, BB) >= 0) /* by construction BB >= L */
        {
          gmp_fprintf (stderr, "Error, unfactored %Zd, please increase NB_MAX_METHODS\n", n);
          return 0;
        }
      if (m++ > 0)
        printf (",");
      gmp_printf ("%Zx", n);
    }
  return 1;
}

/* Given a list L of bi-smooth cofactors, print the corresponding relations
   on stdout. */
void
factor (cofac_list L, cado_poly pol, int lpb0, int lpb1, int verbose)
{
  unsigned long n = L->size, i, pmax;
  int nb_methods;
  facul_method_t *methods;
#define MAX_FACTORS 16
  mpz_t factors[MAX_FACTORS];
  struct modset_t fm, cfm;
  double BB, BBB;
  mpz_list SP;
  prime_info pi;
  double start = seconds ();

  /* we trial divide by all primes < L^(1/2), so that any factor < L
     is necessarily prime */
  pmax = (lpb1 > lpb0) ? lpb1 : lpb0;
  pmax = (unsigned long) ceil (pow (2.0, (double) pmax / 2.0));
  mpz_list_init (SP);
  prime_info_init (pi);
  pmax = 4 * pmax;
  prime_list (SP, pi, 2, pmax);

  nb_methods = 30;
  if (nb_methods >= NB_MAX_METHODS)
    nb_methods = NB_MAX_METHODS - 1;
  methods = facul_make_default_strategy (nb_methods - 3, 0);

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
      if (print_smooth (factors, L->R[i], methods,
                        &fm, &cfm, lpb0, BB, BBB, SP, 0) == 0)
        {
          gmp_fprintf (stderr, "Error for a=%" PRId64 " b=%" PRIu64 "\n",
                       L->a[i], L->b[i]);
          fflush (stderr);
          exit (1);
        }
      printf (":");

      mpz_poly_homogeneous_eval_siui (L->A[i], pol->pols[1], L->a[i], L->b[i]);
      if (print_smooth (factors, L->A[i], methods,
                        &fm, &cfm, lpb1, BB, BBB, SP, 0) == 0)
        {
          gmp_fprintf (stderr, "Error for a=%" PRId64 " b=%" PRIu64 "\n",
                       L->a[i], L->b[i]);
          fflush (stderr);
          exit (1);
        }
      printf ("\n");
      fflush (stdout);
    }

  for (i = 0; i < MAX_FACTORS; i++)
    mpz_clear (factors[i]);

  prime_info_clear (pi);
  mpz_list_clear (SP);
  facul_clear_aux_methods (methods);

  if (verbose)
    fprintf (stderr, "# factor: %.1fs to print %lu rels\n",
             seconds () - start, n);
}

void
create_batch_file (const char *f, unsigned long B, unsigned long L,
                   mpz_poly_t pol, int split)
{
  FILE *fp;
  prime_info pi;
  unsigned long p, h;
  mpz_t P;
  double s = seconds ();

  ASSERT_ALWAYS (L > B);

  fprintf (stderr, "# creating batch file %s", f);
  fflush (stderr);

  prime_info_init (pi);
  fp = fopen (f, "w");
  ASSERT_ALWAYS(fp != NULL);
  mpz_init (P);

  h = (L - B + split - 1) / split;

  for (p = 2; p < B; p = getprime_mt (pi));

  while (B < L)
    {
      p = prime_product_poly (P, pi, (B + h < L) ? B + h : L, p, pol);
      gmp_fprintf (fp, "%lu %Zx\n", B + h, P);
      B += h;
    }

  fclose (fp);
  prime_info_clear (pi);
  mpz_clear (P);
  fprintf (stderr, " took %.0fs\n", seconds () - s);
  fflush (stderr);
}

