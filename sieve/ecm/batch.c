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
  ASSERT(L->size < L->alloc);
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
  l->perm = NULL;
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
  l->perm = realloc (l->perm, newsize * sizeof (uint32_t));
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
  l->perm[l->size] = l->size;
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
  free (l->perm);
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
product_tree (mpz_t *R, uint32_t *perm, unsigned long n, unsigned long *w)
{
  unsigned long h = tree_height (n), i, j;
  mpz_t **T;

  T = (mpz_t**) malloc ((h + 1) * sizeof (mpz_t*));

  /* initialize tree */
  w[0] = n;
  for (i = 0; i <= h; i++)
    {
      w[i] = 1 + ((n - 1) >> i);
      T[i] = (mpz_t*) malloc (w[i] * sizeof (mpz_t));
      for (j = 0; j < w[i]; j++)
        mpz_init (T[i][j]);
    }

  /* initialize T[0] to R */
  for (j = 0; j < n; j++)
    mpz_set (T[0][j], R[perm[j]]);

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

/* Compute the remainder tree using the "scaled" variant
   (http://cr.yp.to/arith/scaledmod-20040820.pdf).
   At the root, we compute a floating-point
   approximation of P/T[h][0] with m+guard bits, where m = nbits(T[h][0]). */
static void
remainder_tree (mpz_t **T, unsigned long n, unsigned long *w, mpz_t P,
                mpz_t *R, uint32_t *perm)
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
  for (i = h; i > 0; i--)
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

  /* from T[0][j] ~ P/R[j]*2^(nbits[0][j] + guard) mod 2^(nbits[0][j] + guard),
     get T[0][j]*R[j]/2^(nbits[0][j] + guard) ~ P mod R[j] */
  for (j = 0; j < n; j++)
    {
      mpz_mul (T[0][j], T[0][j], R[perm[j]]);
      /* T[0][j] ~ P*2^(nbits[0][j] + guard) mod R[j]*2^(nbits[0][j]+guard) */
      mpz_div_2exp (T[0][j], T[0][j], nbits[0][j]);
      /* T[0][j] ~ P*2^guard mod R[j]*2^guard */
      mpz_add_ui (T[0][j], T[0][j], (1UL << guard) - 1UL);
      mpz_div_2exp (T[0][j], T[0][j], guard);
    }

  mpz_clear (Q);
  for (i = 0; i <= h; i++)
    free (nbits[i]);
  free (nbits);
}

/* Clear the product tree. */
static void
clear_product_tree (mpz_t **T, unsigned long n, unsigned long *w)
{
  unsigned long i, j, h = tree_height (n);

  for (i = 0; i <= h; i++)
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
smoothness_test (mpz_t *R, uint32_t *perm, unsigned long n, mpz_t P,
                 int verbose)
{
  unsigned long h = tree_height (n), j, w[MAX_DEPTH];
  mpz_t **T;
  static double t_rem = 0;

  T = product_tree (R, perm, n, w);

  if (verbose > 1)
    fprintf (stderr, "# cofactor product has %zu bits\n",
             mpz_sizeinbase (T[h][0], 2));

  /* compute remainder tree */
  t_rem -= seconds ();
  remainder_tree (T, n, w, P, R, perm);
  t_rem += seconds ();
  if (verbose > 1)
    fprintf (stderr, "# remainder_tree: %.1fs\n", t_rem);

  /* now T[0][j] = P mod R[j] for 0 <= j < n */
  for (j = 0; j < n; j++)
    {
      mpz_gcd (T[0][j], T[0][j], R[perm[j]]);
      mpz_divexact (R[perm[j]], R[perm[j]], T[0][j]);
    }

  clear_product_tree (T, n, w);
}

#define NB_MILLER_RABIN 1

/* invariant:
   relations 0 to *nb_smooth-1 are smooth
   relations *nb_smooth to *nb_unknown-1 are unknown
   relations >= *nb_unknown are non-smooth (useless)
*/
static void
update_status (mpz_t *R, uint32_t *perm,
               unsigned char *b_status_r, unsigned char *b_status_a,
               unsigned long *nb_smooth, unsigned long *nb_unknown,
               unsigned long lim0, unsigned long lpb0)
{
  mpz_t z_B2, z_BL, z_B3, z_L2;

  unsigned long i, j;
  unsigned long B;
  unsigned long L;

  mpz_init (z_B2);
  mpz_init (z_BL);
  mpz_init (z_B3);
  mpz_init (z_L2); /* set to 0 */

  mpz_set_ui (z_B2, lim0);
  mpz_mul_ui (z_B2, z_B2, lim0);
  mpz_mul_ui (z_B3, z_B2, lim0);
  mpz_setbit (z_L2, 2 * lpb0);
  B = lim0;
  L = 1UL << lpb0;
  mpz_set_ui (z_BL, B);
  mpz_mul_ui (z_BL, z_BL, L);
  ASSERT_ALWAYS(mpz_cmp_ui (z_B2, L) >= 0);

  for (j = *nb_smooth; j < *nb_unknown; j++)
    {
      i = perm[j];
      if (b_status_r[i] == STATUS_UNKNOWN)
      {
        /* if R[i] < L, then R[i] is smooth (we assume L <= B^2) */
        if (mpz_cmp_ui (R[i], L) <= 0)
          goto smooth;
        /* if L^2 < R[i] < B^3 or L < R[i] < B^2, then R[i] cannot be smooth */
        else if ((0 < mpz_cmp (R[i], z_L2) && mpz_cmp (R[i], z_B3) < 0) ||
                 (0 < mpz_cmp_ui (R[i], L) && mpz_cmp (R[i], z_B2) < 0) ||
                 mpz_probab_prime_p (R[i], NB_MILLER_RABIN))
        {
          /* relation j is useless, swap it with relation *nb_unknown - 1 */
          (*nb_unknown)--;
          perm[j] = perm[*nb_unknown];
          perm[*nb_unknown] = i;
          j--;
        }
        /* now B^2 <= R[i] <= L^2 or B^3 <= R[i] and R[i] is composite */
        else if (mpz_cmp (R[i], z_BL) <= 0)
        {
        smooth:
          b_status_r[i] = STATUS_SMOOTH;
          if (b_status_a[i] == STATUS_SMOOTH)
          {
            /* relation j is smooth, swap it with relation *nb_smooth */
            perm[j] = perm[*nb_smooth];
            perm[*nb_smooth] = i;
            (*nb_smooth)++;
          }
        }
      }
    }

  mpz_clear (z_B2);
  mpz_clear (z_BL);
  mpz_clear (z_B3);
  mpz_clear (z_L2);
}

/* return the number n of smooth relations in l,
   which should be at the end in locations perm[0], perm[1], ..., perm[n-1] */
unsigned long
find_smooth (cofac_list l, int lpb[2], unsigned long lim[2],
             FILE *batch[2], int verbose)
{
  unsigned long nb_rel_read = l->size;
  unsigned long nb_smooth;
  unsigned long nb_unknown;
  unsigned long i;
  mpz_t P;
  double s;
  double t_smooth = 0;
  double t_update = 0;
  double t_prime = 0;
  double start;
  unsigned int n0_pass;
  unsigned long lim_new[2];
  unsigned char *b_status_r;
  unsigned char *b_status_a;
  int ret;

  ASSERT_ALWAYS(batch[0] != NULL);
  ASSERT_ALWAYS(batch[1] != NULL);

  start = seconds ();

  mpz_init (P);

  b_status_r = (unsigned char *) malloc (nb_rel_read * sizeof(unsigned char));
  b_status_a = (unsigned char *) malloc (nb_rel_read * sizeof(unsigned char));
  for (i = 0; i < nb_rel_read; i++)
    {
      b_status_r[i] = STATUS_UNKNOWN;
      b_status_a[i] = STATUS_UNKNOWN;
      ASSERT(mpz_cmp_ui (l->R[i], 0) > 0);
      ASSERT(mpz_cmp_ui (l->A[i], 0) > 0);
    }

  nb_smooth = 0;
  nb_unknown = nb_rel_read;

  /* the code below assumes lim0 <= 2^lpb0 and lim1 <= 2^lpb1 */
  ASSERT_ALWAYS(lim[0] <= (1UL << lpb[0]));
  ASSERT_ALWAYS(lim[1] <= (1UL << lpb[1]));

  /* Loop */

  /* invariant: the smooth relations are in 0..nb_smooth-1,
     the unknown ones in nb_smooth..nb_unknown-1,
     the remaining ones are not smooth */

  n0_pass = 0;
  while ( (lim[0] < (1UL << lpb[0])) || (lim[1] < (1UL << lpb[1])) )
  {
    n0_pass++;

    if (verbose)
      fprintf (stderr, "# Starting pass %u at %.1fs\n",
               n0_pass, seconds() - start);

    /* it seems faster to start from the algebraic side */
    for (int z = 1; z >= 0; z--)
      {
        s = seconds ();
        ret = gmp_fscanf (batch[z], "%lu %Zx\n", &(lim_new[z]), P);
        ASSERT_ALWAYS(ret == 2);
        s = seconds () - s;
        t_prime += s;
        if (verbose > 1)
          fprintf (stderr, "# Reading prime product of %zu bits took %.0fs (total %.0fs so far)\n",
                   mpz_sizeinbase (P, 2), s, t_prime);

        s = seconds ();
        t_smooth -= seconds();
        if (z == 0)
          smoothness_test (l->R, l->perm + nb_smooth, nb_unknown - nb_smooth, P, verbose);
        else
          smoothness_test (l->A, l->perm + nb_smooth, nb_unknown - nb_smooth, P, verbose);
        t_smooth += seconds();
        if (verbose > 1)
          fprintf (stderr, "# smoothness_test (%lu cofactors) took %.0f seconds"
                   " (total %.0f so far)\n", nb_unknown - nb_smooth, seconds () - s, t_smooth);

        lim[z] = lim_new[z];
        t_update -= seconds();
        /* we only need to update relations in [nb_smooth, nb_unknown-1] */
        if (z == 0)
          update_status (l->R, l->perm, b_status_r, b_status_a,
                         &nb_smooth, &nb_unknown, lim[z], lpb[z]);
        else
          update_status (l->A, l->perm, b_status_a, b_status_r,
                         &nb_smooth, &nb_unknown, lim[z], lpb[z]);
        t_update += seconds();
        if (verbose)
          fprintf (stderr, "# rel_smooth: %lu t_update: %.1f seconds\n",
                   nb_smooth, t_update);
      }
  }

  mpz_clear (P);
  free (b_status_r);
  free (b_status_a);

  if (verbose)
    {
      fprintf (stderr, "# find_smooth %.1fs (t_prime %.1fs, t_smooth %.1fs, t_update %.1fs)\n",
               seconds () - start, t_prime, t_smooth, t_update);
      fprintf (stderr, "# out of %lu rels, found %lu smooth\n",
               nb_rel_read, nb_smooth);
    }

  return nb_smooth;
}

/* print all factors of n, and return the new value of m */
static int
print_smooth_aux (mpz_t *factors, mpz_t n, facul_method_t *methods,
                  struct modset_t *fm, struct modset_t *cfm,
                  int lpb, double BB, double BBB, int m)
{
  unsigned long i;
  int j, res_fac;

  /* any factor < B^2 is necessarily prime */
  for (i = 0; methods[i].method != 0 && mpz_cmp_d (n, BB) >= 0; i++)
    {
      res_fac = facul_doit_onefm_mpz (factors, n, methods[i], fm, cfm,
                                      lpb, BB, BBB);

      ASSERT_ALWAYS(res_fac != FACUL_NOT_SMOOTH);

      /* factors[0..res_fac-1] are prime factors of n, 0 <= res_fac <= 2 */

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
          /* t should be composite, i.e., t >= BB */
          ASSERT(mpz_cmp_d (t, BB) >= 0);
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
          /* t should be composite, i.e., t >= BB */
          ASSERT(mpz_cmp_d (t, BB) >= 0);
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
      ASSERT_ALWAYS (mpz_cmp_d (n, BB) < 0);
      if (m++ > 0)
        printf (",");
      gmp_printf ("%Zx", n);
    }

  return m;
}

double td_time = 0;

static int
trial_divide (mpz_t n, unsigned long *sp, unsigned long spsize)
{
  int m = 0; /* number of already printed factors */
  unsigned long i;

  for (i = 0; i < spsize; i++)
    {
      while (mpz_divisible_ui_p (n, sp[i]))
        {
          if (m++ > 0)
            printf (",");
          printf ("%lx", sp[i]);
          mpz_divexact_ui (n, n, sp[i]);
        }
    }

  return m;
}

/* Print the prime factor of the input 'n' separated by spaces.
   The list SP (small primes) contains all primes < B.
   BB is the prime bound: any factor < BB is necessarily prime.
   'hint' is either 1 or a product of large primes.
*/
static void
print_smooth (mpz_t *factors, mpz_t n, facul_method_t *methods,
              struct modset_t *fm, struct modset_t *cfm,
              unsigned int lpb, double BB, double BBB, unsigned long *sp,
              unsigned long spsize, mpz_t hint)
{
  int m; /* number of already printed factors */

  ASSERT_ALWAYS(mpz_divisible_p (n, hint));
  mpz_divexact (n, n, hint);

  td_time -= seconds ();
  /* remove small primes */
  m = trial_divide (n, sp, spsize);
  td_time += seconds ();

  /* use hint */
  m = print_smooth_aux (factors, hint, methods, fm, cfm, lpb, BB, BBB, m);

  print_smooth_aux (factors, n, methods, fm, cfm, lpb, BB, BBB, m);
}

/* strip integers in l[0..n-1] which do not divide P */
unsigned long
strip (unsigned long *l, unsigned long n, mpz_t P)
{
  unsigned long i, j;

  for (i = j = 0; i < n; i++)
    if (mpz_divisible_ui_p (P, l[i]))
      l[j++] = l[i];
  return j;
}

/* Given a list L of bi-smooth cofactors, print the corresponding relations
   on stdout.
   n is the number of bi-smooth cofactors in L.
*/
void
factor (cofac_list L, unsigned long n, cado_poly pol, int lpb0,
        int lpb1, int verbose)
{
  unsigned long i, pmax, *sp0, *sp1, spsize[2];
  int nb_methods;
  facul_method_t *methods;
  mpz_t Q[2];
  struct modset_t fm, cfm;
  double BB, BBB;
  mpz_list SP;
  prime_info pi;
  double start = seconds ();
  mpz_t *norm0, *norm1;
  uint32_t *perm = L->perm;

  /* we trial divide by all primes < L^(1/2), so that any factor < L
     is necessarily prime */
  pmax = (lpb1 > lpb0) ? lpb1 : lpb0;
  pmax = (unsigned long) ceil (pow (2.0, (double) pmax / 2.0));
  pmax = 16 * pmax;

  mpz_init (Q[0]);
  mpz_init (Q[1]);

  mpz_list_init (SP);
  prime_info_init (pi);
  prime_list_poly (SP, pi, 2, pmax, pol->pols[0]);
  spsize[0] = SP->size;
  sp0 = malloc (spsize[0] * sizeof (unsigned long));
  for (i = 0; i < spsize[0]; i++)
    sp0[i] = mpz_get_ui (SP->l[i]);
  prime_info_clear (pi);
  mpz_list_clear (SP);

  mpz_list_init (SP);
  prime_info_init (pi);
  prime_list_poly (SP, pi, 2, pmax + pmax / 2, pol->pols[1]);
  spsize[1] = SP->size;
  sp1 = malloc (spsize[1] * sizeof (unsigned long));
  for (i = 0; i < spsize[1]; i++)
    sp1[i] = mpz_get_ui (SP->l[i]);
  prime_info_clear (pi);
  mpz_list_clear (SP);

  nb_methods = 30;
  if (nb_methods >= NB_MAX_METHODS)
    nb_methods = NB_MAX_METHODS - 1;
  methods = facul_make_default_strategy (nb_methods - 3, 0);

  BB = (double) pmax * (double) pmax;
  BBB = BB * (double) pmax;

  /* compute all norms F(a,b) and G(a,b) */
  norm0 = malloc (n * sizeof (mpz_t));
  norm1 = malloc (n * sizeof (mpz_t));
  for (i = 0; i < n; i++)
    {
      mpz_init (norm0[i]);
      mpz_init (norm1[i]);
      mpz_poly_homogeneous_eval_siui (norm0[i], pol->pols[0],
                                      L->a[perm[i]], L->b[perm[i]]);
      mpz_poly_homogeneous_eval_siui (norm1[i], pol->pols[1],
                                      L->a[perm[i]], L->b[perm[i]]);
    }

  for (i = 0; i < n; i++)
    {
      printf ("%" PRId64 ",%" PRIu64 ":", L->a[perm[i]], L->b[perm[i]]);

      print_smooth (Q, norm0[i], methods, &fm, &cfm, lpb0, BB, BBB,
                    sp0, spsize[0], L->R[perm[i]]);
      printf (":");

      print_smooth (Q, norm1[i], methods, &fm, &cfm, lpb1, BB, BBB,
                    sp1, spsize[1], L->A[perm[i]]);
      printf ("\n");
      fflush (stdout);
      mpz_clear (norm0[i]);
      mpz_clear (norm1[i]);
    }
  free (norm0);
  free (norm1);
  fprintf (stderr, "td_time=%.1fs sp0:%lu(%lu) sp1:%lu(%lu)\n", td_time, spsize[0], sp0[spsize[0]-1], spsize[1], sp1[spsize[1]-1]);

  mpz_clear (Q[0]);
  mpz_clear (Q[1]);

  free (sp0);
  free (sp1);
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

