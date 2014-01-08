/*

We want a function that checks whether a relation is, with high probability, 
a duplicate of an earlier relation. 

"Earlier" requires that we define some ordering on relations. Note that we may
sieve special-q on both sides for some factorisations, in particular SNFS.
It is possible that the same prime occurs on both sides if it divides the 
two polynomials' resultant, which likewise may occur for SNFS.

We define the sieving order as sieving special-q on the algebraic side of 
norm n if n is in the special-q range for the algebraic side, 
then on the rational side of norm n if n is in the special-q range for the 
rational side, then increasing n by 1.

Let Q_0 be the primes in the special-q range on the algebraic side, 
and Q_1 on the rational side. If we sieve only one side i \in {0,1}, 
then the set Q_{1-i} is empty.

When we sieve special-q on only one side i, a relation is "earlier" 
simply if it contains a prime on side i that is in Q_i and smaller than 
the current q. (1)

If we sieve both sides, a relation is earlier if
- (1) holds, or
- we sieve q on the algebraic side and it contains a prime q' on the rational 
  side in Q_1 less than q
- we sieve q on the rational side and it contains a prime q' on the algebraic 
  side in Q_0 less than or equal to q

Then we must go through the various steps the siever does to identify 
relations and check the conditions that cause the relation to be found or 
missed.

The relation is a dupe if:
- It may have been found by an "earlier" special-q $q'$
- The relation must be in the I,J sieving area that was used when sieving in 
  the $q'$-lattice
- The estimated norm, minus the log of the FB primes, is within the report 
  threshold. This is tricky, as there are many shortcuts to norm estimation,
  and the contribution of FB primes depends on the choice of log scale and 
  on limits on prime powers that are sieved
- The cofactors of the two norms, without q' but with q on the special-q side, 
  are within the cofactor bounds mfb[ar]
- The cofactors on both sides can be factored by the cofactorization routines 
  that were used when sieving q'.

Thus the function to check for duplicates needs the following information:

- The relation with fully factored norms
- The special-q sieve range on the two sides
- For each q' that was sieved "earlier"
  - The parameters used for generating the special-q' lattice (skewness)
  - The I,J values
  - The log scale
  - The cofactor bounds mfb[ar]
  - The cofactorization parameters

*/


#if 0
    From relation.h, copy-pasted for easy reference:

    typedef struct {
      unsigned long p;      /* rational prime */
      int e;                /* exponent (may want negative exponent in sqrt) */
    } rat_prime_t;

    typedef struct {
      unsigned long p;      /* algebraic prime */
      unsigned long r;      /* corresponding root: r = a/b mod p */
      int e;                /* exponent (may want negative exponent in sqrt) */
    } alg_prime_t;

    typedef struct {
      int64_t a;		/* only a is allowed to be negative */
      uint64_t b;
      rat_prime_t *rp;	/* array of rational primes */
      alg_prime_t *ap;	/* array of algebraic primes */
      uint8_t nb_rp;	/* number of rational primes */
      uint8_t nb_ap;        /* number of algebraic primes */
      uint8_t nb_rp_alloc;	/* allocated space for rp */
      uint8_t nb_ap_alloc;	/* allocated space for ap */
    } relation_t;
#endif




#include "cado.h"
#include <stdio.h>
#include <gmp.h>
#include <string.h>
#include "gmp_aux.h"
#include "mpz_array.h"
#include "las-duplicate.h"
#include "las-qlattice.h"
#include "las-coordinates.h"
#include "las-norms.h"


static const int verbose = 1;

static void
fill_in_sieve_info(sieve_info_ptr new_si, const unsigned long p,
                   const int64_t a, const uint64_t b, const double skewness,
                   sieve_info_srcptr old_si)
{
  // memset(new_si, 0, sizeof(sieve_info));
  new_si->I = old_si->I;
  new_si->J = old_si->J;

  new_si->cpoly = old_si->cpoly; /* A pointer, and polynomial not get modified */

  new_si->conf->side = old_si->conf->side;
  new_si->conf->logI = old_si->conf->logI;
  for(int side = 0; side < 2; side++) {
    new_si->conf->sides[side]->lim = old_si->conf->sides[side]->lim;
    new_si->conf->sides[side]->lpb = old_si->conf->sides[side]->lpb;
    new_si->conf->sides[side]->mfb = old_si->conf->sides[side]->mfb;
    new_si->conf->sides[side]->lambda = old_si->conf->sides[side]->lambda;
    
    new_si->sides[side]->logmax = old_si->sides[side]->logmax;
    new_si->sides[side]->strategy = old_si->sides[side]->strategy;
    int d = new_si->cpoly->pols[side]->deg;
    mpz_poly_init (new_si->sides[side]->fij, d);
    new_si->sides[side]->fijd = (double *) malloc((d + 1) * sizeof(double));
  }

  /* Compute the root a/b (mod p) */
  mpz_init(new_si->doing->p);
  mpz_init(new_si->doing->r);
  mpz_set_uint64(new_si->doing->r, b);
  mpz_set_uint64(new_si->doing->p, p);
  mpz_invert(new_si->doing->r, new_si->doing->r, new_si->doing->p);
  mpz_mul_int64(new_si->doing->r, new_si->doing->r, a);
  mpz_mod(new_si->doing->r, new_si->doing->r, new_si->doing->p);

  for (int side = 0; side < 2; side++) {
    mpz_init_set(new_si->BB[side], old_si->BB[side]);
    mpz_init_set(new_si->BBB[side], old_si->BBB[side]);
    mpz_init_set(new_si->BBBB[side], old_si->BBBB[side]);
  }
  SkewGauss(new_si, skewness);  
}

static void
clear_sieve_info(sieve_info_ptr new_si)
{
  for(int side = 0; side < 2; side++)
    mpz_poly_clear (new_si->sides[side]->fij);
  for (int side = 0; side < 2; side++) {
    mpz_clear(new_si->BB[side]);
    mpz_clear(new_si->BBB[side]);
    mpz_clear(new_si->BBBB[side]);
    free(new_si->sides[side]->fijd);
  }
  mpz_clear(new_si->doing->r);
  mpz_clear(new_si->doing->p);
}

/* Compute the product of the nr_lp large primes in large_primes,
   skipping sq. If sq occurs multiple times in large_primes, it is
   skipped only once. */
static void
compute_cofactor(mpz_t cof, const unsigned long sq, 
                 const unsigned long *large_primes, const int nr_lp)
{
  mpz_set_ui (cof, 1UL);
  int saw_sq = 0;
  for (int i = 0; i < nr_lp; i++) {
    const unsigned long p = large_primes[i];
    if (!saw_sq && p == sq) {
      saw_sq = 1;
      continue;
    }
    mpz_mul_ui(cof, cof, p);
  }
}

#define WRAP(x) x

/* Return 1 if the relation is probably a duplicate of an relation found
   when sieving the sq described by si. Return 0 if it is probably not a
   duplicate */
int
check_one_prime(const unsigned long sq MAYBE_UNUSED, 
                const int64_t a, const uint64_t b,
                const double skewness, const int nb_threads,
                mpz_t cof[2], int sq_side, sieve_info_ptr si)
{
  const unsigned long p = mpz_get_ui(si->doing->p);
  const unsigned long r = mpz_get_ui(si->doing->r);
  mpz_array_t *f[2] = { NULL, };
  uint32_array_t *m[2] = { NULL, }; /* corresponding multiplicities */

  /* Compute i,j-coordinates of this relation in the special-q lattice when
     p was used as the special-q value. */

  if (verbose) {
    printf("# DUPECHECK Checking if relation (a,b) = (%" PRId64 ",%" PRIu64 ") is a dupe of sieving special-q q=%lu; rho=%lu\n", a, b, p, r);
    printf("# DUPECHECK Using special-q basis a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 "\n", si->a0, si->b0, si->a1, si->b1);
  }

  const uint32_t oldI = si->I, oldJ = si->J;
  /* If resulting optimal J is so small that it's not worth sieving,
     this special-q gets skipped, so relation is not a duplicate */
  if (sieve_info_adjust_IJ(si, skewness, nb_threads) == 0) {
    if (verbose) {
      printf("# DUPECHECK J too small\n");
    }
    return 0;
  }

  sieve_info_update_norm_data(si, nb_threads);

  const uint32_t I = si->I, J = si->J;
  int i;
  unsigned int j;

  if (verbose && (oldI != I || oldJ != J)) {
    printf ("# DUPECHECK oldI = %u, I = %u, oldJ = %u, J = %u\n", oldI, I, oldJ, J);
  }
  
  int ok = ABToIJ(&i, &j, a, b, si);

  if (!ok)
    abort();

  /* If the coordinate is outside the i,j-region used when sieving
     the special-q described in si, then it's not a duplicate */
  if ((i < 0 && (uint32_t)(-i) > I/2) || (i > 0 && (uint32_t)i > I/2-1) || (j >= J)) {
    if (verbose) {
      printf("# DUPECHECK (i,j) = (%d, %u) is outside sieve region\n", i, j);
    }
    return 0;
  }

  /* Check that the cofactor is within the mfb bound */
  if (!check_leftover_norm (cof[sq_side], si, sq_side)) {
    if (verbose) {
      gmp_printf("# DUPECHECK cofactor %Zd is outside bounds\n", cof);
    }
    return 0;
  }

  mpz_t BLPrat;
  mpz_init (BLPrat);
  mpz_set_ui (BLPrat, si->conf->sides[RATIONAL_SIDE]->lim);
  mpz_mul_2exp (BLPrat, BLPrat, si->conf->sides[RATIONAL_SIDE]->lpb); /* fb bound * lp bound */
  for(int side = 0 ; side < 2 ; side++) {
      f[side] = alloc_mpz_array (8);
      m[side] = alloc_uint32_array (8);
  }

#if 0
  int
  factor_both_leftover_norms(mpz_t *norm, const mpz_t BLPrat, mpz_array_t **f,
                             uint32_array_t **m, sieve_info_srcptr si)
#endif                       
  int pass = factor_both_leftover_norms(cof, BLPrat, f, m, si);

  if (pass <= 0) {
    if (verbose) {
      gmp_printf("# DUPECHECK norms not both smooth, left over factors: %Zd, %Zd\n", cof[0], cof[1]);
    }
  }

  mpz_clear(BLPrat);
  for(int side = 0 ; side < 2 ; side++) {
      clear_uint32_array (m[side]);
      clear_mpz_array (f[side]);
  }

  if (pass <= 0) {
    return 0;
  }

  return 1;
}

/* Return 1 if the relation is probably a duplicate of a relation found
   "earlier", and 0 if it is probably not a duplicate */
int
relation_is_duplicate(relation_t *relation, const double skewness, 
                      const int nb_threads, sieve_info_srcptr si)
{
  int is_dupe = 0;
  /* If the special-q does not fit in an unsigned long, we assume it's not a
     duplicate and just move on */
  if (!mpz_fits_ulong_p(si->doing->p)) {
    return 0;
  }

  const unsigned long sq = mpz_get_ui(si->doing->p);
  const int sq_side = si->doing->side;
  unsigned long large_primes[2][10];
  unsigned int nr_lp[2] = {0, 0};
  mpz_t cof[2];

  ASSERT_ALWAYS(0 <= RATIONAL_SIDE && RATIONAL_SIDE < 2);
  ASSERT_ALWAYS(0 <= ALGEBRAIC_SIDE && ALGEBRAIC_SIDE < 2);
  ASSERT_ALWAYS(RATIONAL_SIDE != ALGEBRAIC_SIDE);
  ASSERT_ALWAYS(sq_side == RATIONAL_SIDE || sq_side == ALGEBRAIC_SIDE);

  /* fbb is the factor base bound on the special-q side */
  const int64_t a = relation->a;
  const uint64_t b = relation->b;

  mpz_init(cof[0]);
  mpz_init(cof[1]);

  /* Extract the list of large primes for the rational side */
  for (int i = 0; i < relation->nb_rp; i++) {
    const unsigned long fbb = si->conf->sides[RATIONAL_SIDE]->lim;
    const unsigned long p = relation->rp[i].p;
    const int e = relation->rp[i].e;
    ASSERT_ALWAYS(e > 0);
    if (p > fbb) {
      for (int i = 0; i < e; i++) {
        large_primes[RATIONAL_SIDE][nr_lp[RATIONAL_SIDE]++] = p;
      }
    }
  }
  /* Extract the list of large primes for the algebraic side */
  for (int i = 0; i < relation->nb_ap; i++) {
    const unsigned long fbb = si->conf->sides[ALGEBRAIC_SIDE]->lim;
    unsigned long p = relation->ap[i].p;
    const int e = relation->ap[i].e;
    ASSERT_ALWAYS(e > 0);
    if (p > fbb) {
      for (int i = 0; i < e; i++)
        large_primes[ALGEBRAIC_SIDE][nr_lp[ALGEBRAIC_SIDE]++] = p;
    }
  }

  for (unsigned int i = 0; i < nr_lp[sq_side] && !is_dupe; i++) {
    const unsigned long p = large_primes[sq_side][i];

    /* Test the ordering of the relations. If the current special-q is sieved
       "earlier" than the large-prime p, then this relation is not a duplicate
       of what might be found when sieving p */
    /* TODO: handle cases where we sieve both sides */
    if (p >= sq)
      continue;

    /* Compute the cofactor for the non-special-q side. The facul library
       modifies the mpz_t, so we have to re-compute it again each time */
    compute_cofactor(cof[1-sq_side], 0, large_primes[1-sq_side], nr_lp[1-sq_side]);

    /* Create a dummy sieve_info struct with just enough info to let us use
       the lattice-reduction and coordinate-conversion functions */
    sieve_info new_si;
    fill_in_sieve_info (new_si, p, a, b, skewness, si);
    compute_cofactor(cof[sq_side], p, large_primes[sq_side], nr_lp[sq_side]);
    is_dupe = check_one_prime(sq, a, b, skewness, nb_threads, cof, sq_side, new_si);
    if (verbose) {
      printf("# DUPECHECK relation is probably%s a dupe\n", is_dupe ? "" : " not");
    }
    clear_sieve_info(new_si);
  }
  mpz_clear(cof[0]);
  mpz_clear(cof[1]);

  return is_dupe;
}
