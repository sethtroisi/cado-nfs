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


#include "cado.h"
#include <stdio.h>
#include <stdarg.h>
#include <gmp.h>
#include <string.h>
#include "gmp_aux.h"
#include "mpz_array.h"
#include "las-duplicate.h"
#include "las-qlattice.h"
#include "las-coordinates.h"
#include "las-norms.h"
#include "las-cofactor.h"


static const int verbose = 1;

static void
compute_a_over_b_mod_p(mpz_t r, const int64_t a, const uint64_t b, const mpz_t p)
{
  mpz_set_uint64(r, b);
  mpz_invert(r, r, p);
  mpz_mul_int64(r, r, a);
  mpz_mod(r, r, p);
}

sieve_info_ptr
fill_in_sieve_info(const mpz_t q, const mpz_t rho,
                   const int sq_side, const uint32_t I, const uint32_t J,
                   const unsigned long limits[2], facul_strategy_t *strategy[2],
                   cado_poly_ptr cpoly, siever_config_srcptr conf)
{
  sieve_info_ptr new_si;

  new_si = (sieve_info_ptr) malloc(sizeof(sieve_info));
  ASSERT_ALWAYS(new_si != NULL);
  new_si->I = I;
  new_si->J = J;

  // memset(new_si, 0, sizeof(sieve_info));
  new_si->cpoly = cpoly; /* A pointer, and polynomial not get modified */

  memmove(new_si->conf, conf, sizeof(new_si->conf));

  /* Allocate memory */
  sieve_info_init_norm_data(new_si);

  mpz_init_set(new_si->doing->p, q);
  mpz_init_set(new_si->doing->r, rho);
  new_si->doing->side = sq_side;

  for (int side = 0; side < 2; side++) {
    const unsigned long lim = limits[side];
    mpz_init (new_si->BB[side]);
    mpz_init (new_si->BBB[side]);
    mpz_init (new_si->BBBB[side]);
    mpz_ui_pow_ui (new_si->BB[side], lim, 2);
    mpz_mul_ui (new_si->BBB[side], new_si->BB[side], lim);
    mpz_mul_ui (new_si->BBBB[side], new_si->BBB[side], lim);

    new_si->sides[side]->strategy = strategy[side];
  }
  SkewGauss(new_si);  
  return new_si;
}


static sieve_info_ptr
fill_in_sieve_info_from_si(const unsigned long p, const int64_t a, const uint64_t b,
                           sieve_info_srcptr old_si)
{
  const unsigned long lim[2] = {old_si->conf->sides[0]->lim, old_si->conf->sides[1]->lim};
  facul_strategy_t *strategies[2] = {old_si->sides[0]->strategy, old_si->sides[1]->strategy};
  mpz_t sq, rho;
  mpz_init_set_ui(sq, p);
  mpz_init(rho);
  compute_a_over_b_mod_p(rho, a, b, sq);
  return fill_in_sieve_info(sq, rho, old_si->doing->side, old_si->I, old_si->J,
                            lim, strategies, old_si->cpoly, old_si->conf);
  mpz_clear(sq);
  mpz_clear(rho);
}

void
clear_sieve_info(sieve_info_ptr new_si)
{
  for (int side = 0; side < 2; side++) {
    mpz_clear(new_si->BB[side]);
    mpz_clear(new_si->BBB[side]);
    mpz_clear(new_si->BBBB[side]);
  }
  sieve_info_clear_norm_data(new_si);
  mpz_clear(new_si->doing->r);
  mpz_clear(new_si->doing->p);
  free (new_si);
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
    if (sq != 0 && !saw_sq && p == sq) {
      saw_sq = 1;
      continue;
    }
    mpz_mul_ui(cof, cof, p);
  }
  ASSERT_ALWAYS(sq == 0 || saw_sq);
}

/* Here, we assume for now that the log norm estimate is exact, i.e., that
  it produces the correctly rounded l = log_b(F(a,b)).
  To make the estimate more precise, we should call the norm calculation
  function for the correct bucket region and pick the entry corresponding to
  this i,j-coordinate. */
static unsigned char
estimate_lognorm(sieve_info_srcptr si, const int i, const unsigned int j,
                  const int side)
{
  mpz_t norm;
  mpz_init (norm);
  mpz_poly_homogeneous_eval_siui(norm, si->sides[side]->fij, i, j);
  return fb_log(mpz_get_d(norm), si->sides[side]->scale * LOG_SCALE, 0.) + GUARD;
  mpz_clear (norm);
}

/* If e == 0, returns 1. Otherwise, if b > lim, returns b.
   Otherwise returns the largest b^k with k <= e and b^k <= lim. */
static unsigned long
bounded_pow(const unsigned long b, const unsigned long e, const unsigned long lim)
{
  if (e == 0) {
    return 1;
  }
  if (b > lim) {
    return b;
  }
  unsigned long r = b;
  /* overflow is the largest ulong that can be multiplied by b without
     overflowing */
  unsigned long overflow = ULONG_MAX / b;
  for (unsigned long i = 1; i < e && r <= overflow; i++) {
    unsigned long t = r * b;
    if (t > lim) {
      break;
    }
    r = t;
  }
  return r;
}

/*
   We subtract the sieve contribution of the factor base primes.
   I.e., for each p^k || F(a,b) with p <= fbb:
     if k > 1 then we reduce k if necessary s.t. p^k <= powlim, but never to k < 1
     set l := l - log_b(p^k)
*/
static unsigned char
subtract_fb_log(const unsigned char lognorm, const relation_t *relation,
                sieve_info_srcptr si, const int side)
{
  const unsigned long fbb = si->conf->sides[side]->lim;
  const unsigned int nb_p = relation_get_nb_p(relation, side);
  unsigned char new_lognorm = lognorm;

  for (unsigned int i = 0; i < nb_p; i++) {
    const unsigned long p = relation_get_p(relation, side, i);
    const int e = relation_get_e(relation, side, i);
    ASSERT_ALWAYS(e > 0);
    if (p <= fbb) {
      const unsigned long p_pow = bounded_pow(p, e, si->conf->sides[side]->powlim);
      const unsigned char p_pow_log = fb_log(p_pow, si->sides[side]->scale * LOG_SCALE, 0.);
      if (p_pow_log > new_lognorm) {
        if (0)
          fprintf(stderr, "Warning: lognorm underflow for relation a,b = %" PRId64 ", %" PRIu64 "\n",
                  relation->a, relation->b);
        new_lognorm = 0;
      } else {
        new_lognorm -= p_pow_log;
      }
    }
  }
  return new_lognorm;
}

/* Return 1 if the relation is probably a duplicate of an relation found
   when sieving the sq described by si. Return 0 if it is probably not a
   duplicate */
int
sq_finds_relation(const unsigned long sq, const int sq_side,
                  const relation_t *relation,
                  const int nb_threads, sieve_info_srcptr old_si)
{
  mpz_array_t *f[2] = { NULL, }; /* Factors of the relation's norms */
  uint32_array_t *m[2] = { NULL, }; /* corresponding multiplicities */
  mpz_t cof[2];
  int is_dupe = 1; /* Assumed dupe until proven innocent */
  sieve_info_ptr si;
  const size_t max_large_primes = 10;
  unsigned long large_primes[2][max_large_primes];
  unsigned int nr_lp[2] = {0, 0};

  /* Extract the list of large primes for each side */
  for (int side = 0; side < 2; side++) {
    unsigned int nb_p = relation_get_nb_p(relation, side);
    for (unsigned int i = 0; i < nb_p; i++) {
      const unsigned long fbb = old_si->conf->sides[side]->lim;
      const unsigned long p = relation_get_p(relation, side, i);
      const int e = relation_get_e(relation, side, i);
      ASSERT_ALWAYS(e > 0);
      if (p > fbb) {
        for (int i = 0; i < e; i++) {
          if (nr_lp[side] < max_large_primes)
            large_primes[side][nr_lp[side]++] = p;
        }
      }
    }
    mpz_init(cof[side]);
    /* Compute the cofactor, dividing by sq on the special-q side */
    compute_cofactor(cof[side], side == sq_side ? sq : 0, large_primes[side], nr_lp[side]);
  }

  // si = get_sieve_info_from_config(las, sc, pl);
  /* Create a dummy sieve_info struct with just enough info to let us use
     the lattice-reduction and coordinate-conversion functions */
  si = fill_in_sieve_info_from_si (sq, relation->a, relation->b, old_si);
  const unsigned long r = mpz_get_ui(si->doing->r);

  const uint32_t oldI = si->I, oldJ = si->J;
  /* If resulting optimal J is so small that it's not worth sieving,
     this special-q gets skipped, so relation is not a duplicate */
  if (sieve_info_adjust_IJ(si, nb_threads) == 0) {
    if (verbose) {
      verbose_output_print(0, 1, "# DUPECHECK J too small\n");
    }
    is_dupe = 0;
    goto clear_and_exit;
  }

  sieve_info_update_norm_data(si, nb_threads);

  if (verbose) {
    verbose_output_print(0, 1, "# DUPECHECK Checking if relation (a,b) = (%" PRId64 ",%" PRIu64 ") is a dupe of sieving special-q -q0 %lu -rho %lu\n", relation->a, relation->b, sq, r);
    verbose_output_print(0, 1, "# DUPECHECK Using special-q basis a0=%" PRId64 "; b0=%" PRId64 "; a1=%" PRId64 "; b1=%" PRId64 "\n", si->a0, si->b0, si->a1, si->b1);
  }

  const uint32_t I = si->I, J = si->J;
  int i;
  unsigned int j;

  if (verbose && (oldI != I || oldJ != J)) {
    verbose_output_print(0, 1, "# DUPECHECK oldI = %u, I = %u, oldJ = %u, J = %u\n", oldI, I, oldJ, J);
  }
  
  /* Compute i,j-coordinates of this relation in the special-q lattice when
     p was used as the special-q value. */
  int ok = ABToIJ(&i, &j, relation->a, relation->b, si);

  if (!ok)
    abort();

  /* If the coordinate is outside the i,j-region used when sieving
     the special-q described in si, then it's not a duplicate */
  if ((i < 0 && (uint32_t)(-i) > I/2) || (i > 0 && (uint32_t)i > I/2-1) || (j >= J)) {
    if (verbose) {
      verbose_output_print(0, 1, "# DUPECHECK (i,j) = (%d, %u) is outside sieve region\n", i, j);
    }
    is_dupe = 0;
    goto clear_and_exit;
  }

  unsigned char remaining_lognorm[2];
  for (int side = 0; side < 2; side++) {
    const unsigned char lognorm = estimate_lognorm(si, i, j, side);
    remaining_lognorm[side] = subtract_fb_log(lognorm, relation, si, side);
    if (remaining_lognorm[side] > si->sides[side]->bound) {
      verbose_output_print(0, 1, "# DUPECHECK On side %d, remaining lognorm = %hhu > bound = %hhu\n",
              side, remaining_lognorm[side], si->sides[side]->bound);
      is_dupe = 0;
    }
  }
  verbose_output_print(0, 1, "# DUPECHECK relation had i=%d, j=%u, remaining lognorms %hhu, %hhu\n",
           i, j, remaining_lognorm[0], remaining_lognorm[1]);
  if (!is_dupe) {
    goto clear_and_exit;
  }

  /* Check that the cofactor is within the mfb bound */
  if (!check_leftover_norm (cof[sq_side], si, sq_side)) {
    if (verbose) {
      verbose_output_vfprint(0, 1, gmp_vfprintf, "# DUPECHECK cofactor %Zd is outside bounds\n", cof);
    }
    is_dupe = 0;
    goto clear_and_exit;
  }

  mpz_t BLPrat;
  mpz_init (BLPrat);
  mpz_set_ui (BLPrat, si->conf->sides[RATIONAL_SIDE]->lim);
  mpz_mul_2exp (BLPrat, BLPrat, si->conf->sides[RATIONAL_SIDE]->lpb); /* fb bound * lp bound */
  for(int side = 0 ; side < 2 ; side++) {
      f[side] = alloc_mpz_array (1);
      m[side] = alloc_uint32_array (1);
  }

  int pass = factor_both_leftover_norms(cof, BLPrat, f, m, si);

  if (pass <= 0) {
    if (verbose) {
      verbose_output_vfprint(0, 1, gmp_vfprintf, "# DUPECHECK norms not both smooth, left over factors: %Zd, %Zd\n", cof[0], cof[1]);
    }
  }

  mpz_clear(BLPrat);
  for(int side = 0 ; side < 2 ; side++) {
      clear_uint32_array (m[side]);
      clear_mpz_array (f[side]);
  }

  if (pass <= 0) {
    is_dupe = 0;
    goto clear_and_exit;
  }

clear_and_exit:
  clear_sieve_info(si);
  mpz_clear(cof[0]);
  mpz_clear(cof[1]);

  return is_dupe;
}


/* Return whether a given prime that divides the relation on a given side is
   a special-q-sieved prime. This currently tests simply that the side is the
   same as the "doing" side in the sieve_info (FIXME for sieving both sides),
   and that the prime is greater than the factor base bound. */
static int
sq_is_sieved(const unsigned long sq, int side, sieve_info_srcptr si)
{
  return side == si->doing->side && sq > si->conf->sides[side]->lim;
}

/* Return ordering between the special-q "sq" on side "side" and the special-q
   specified in si->doing: -1 if (sq,side) is earlier, 0 if they are identical,
   and 1 if (sq, side) is later */
static int
sq_cmp(const unsigned long sq, const int side, sieve_info_srcptr si)
{
  /* We assume sq are sieved in lexicographical ordering of (sq, side) */
  int cmp = -mpz_cmp_ui(si->doing->p, sq); /* negate because of swapped
                                              operands */
  if (cmp == 0)
    cmp = (side < si->doing->side) ? -1 : (side == si->doing->side) ? 0 : 1;
  return cmp;
}

/* For one special-q identified by (sq, side) (the root r is given
   implicitly by the a,b-coordinate in relation), check whether the
   relation is a duplicate of sieving that special-q: check whether 
   (sq, side) was sieved before the special-q specified in si, and
   whether sieving (sq, side) can find the relation with the sieving
   parameters specified in si->conf. */ 

int
check_one_prime(const unsigned long sq, const int side,
                const relation_t *relation, const int nb_threads,
                sieve_info_srcptr si)
{
  int is_dupe = 0;
  if (sq_is_sieved(sq, side, si) && sq_cmp(sq, side, si) < 0) {
    is_dupe = sq_finds_relation(sq, side, relation, nb_threads, si);
    if (verbose) {
      verbose_output_print(0, 1, "# DUPECHECK relation is probably%s a dupe\n",
             is_dupe ? "" : " not");
    }
  }
  return is_dupe;
}

/* Return 1 if the relation is probably a duplicate of a relation found
   "earlier", and 0 if it is probably not a duplicate */
int
relation_is_duplicate(relation_t *relation, const int nb_threads,
                      sieve_info_srcptr si)
{
  /* If the special-q does not fit in an unsigned long, we assume it's not a
     duplicate and just move on */
  if (!mpz_fits_ulong_p(si->doing->p)) {
    return 0;
  }

  int is_dupe = 0;
  /* Test the prime factors on the rational side. It would be nice to test
     both sides in a loop (int side=0; side<2; side++) but the relation_t
     stores primes on the rational and algebraic sides differently... */
  for (unsigned int i = 0; i < relation->nb_rp && !is_dupe; i++) {
    is_dupe = check_one_prime(relation->rp[i].p, RATIONAL_SIDE,
                              relation, nb_threads, si);
  }

  /* Now the primes on the algebraic side */
  for (unsigned int i = 0; i < relation->nb_ap && !is_dupe; i++) {
    is_dupe = check_one_prime(relation->ap[i].p, ALGEBRAIC_SIDE,
                              relation, nb_threads, si);
  }

  return is_dupe;
}
