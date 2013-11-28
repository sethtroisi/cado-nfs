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
#include <gmp.h>
#include <stdio.h>
#include "las-duplicate.h"
#include "las-qlattice.h"
#include "las-coordinates.h"
#include "gmp_aux.h"


static int intlog2(uint32_t n)
{
  int l = 0;
  while (n > 1) {n >>= 1; l++;}
  return l;
}

static void
fill_in_sieve_info(sieve_info_ptr si, const uint32_t I, const uint32_t J, 
                   const unsigned long p, const int64_t a, const uint64_t b)
{
  si->I = I;
  si->J = J;
  si->conf->logI = intlog2(I);

  /* First compute the root a/b (mod p) */
  mpz_init(si->doing->p);
  mpz_init(si->doing->r);
  mpz_set_uint64(si->doing->r, b);
  mpz_set_uint64(si->doing->p, p);
  mpz_invert(si->doing->r, si->doing->r, si->doing->p);
  mpz_mul_int64(si->doing->r, si->doing->r, a);
  mpz_mod(si->doing->r, si->doing->r, si->doing->p);

}

/* Return 1 if the relation is probably a duplicate of an relation found
   when sieving the sq described by si. Return 0 if it is probably not a
   duplicate */
static int
check_one_prime(const unsigned long sq, 
                const int64_t a, const uint64_t b, const double skewness,
                const int nb_threads, sieve_info_ptr si)
{
  const unsigned long p = mpz_get_ui(si->doing->p);

  if (p == sq) /* Dummy to get rid of "unused" warning */
    return 0;

  /* Compute i,j-coordinates of this relation in the special-q lattice when
     p was used as the special-q value. */

  SkewGauss(si, skewness);  

  const uint32_t oldI = si->I, oldJ = si->J;
  /* If resulting optimal J is so small that it's not worth sieving,
     this special-q gets skipped, so relation is not a duplicate */
  if (sieve_info_adjust_IJ(si, skewness, nb_threads) == 0) {
    // fprintf(stderr, "sq = %lu, p = %lu discarded\n", sq, p);
    return 0;
  }
  
  const uint32_t I = si->I, J = si->J;
  int i;
  unsigned int j;

  if (oldI != I || oldJ != J) {
    // fprintf (stderr, "oldI = %u, I = %u, oldJ = %u, J = %u\n", oldI, I, oldJ, J);
  }
  
  int ok = ABToIJ(&i, &j, a, b, si);

  if (!ok)
    abort();

  /* If the coordinate is outside the i,j-region used when sieving
     the special-q described in si, then it's not a duplicate */
  if ((i < 0 && (uint32_t)(-i) > I) || (i > 0 && (uint32_t)i > I-1) || (j >= J))
    return 0;

  return 1;
}

/* Return 1 if the relation is probably a duplicate of a relation found
   "earlier", and 0 if it is probably not a duplicate */
int
relation_is_duplicate(relation_t *relation, const double skewness, 
                      const int nb_threads, sieve_info_srcptr si)
{
  /* If the special-q does not fit in an unsigned long, we assume it's not a
     duplicate and just move on */
  if (!mpz_fits_ulong_p(si->doing->p)) {
    return 0;
  }

  const unsigned long sq = mpz_get_ui(si->doing->p);
  const int sq_side = si->doing->side;
  unsigned long large_primes[10];
  unsigned int nr_lp = 0;
  
  const unsigned int I = si->I;
  const unsigned int J = si->J;
  /* fbb is the factor base bound on the special-q side */
  const unsigned long fbb = si->conf->sides[sq_side]->lim;
  const int64_t a = relation->a;
  const uint64_t b = relation->b;
  
  if (sq_side == RATIONAL_SIDE) {
    for (int i = 0; i < relation->nb_rp; i++) {
      unsigned long p = relation->rp[i].p;
      if (p > fbb) {
        large_primes[nr_lp++] = p;
      }
    }
  } else if (sq_side == ALGEBRAIC_SIDE) {
    for (int i = 0; i < relation->nb_ap; i++) {
      unsigned long p = relation->ap[i].p;
      if (p > fbb) {
        large_primes[nr_lp++] = p;
      }
    }
  } else
    abort();

  for (unsigned int i = 0; i < nr_lp; i++) {
    const unsigned long p = large_primes[i];

    /* Test the ordering of the relations. If the current special-q is sieved
       "later" than the large-prime p, then this relation is not a duplicate
       of what might be found when sieving p */
    /* TODO: handle cases where we sieve both sides */
    if (p >= sq)
      return 0;

    /* Create a dummy sieve_info struct with just enough info to let us use
       the lattice-reduction and coordinate-conversion functions */
    sieve_info new_si;
    fill_in_sieve_info (new_si, I, J, p, a, b);
    if (check_one_prime(sq, a, b, skewness, nb_threads, new_si))
      return 1;

    mpz_clear(new_si->doing->r);
    mpz_clear(new_si->doing->p);
  }

  return 0;
}
