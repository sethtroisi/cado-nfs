#include "cado.h"
#include <pthread.h>

#include "utils.h"
#include "las-cofactor.h"
#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"

/* {{{ factor_leftover_norm */

#define NMILLER_RABIN 1 /* in the worst case, what can happen is that a
                           composite number is declared as prime, thus
                           a relation might be missed, but this will not
                           affect correctness */
#define IS_PROBAB_PRIME(X) (0 != mpz_probab_prime_p((X), NMILLER_RABIN))
#define BITSIZE(X)      (mpz_sizeinbase((X), 2))
#define NFACTORS        8 /* maximal number of large primes */

/************************ cofactorization ********************************/

/* {{{ cofactoring area */

/* Return 0 if the leftover norm n cannot yield a relation.

   Possible cases, where qj represents a prime in [B,L], and rj a prime > L
   (assuming L < B^2, which might be false for the DLP descent):
   (0) n >= 2^mfb
   (a) n < L:           1 or q1
   (b) L < n < B^2:     r1 -> cannot yield a relation
   (c) B^2 < n < B*L:   r1 or q1*q2
   (d) B*L < n < L^2:   r1 or q1*q2 or q1*r2
   (e) L^2 < n < B^3:   r1 or q1*r2 or r1*r2 -> cannot yield a relation
   (f) B^3 < n < B^2*L: r1 or q1*r2 or r1*r2 or q1*q2*q3
   (g) B^2*L < n < L^3: r1 or q1*r2 or r1*r2
   (h) L^3 < n < B^4:   r1 or q1*r2, r1*r2 or q1*q2*r3 or q1*r2*r3 or r1*r2*r3
                        -> cannot yield a relation
*/
int
check_leftover_norm (const mpz_t n, sieve_info_srcptr si, int side)
{
  size_t s = mpz_sizeinbase (n, 2);
  unsigned int lpb = si->conf->sides[side]->lpb;
  unsigned int mfb = si->conf->sides[side]->mfb;
  unsigned int klpb;
  double nd, kB, B;

  if (s > mfb)
    return 0; /* n has more than mfb bits, which is the given limit */

  if (s <= lpb)
    return 1; /* case (a) */
  /* Note that in the case where L > B^2, if we're below L it's still fine of
     course, but we have no guarantee that our cofactor is prime... */

  nd = mpz_get_d (n);
  B = (double) si->conf->sides[side]->lim;
  kB = B * B;
  for (klpb = lpb; klpb < s; klpb += lpb, kB *= B)
    {
      /* invariant: klpb = k * lpb, kB = B^(k+1) */
      if (nd < kB) /* L^k < n < B^(k+1) */
	return 0;
    }

  // TODO: maybe we should pass the modulus to the facul machinery
  // instead of reconstructing it.
  int prime=0;
  if (s <= MODREDCUL_MAXBITS) {
      modulusredcul_t m;
      ASSERT(mpz_fits_ulong_p(n));
      modredcul_initmod_ul (m, mpz_get_ui(n));
      prime = modredcul_sprp2(m);
      modredcul_clearmod (m);
  } else if (s <= MODREDC15UL_MAXBITS) {
      modulusredc15ul_t m;
      unsigned long t[2];
      modintredc15ul_t nn;
      size_t written;
      mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, n);
      ASSERT_ALWAYS(written <= 2);
      modredc15ul_intset_uls (nn, t, written);
      modredc15ul_initmod_int (m, nn);
      prime = modredc15ul_sprp2(m);
      modredc15ul_clearmod (m);
  } else if (s <= MODREDC2UL2_MAXBITS) {
      modulusredc2ul2_t m;
      unsigned long t[2];
      modintredc2ul2_t nn;
      size_t written;
      mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, n);
      ASSERT_ALWAYS(written <= 2);
      modredc2ul2_intset_uls (nn, t, written);
      modredc2ul2_initmod_int (m, nn);
      prime = modredc2ul2_sprp2(m);
      modredc2ul2_clearmod (m);
  } else {
      prime = mpz_probab_prime_p (n, 1);
  }
  if (prime)
    return 0; /* n is a pseudo-prime larger than L */
  return 1;
}

/* This is the header-comment for the old factor_leftover_norm()
 * function, that is now deleted */
/* This function was contributed by Jerome Milan (and bugs were introduced
   by Paul Zimmermann :-).
   Input: n - the number to be factored (leftover norm). Must be composite!
              Assumed to have no factor < B (factor base bound).
          L - large prime bound is L=2^l
   Assumes n > 0.
   Return value:
          -1 if n has a prime factor larger than L
          1 if all prime factors of n are < L
          0 if n could not be completely factored
   Output:
          the prime factors of n are factors->data[0..factors->length-1],
          with corresponding multiplicities multis[0..factors->length-1].
*/

/* This is the same function as factor_leftover_norm() but it works
   with both norms! It is used when we want to factor these norms
   simultaneously and not one after the other.
*/

int
factor_both_leftover_norms(mpz_t* n, mpz_array_t** factors,
			       uint32_array_t** multis,
			       sieve_info_srcptr si)
{
  int is_smooth[2] = {FACUL_MAYBE, FACUL_MAYBE};
  /* To remember if a cofactor is already factored.*/

  mpz_t ** mpz_factors = (mpz_t **) malloc(sizeof (mpz_t*) * 2);
  for (int i = 0; i < 2; i++) {
    mpz_factors[i] = (mpz_t *) calloc(16, sizeof (*mpz_factors[i]));
    for (int j = 0; j < 16; ++j) {
        mpz_init(mpz_factors[i][j]);
    }
  }
#define FREE_MPZ_FACTOR do {          \
    for (int i = 0; i < 2; ++i) {       \
        for (int j = 0; j < 16; ++j)    \
            mpz_clear(mpz_factors[i][j]);\
        free(mpz_factors[i]);           \
    }                                   \
    free(mpz_factors);                  \
} while (0)

  for (int side = 0; side < 2; side++)
    {
      unsigned int lpb = si->strategies->lpb[side];
      double B = (double) si->conf->sides[side]->lim;
      factors[side]->length = 0;
      multis[side]->length = 0;

      /* If n < B^2, then n is prime, since all primes < B have been removed */
      if (mpz_get_d (n[side]) < B * B)
	{
	  /* if n > L, return -1 */
	  if (mpz_sizeinbase (n[side], 2) > lpb)
	    {
              FREE_MPZ_FACTOR;
	      return -1;
	    }
	  if (mpz_cmp_ui (n[side], 1) > 0) /* 1 is special */
	    {
	      append_mpz_to_array (factors[side], n[side]);
	      append_uint32_to_array (multis[side], 1);
	    }
	  is_smooth[side] = FACUL_SMOOTH;
	}
    }

  /* use the facul library */
  //gmp_printf ("facul: %Zd, %Zd\n", n[0], n[1]);
  int* facul_code = facul_both (mpz_factors, n, si->strategies, is_smooth);

  if (facul_code[0] == FACUL_NOT_SMOOTH ||
      facul_code[1] == FACUL_NOT_SMOOTH)
    {
      //free ul
      FREE_MPZ_FACTOR;
      free (facul_code);
      return -1;
    }

  ASSERT (facul_code[0] == 0 || mpz_cmp (n[0], mpz_factors[0][0]) != 0);
  ASSERT (facul_code[1] == 0 || mpz_cmp (n[1], mpz_factors[1][0]) != 0);
  int ret = is_smooth[0] == 1 && is_smooth[1] == 1;
  for (int side = 0; ret == 1 && side < 2; side++)
    {
      unsigned int lpb = si->strategies->lpb[side];
      uint32_t i, nr_factors;
      if (facul_code[side] > 0)
	{
	  nr_factors = facul_code[side];
	  for (i = 0; i < nr_factors; i++)
	    {
              if (mpz_sizeinbase(mpz_factors[side][i], 2) > lpb)
		/* Larger than large prime bound? */
		{
		  ret = -1;
		  break;
		}
	      mpz_divexact (n[side], n[side], mpz_factors[side][i]);
	      append_mpz_to_array (factors[side], mpz_factors[side][i]);
	      append_uint32_to_array (multis[side], 1);
	      /* FIXME, deal with repeated
		 factors correctly */
	    }
	  if (ret == -1)
	    break;

	  double B = (double) si->conf->sides[side]->lim;
	  if (mpz_get_d (n[side]) < B * B)
	    {
	      if (mpz_sizeinbase (n[side], 2) > lpb)
		{
		  ret = -1;
		  break;
		}

	      else if (mpz_cmp_ui (n[side], 1) > 0) /* 1 is special */
		{
		  append_mpz_to_array (factors[side], n[side]);
		  append_uint32_to_array (multis[side], 1);
		  ret = 1;
		}
	    }

	  if (check_leftover_norm (n[side], si, side) == 0)
	    {
	      ret = -1;
	      break;
	    }
	}
    }
  //free
  FREE_MPZ_FACTOR;
  free (facul_code);
  /* ret = 0  => unable to completely factor n */
  return ret; 
}


/*}}}*/


