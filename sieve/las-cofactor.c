#include "cado.h"
#include <pthread.h>

#include "utils.h"
#include "las-cofactor.h"

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
   FIXME: need to check L^k < n < B^(k+1) too.
   XXX: In doing this, pay attention to the fact that for the descent,
   we might have B^2<L.

   Possible cases, where qj represents a prime in [B,L], and rj a prime > L:
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


/*
 * This is a pre-check, based on the log-approximation in the sieve
 * array. In principle, this could save a lot of trial-divisions and
 * primality tests, if done early enough. However, this is costly and
 * since the log-approximation is not very good, we lose a lot of
 * relations.
 * Ideas to test for the future:
 *   - use it for the large holes, but not the one between 1 and 2 LP
 *   - move this check inside search_survivor(), so that the scanning is
 *   shared.
 */

void pre_check_leftover_norm(unsigned char * const SS[2], sieve_info_srcptr si)
{
    unsigned char L[2], LL[2], LLL[2], BB[2], BBB[2], BBBB[2];
    // TODO: precompute those and store them in si.
    for (int side = 0; side < 2; ++side) {
        L[side] = (unsigned char)(si->conf->sides[side]->lpb * si->sides[side]->scale);
        LL[side] = (unsigned char)(2*si->conf->sides[side]->lpb * si->sides[side]->scale);
        LLL[side] = (unsigned char)(3*si->conf->sides[side]->lpb * si->sides[side]->scale);
        size_t bb = mpz_sizeinbase (si->BB[side], 2);
        BB[side] = (unsigned char)(bb*si->sides[side]->scale);
        size_t bbb = mpz_sizeinbase (si->BBB[side], 2);
        BBB[side] = (unsigned char)(bbb*si->sides[side]->scale);
        size_t bbbb = mpz_sizeinbase (si->BBBB[side], 2);
        BBBB[side] = (unsigned char)(bbbb*si->sides[side]->scale);
    }

    const unsigned char GG=2;
    for (int i = 0; i < 1 << LOG_BUCKET_REGION; i++) {
        // TODO: fast scan
        if (SS[0][i] == 255)
            continue;
        for (int side = 0; side < 2; side++) {
            if ( (SS[side][i] > L[side]+GG && SS[side][i] < BB[side]-GG) ||
                 (SS[side][i] > LL[side]+GG && SS[side][i] < BBB[side]-GG) ||
                 (SS[side][i] > LLL[side]+GG && SS[side][i] < BBBB[side]-GG)) {
                SS[0][i] = 255;
                break;
            }
        }
    }
}


int
check_leftover_norm (const mpz_t n, sieve_info_srcptr si, int side)
{
  size_t s = mpz_sizeinbase (n, 2);
  unsigned int lpb = si->conf->sides[side]->lpb;
  unsigned int mfb = si->conf->sides[side]->mfb;

  if (s > mfb)
    return 0; /* n has more than mfb bits, which is the given limit */
  /* now n < 2^mfb */
  if (s <= lpb)
    return 1; /* case (a) */
    /* Note also that in the descent case where L > B^2, if we're below L
     * it's still fine of course, but we have no guarantee that our
     * cofactor is prime... */
  /* now n >= L=2^lpb */
  if (mpz_cmp (n, si->BB[side]) < 0)
    return 0; /* case (b) */
  /* now n >= B^2 */
  if (2 * lpb < s)
    {
      if (mpz_cmp (n, si->BBB[side]) < 0)
        return 0; /* case (e) */
      if (3 * lpb < s && mpz_cmp (n, si->BBBB[side]) < 0)
        return 0; /* case (h) */
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
static int
factor_leftover_norm (mpz_t n, mpz_array_t* const factors,
                      uint32_array_t* const multis,
                      sieve_info_srcptr si, int side)
{
  unsigned int lpb = si->conf->sides[side]->lpb;
  facul_strategy_t *strategy = si->sides[side]->strategy;
  uint32_t i, nr_factors;
  unsigned long ul_factors[16];
  int facul_code;

  /* For the moment this code can't cope with too large factors */
  ASSERT(lpb <= ULONG_BITS);

  factors->length = 0;
  multis->length = 0;

  /* If n < B^2, then n is prime, since all primes < B have been removed */
  if (mpz_cmp (n, si->BB[side]) < 0)
    {
      /* if n > L, return -1 */
      if (mpz_sizeinbase (n, 2) > lpb)
        return -1;

      if (mpz_cmp_ui (n, 1) > 0) /* 1 is special */
        {
          append_mpz_to_array (factors, n);
          append_uint32_to_array (multis, 1);
        }
      return 1;
    }

  /* use the facul library */
  //gmp_printf ("facul: %Zd\n", n);
  facul_code = facul (ul_factors, n, strategy);

  if (facul_code == FACUL_NOT_SMOOTH)
    return -1;

  ASSERT (facul_code == 0 || mpz_cmp_ui (n, ul_factors[0]) != 0);

  /* we use this mask to trap prime factors above bound */
  unsigned long oversize_mask = (-1UL) << lpb;
  if (lpb == ULONG_BITS) oversize_mask = 0;

  if (facul_code > 0)
    {
      nr_factors = facul_code;
      for (i = 0; i < nr_factors; i++)
	{
	  unsigned long r;
	  mpz_t t;
	 if (ul_factors[i] & oversize_mask) /* Larger than large prime bound? */
            return -1;
	  r = mpz_tdiv_q_ui (n, n, ul_factors[i]);
	  ASSERT_ALWAYS (r == 0UL);
	  mpz_init (t);
	  mpz_set_ui (t, ul_factors[i]);
	  append_mpz_to_array (factors, t);
	  mpz_clear (t);
	  append_uint32_to_array (multis, 1); /* FIXME, deal with repeated
						 factors correctly */
	}

      if (mpz_cmp (n, si->BB[side]) < 0)
        {
          if (mpz_sizeinbase (n, 2) > lpb)
            return -1;

          if (mpz_cmp_ui (n, 1) > 0) /* 1 is special */
            {
              append_mpz_to_array (factors, n);
              append_uint32_to_array (multis, 1);
            }
          return 1;
        }

      if (check_leftover_norm (n, si, side) == 0)
        return -1;
    }
  return 0; /* unable to completely factor n */
}


/*
This is the same function that factor_leftover_norm() but we will work
with the both norms! It's used when we want to factor these norms in the
same time and not one after the other.
*/

static int
factor_both_leftover_norms_src (mpz_t* n, mpz_array_t** factors,
			       uint32_array_t** multis,
			       sieve_info_srcptr si)
{
  int is_smooth[2] = {FACUL_MAYBE, FACUL_MAYBE};
  /* To remember if a cofactor is already factored.*/

  unsigned long** ul_factors = malloc(sizeof (*ul_factors) * 2);
  for (int i = 0; i < 2; i++)
    ul_factors[i] = calloc(16, sizeof (*ul_factors[i]));


  for (int side = 0; side < 2; side++)
    {
      unsigned int lpb = si->strategies->lpb[side];
      /* For the moment this code can't cope with too large factors */
      ASSERT(lpb <= ULONG_BITS);

      factors[side]->length = 0;
      multis[side]->length = 0;

      /* If n < B^2, then n is prime, since all primes < B have been removed */
      if (mpz_cmp (n[side], si->BB[side]) < 0)
	{
	  /* if n > L, return -1 */
	  if (mpz_sizeinbase (n[side], 2) > lpb)
	    {
	      //free
	      free (ul_factors[0]);
	      free (ul_factors[1]);
	      free (ul_factors);
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
  int* facul_code = facul_both (ul_factors, n, si->strategies, is_smooth);

  if (facul_code[0] == FACUL_NOT_SMOOTH ||
      facul_code[1] == FACUL_NOT_SMOOTH)
    {
      //free ul
      free (ul_factors[0]);
      free (ul_factors[1]);
      free (ul_factors);
      free (facul_code);
      return -1;
    }

  ASSERT (facul_code[0] == 0 || mpz_cmp_ui (n[0], ul_factors[0][0]) != 0);
  ASSERT (facul_code[1] == 0 || mpz_cmp_ui (n[1], ul_factors[1][0]) != 0);
  int ret = is_smooth[0] == 1 && is_smooth[1] == 1;
  for (int side = 0; ret == 1 && side < 2; side++)
    {
      unsigned int lpb = si->strategies->lpb[side];
      uint32_t i, nr_factors;
      /* we use this mask to trap prime factors above bound */
      unsigned long oversize_mask = (-1UL) << lpb;
      if (lpb == ULONG_BITS) oversize_mask = 0;

      if (facul_code[side] > 0)
	{
	  nr_factors = facul_code[side];
	  for (i = 0; i < nr_factors; i++)
	    {
	      unsigned long r;
	      mpz_t t;
	      if (ul_factors[side][i] & oversize_mask)
		/* Larger than large prime bound? */
		{
		  ret = -1;
		  break;
		}
	      r = mpz_tdiv_q_ui (n[side], n[side], ul_factors[side][i]);
	      ASSERT_ALWAYS (r == 0UL);
	      mpz_init (t);
	      mpz_set_ui (t, ul_factors[side][i]);
	      append_mpz_to_array (factors[side], t);
	      mpz_clear (t);
	      append_uint32_to_array (multis[side], 1);
	      /* FIXME, deal with repeated
		 factors correctly */
	    }
	  if (ret == -1)
	    break;

	  if (mpz_cmp (n[side], si->BB[side]) < 0)
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
  free (ul_factors[0]);
  free (ul_factors[1]);
  free (ul_factors);
  free (facul_code);
  /* ret = 0  => unable to completely factor n */
  return ret; 
}




/*}}}*/

/* {{{ factor_both_leftover_norms */
/*
    This function factors the leftover norms on both sides. Two
    possibilities: first, it simply calls factor_leftover_norm()
    twice, first on the side more likely to fail the smoothness test
    so that the second call can be omitted. Second, it calls
    factor_leftover_norms_src() to factor our norms on both sides in
    the same time.  
    The long-term plan is to use a more elaborate factoring strategy, 
    passing both composites to facul_*() so it can try individual 
    factoring methods on each side until smoothness or
    (likely) non-smoothness is decided.
*/

double cof_calls[2][256] = {{0},{0}};
double cof_fails[2][256] = {{0},{0}};

int
factor_both_leftover_norms(mpz_t *norm, const mpz_t BLPrat, mpz_array_t **f,
			   uint32_array_t **m, sieve_info_srcptr si)
{
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    unsigned int nbits[2];
    int first, pass = 1;

    for (int z = 0; z < 2; z++)
      nbits[z] = mpz_sizeinbase (norm[z], 2);

    /* Thread-local copy of cof_calls/cof_fails */
    double cc[2], cf[2];

    pthread_mutex_lock(&mutex);
    cc[0] = cof_calls[0][nbits[0]];
    cc[1] = cof_calls[1][nbits[1]];
    cf[0] = cof_fails[0][nbits[0]];
    cf[1] = cof_fails[1][nbits[1]];
    pthread_mutex_unlock(&mutex);

    if (cc[0] > 0. && cc[1] > 0.)
      {
	if (cf[0]* cc[1] > cf[1] * cc[0])
	  first = 0;
	else
	  first = 1;
      }
    else
      {
	/* if norm[RATIONAL_SIDE] is above BLPrat, then it might not
	 * be smooth. We factor it first. Otherwise we factor it last. */
	first = mpz_cmp (norm[RATIONAL_SIDE], BLPrat) > 0
	  ? RATIONAL_SIDE : ALGEBRAIC_SIDE;
      }

    cc[0] = cc[1] = cf[0] = cf[1] = 0.;


    /*
      The variable 'have_file' is to know if we have or not a file to
      cofactor our norms, i.e. use si->strategies or
      si->sides[*]->strategy.
    */

    int have_file = (si->strategies != NULL);
    if (have_file)
	pass = factor_both_leftover_norms_src (norm, f, m, si);
    else {
        for (int z = 0 ; pass > 0 && z < 2 ; z++){
	    int side = first ^ z;
	    cc[side] ++;
	    pass = factor_leftover_norm (norm[side], f[side], m[side],
					 si, side);
	    if (pass <= 0)
	      cf[side] ++;
	}
    }
    /* Add thread-local statistics to the global arrays */
    pthread_mutex_lock(&mutex);
    cof_calls[0][nbits[0]] += cc[0];
    cof_calls[1][nbits[1]] += cc[1];
    cof_fails[0][nbits[0]] += cf[0];
    cof_fails[1][nbits[1]] += cf[1];
    pthread_mutex_unlock(&mutex);

    return pass;
}
/*}}}*/
