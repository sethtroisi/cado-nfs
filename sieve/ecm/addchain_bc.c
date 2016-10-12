#include "cado.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include <gmp.h>

#include "portability.h"
#include "macros.h"
#include "getprime.h"
#include "gmp_aux.h" /* for nbits */
#include "prac_bc.h"
#include "addchain_bc.h"

char
addchain_find_best_r (mpz_srcptr k, const unsigned char m)
{
  int l = nbits(2*m);
  unsigned int r = k->_mp_d[0] & ((1 << l) - 1);
  if (r < m)
    return (char) r;
  else if (m - r < m)
    return (char) ((int) r) - ((int) m);
  else
  {
    l--;
    r = k->_mp_d[0] & ((1 << l) - 1);
    if (r < m)
      return (char) r;
    else /* we are sure that (m - r < m) */
    return (char) ((int) r) - ((int) m);
  }
}


static double
addchain_rec (mpz_t k, const unsigned char m, addchain_cost_srcptr opcost,
              bc_state_t *state)
{
  /* Case 0: k is <= m and (2 or odd) [ i.e., one of the precomputed points ] */
  if (mpz_cmp_ui (k, m) <= 0 && (mpz_cmp_ui (k, 2) == 0 || mpz_odd_p (k)))
  {
    /* starting point is kP bc[1] = code(k); bc[2] = 0x00 */
    if (state)
    {
      bytecoder ((literal_t) (mpz_get_ui (k) >> 1), state);
      bytecoder ((literal_t) 0x00, state);
    }
    return 0.;
  }
  /* Case 1: k == m + 2 */
  else if (mpz_cmp_ui (k, m+2) == 0)
  {
    /* starting point is kP + 2P: bc[1] = code(2); bc[2] = code(m) */
    if (state)
    {
      bytecoder ((literal_t) 0xff, state); /* special code for 2 */
      bytecoder ((literal_t) (m >> 1), state);
    }
    return opcost->add;
  }
  /* Case 2:  m+4 <= k <= 3*m-2 and k % 6 == 1 */
  else if (0 <= mpz_cmp_ui (k, m+4) && mpz_cmp_ui (k, 3*m-2)
                                    && mpz_congruent_ui_p (k, 1, 6))
  {
    unsigned int kui = mpz_get_ui (k);
    unsigned char s = (unsigned char) ((kui+2)/3);
    unsigned char r = (unsigned char) ((kui-4)/3); 
    /* starting point is sP bc[1] = code(s); bc[2] = 0x00 */
    if (state)
    {
      bytecoder ((literal_t) (s >> 1), state);
      bytecoder ((literal_t) 0x00, state);
      /* k = 2*s + r */
      bytecoder ((literal_t) (0x7f & (r >> 1)), state); /* dbladd with rP */
    }
    return opcost->dbladd;
  }
  /* Case 3:  m+4 <= k <= 3*m and k % 6 == 3 */
  else if (0 <= mpz_cmp_ui (k, m+4) && mpz_cmp_ui (k, 3*m)
                                    && mpz_congruent_ui_p (k, 3, 6))
  {
    unsigned int kui = mpz_get_ui (k);
    unsigned char s = (unsigned char) (kui/3);
    /* starting point is sP bc[1] = code(s); bc[2] = 0x00 */
    if (state)
    {
      bytecoder ((literal_t) (s >> 1), state);
      bytecoder ((literal_t) 0x00, state);
      /* k = 2*s + s = 3*s */
      bytecoder ((literal_t) (0x7f & (s >> 1)), state); /* dbladd with sP */
      /* XXX If one day, we have a tripling more efficient that dbladd, we could
       * use it here */
    }
    return opcost->dbladd;
  }
  /* Case 4:  m+4 <= k <= 3*m-4 and k % 6 == 5 */
  else if (0 <= mpz_cmp_ui (k, m+4) && mpz_cmp_ui (k, 3*m-4)
                                    && mpz_congruent_ui_p (k, 5, 6))
  {
    unsigned int kui = mpz_get_ui (k);
    unsigned char s = (unsigned char) ((kui-2)/3);
    unsigned char r = (unsigned char) ((kui+4)/3); 
    /* starting point is sP bc[1] = code(s); bc[2] = 0x00 */
    if (state)
    {
      bytecoder ((literal_t) (s >> 1), state);
      bytecoder ((literal_t) 0x00, state);
      /* k = 2*s + r */
      bytecoder ((literal_t) (0x7f & (r >> 1)), state); /* dbladd with rP */
    }
    return opcost->dbladd;
  }
  /* Case 5: 4 <= k <= 2*m-2 and k % 4 == 0 */
  else if (0 <= mpz_cmp_ui (k, m+4) && mpz_cmp_ui (k, 3*m-4)
                                    && mpz_congruent_2exp_p (k, 0, 2))
  {
    unsigned int kui = mpz_get_ui (k);
    unsigned char r1 = (unsigned char) ((kui-2)/2);
    unsigned char r2 = (unsigned char) ((kui+2)/2); 
    /* starting point is r1P + r2P: bc[1] = 0x80 | code(r1); bc[2] = code(r2) */
    if (state)
    {
      bytecoder ((literal_t) (0x80 | (r1 >> 2)), state);
      bytecoder ((literal_t) (r2 >> 1), state);
    }
    return opcost->add;
  }
  /* Even case */
  else if (mpz_even_p (k))
  {
    mpz_fdiv_q_2exp (k, k, 1); /* divides by 2 */
    double cost = addchain_rec (k, m, opcost, state);
    if (state)
      bytecoder (ADDCHAIN_DBL, state);
    return cost + opcost->dbl;
  }
  else
  {
    char r = addchain_find_best_r (k, m);
    if (r > 0)
      mpz_sub_ui (k, k, r);
    else
      mpz_add_ui (k, k, -r);
    mpz_fdiv_q_2exp (k, k, 1); /* divides by 2 */
    double cost = addchain_rec (k, m, opcost, state);
    if (state)
    {
      if (r > 0)
        bytecoder ((literal_t) (0x7f & (r >> 1)), state); /* dbladd with rP */
      else
        bytecoder ((literal_t) (0x80 | ((-r) >> 1)), state);/* dblsub with rP */
    }
    return cost + opcost->dbladd;
  }
}

double
addchain (mpz_srcptr E, const unsigned char m, addchain_cost_srcptr opcost,
          bc_state_t *state)
{
  double cost;
  mpz_t k;

  FATAL_ERROR_CHECK (m > ADDCHAIN_M_MAX,
                                "m must be smaller or equal to ADDCHAIN_M_MAX");
  FATAL_ERROR_CHECK (m % 2 == 0, "m must be odd");

  mpz_init_set (k, E);

  /* We put m as the first byte in the bytecode */
  if (state)
    bytecoder ((literal_t) m, state);

  /* Cost of the precomputation: 1 DBL (to compute 2P) and (m-1)/2 ADD (to
   * compute (2k+1)P from (2k-1)P and 2P, for k in [1..(m-1)/2])
   */
  cost = opcost->dbl_precomp + (m >> 2) * opcost->add_precomp;
  
  /* Recursively compute the cost of the addition chain */
  cost += addchain_rec (k, m, opcost, state);

  FATAL_ERROR_CHECK (mpz_cmp_ui (k, 1) != 0, "Error, k != 1 at the end");
  mpz_clear (k);
  return cost;
}

/* Bytecode an addition chain for
 *        E = prod (p^floor(log(B1)/log(p)) for f odd prime <= B1)
 */
double
addchain_bytecode (const unsigned int B1, addchain_cost_srcptr opcost,
                   bc_state_t *state)
{
  mpz_t E;

  /* E = product of all primes between 3 and B1 */
  mpz_init (E);
  mpz_set_ui (E, 1UL);
  prime_info pi;
  prime_info_init (pi);
  unsigned long p = getprime_mt (pi);
  ASSERT (p == 3);
  
  for ( ; p <= B1; p = (unsigned int) getprime_mt (pi))
  {
    for (unsigned long q = 1; q <= B1 / p; q *= p)
      mpz_mul_ui (E, E, p);
  }
  prime_info_clear (pi);
  
  /* Computes the bytecode for E: try every odd m up to ADDCHAIN_M_MAX and keep
   * the one with the lowest cost.
   */
  double chaincost, mincost = DBL_MAX;
  unsigned char best_m = 0;
  for (unsigned char m = 1 ; m <= ADDCHAIN_M_MAX ; m += 2)
  {
    double cost = addchain (E, m, opcost, NULL);
    if (cost < mincost)
	  {
	    mincost = cost;
	    best_m = m;
	  }
  }

  /* Put the best addchain into bytecode */
  chaincost = addchain (E, best_m, opcost, state);
 
  mpz_clear (E);
  return chaincost;
}
