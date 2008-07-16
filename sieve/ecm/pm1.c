#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "utils.h"
#include "pp1.h"
#include "pm1.h"

/* #define PARI */

void
pm1_stage1 (residue_t x, const unsigned long *E, const int E_nrwords, 
            const modulus_t m)
{
  mod_pow_mp (x, x, E, E_nrwords, m);
}


/* Looks for a factor of the modulus m, using the P-1 algorithm.
   The parameteres of P-1 are given in plan.
   If a factor is found, returns 1 and the factor in f, otherwise
   returns 0 and f is undefined.
   Upon return, x contains the end-of-stage 1 residue, so it can be
   resumed, if desired. */

unsigned long
pm1 (residue_t x, const modulus_t m, const pm1_plan_t *plan)
{
  residue_t t, X, one, two;
  unsigned long f;
  
  mod_init_noset0 (one, m);
  mod_init_noset0 (two, m);
  mod_set_ul_reduced (one, 1UL, m);
  mod_add (two, one, one, m);
  
  /* Stage 1, a simple exponentiation */
  mod_2pow_mp (x, two, plan->E, plan->E_nrwords, plan->E_mask, m);
  
#ifdef PARI
  printf ("E = B1_exponent (%u); x = Mod(2, %lu)^E; x == %lu /* PARI */\n", 
          plan->B1, mod_getmod_ul (m), mod_get_ul (x, m));
#endif

  mod_init_noset0 (t, m);
  mod_sub (t, x, one, m);
  mod_gcd (&f, t, m);

  if (f > 1UL || plan->B1 >= plan->stage2.B2)
    {
      mod_clear (one, m);
      mod_clear (two, m);
      mod_clear (t, m);
      return f;
    }

  /* Compute X = x + 1/x. TODO: Speed this up, use precomputed 2^{3w} % m? */
  mod_init_noset0 (X, m);
  mod_inv (X, x, m);
  mod_add (X, X, x, m);
  
#ifdef PARI
  printf ("X = x+1/x; X == %lu /* PARI */\n", mod_get_ul (X, m));
#endif
  
  pp1_stage2 (t, X, &(plan->stage2), two, m);
  mod_gcd (&f, t, m);
  
  mod_clear (one, m);
  mod_clear (two, m);
  mod_clear (t, m);
  mod_clear (X, m);
  return f;
}



void 
pm1_make_plan (pm1_plan_t *plan, const unsigned int B1, const unsigned int B2,
	       int verbose)
{
  mpz_t E;
  unsigned int p;
  size_t tmp_E_nrwords;
  
  /* Generate the exponent for stage 1 */
  plan->B1 = B1;
  mpz_init (E);
  mpz_set_ui (E, 1UL);
  for (p = 2; p <= B1; p = (unsigned int) getprime (p))
    {
      unsigned long q;
      /* FIXME: use p^k s.t. (p-1)p^(k-1) <= B1 instead, except for
         p=2 because our base 2 is a QR for primes == 1 (mod 8) already */
      for (q = p; q * p < B1; q *= p);
      mpz_mul_ui (E, E, q);
    }
  
  if (verbose)
    gmp_printf ("pm1_make_plan: E = %Zd;\n", E);
  
  plan->E = mpz_export (NULL, &tmp_E_nrwords, -1, sizeof(unsigned long),
                        0, 0, E);
  plan->E_nrwords = (unsigned int) tmp_E_nrwords;
  mpz_clear (E);
  /* Find highest set bit in E. */
  ASSERT (plan->E[plan->E_nrwords - 1] != 0);
  plan->E_mask = ~0UL - (~0UL >> 1); /* Only MSB set */
  while ((plan->E[plan->E_nrwords - 1] & plan->E_mask) == 0UL)
    plan->E_mask >>= 1;
  getprime (0);
  
  stage2_make_plan (&(plan->stage2), B1, B2, verbose);
}


void 
pm1_clear_plan (pm1_plan_t *plan)
{
  stage2_clear_plan (&(plan->stage2));

  free (plan->E);
  plan->E = NULL;
  plan->E_nrwords = 0;
  plan->B1 = 0;
}
