#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "utils.h"
#include "pm1.h"

/* #define PARI */

void
pm1_stage1 (residue_t x, const unsigned long *E, const int E_nrwords, 
            const modulus_t m)
{
  mod_pow_mp (x, x, E, E_nrwords, m);
}


unsigned long 
pm1_stage2 (residue_t r, const residue_t X, const stage2_plan_t *plan, 
	    const residue_t two, const modulus_t m)
{
  residue_t Xd, Xid, Xid1, a, t;
  residue_t *Xj;
  unsigned int k, l;
#ifdef PARI
  unsigned int id;
#endif
  
  mod_init_noset0 (t, m);

#define FAST_XJ_INIT
#ifndef FAST_XJ_INIT
  /* Compute Xd = V_d(X) */
  mod_V_ul (Xd, X, two, plan->d, m);

#ifdef PARI
      printf ("d = %u; Xd = V(d, X); Xd == %lu /* PARI */\n",
	      plan->d, mod_get_ul (Xd, m));
#endif
  
  /* Compute V_j(X) for j in S_1. Really slow, use common addition 
     chain below instead */
  Xj = malloc (plan->s1 * sizeof(residue_t));
  ASSERT (Xj != NULL);
  for (k = 0; k < plan->s1; k++)
    {
      mod_init_noset0 (Xj[k], m);
      mod_V_ul (Xj[k], X, two, plan->S1[k], m);
#ifdef PARI
      printf ("V(%u, X) == %lu /* = Xj[%d] */ /* PARI */\n",
	      plan->S1[k], mod_get_ul (Xj[k], m), k);
#endif
    }
#else /* if FAST_XJ_INIT */
  /* Faster way: compute all the j, 1 <= j < d/2, gcd(j,d)=1 with two 
     arithmetic progressions 1+6k and 5+6k (this assumes 6|d).
     We need two values of each progression (1, 7 and 5, 11) and the 
     common difference 6. These can be computed with the Lucas chain
     1, 2, 3, 5, 6, 7, 11 at the cost of 6 multiplies. */
  ASSERT (plan->d % 6 == 0);
  {
    residue_t ap1_0, ap1_1, ap5_0, ap5_1, X2, X6;
    int i1, i5;
    mod_init_noset0 (ap1_0, m);
    mod_init_noset0 (ap1_1, m);
    mod_init_noset0 (ap5_0, m);
    mod_init_noset0 (ap5_1, m);
    mod_init_noset0 (X6, m);
    mod_init_noset0 (X2, m);
    
    /* Init ap1_0 = V_1(X), ap1_1 = V_7(X), ap5_0 = V_5(X), ap5_1 = V_11(X)
       and X6 = V_6(X) */
    mod_set (ap1_0, X, m); /* ap1_0 = V_1(X) = X */
    mod_mul (X2, X, X, m);
    mod_sub (X2, X2, two, m); /* X2 = V_2(X) = X^2 - 2 */
    mod_mul (X6, X2, X, m);
    mod_sub (X6, X6, X, m); /* V_3(X) = V_2(X) * V_1(X) - V_1(X) */
    mod_mul (ap5_0, X6, X2, m);
    mod_sub (ap5_0, ap5_0, X, m); /* V_5(X) = V_3(X) * V_2(X) - V_1(X) */
    mod_mul (X6, X6, X6, m);
    mod_sub (X6, X6, two, m); /* V_6(X) = V_3(X)*V_3(X) - 2 */
    mod_mul (ap1_1, X6, X, m);
    mod_sub (ap1_1, ap1_1, ap5_0, m); /* V_7(X) = V_6(X) * V_1(X) - V_5(X) */
    mod_mul (ap5_1, X6, ap5_0, m);
    mod_sub (ap5_1, ap5_1, X, m); /* V_11(X) = V_6(X) * V_5(X) - V_1(X) */
    
    mod_clear (X2, m);
    
#ifdef PARI
    printf ("V(1, X) == %lu /* PARI */\n", mod_get_ul (ap1_0, m));
    printf ("V(5, X) == %lu /* PARI */\n", mod_get_ul (ap5_0, m));
    printf ("V(7, X) == %lu /* PARI */\n", mod_get_ul (ap1_1, m));
    printf ("V(11, X) == %lu /* PARI */\n", mod_get_ul (ap5_1, m));
#endif
    
    /* Now we generate all the V_j(X) for j in S_1 */
    Xj = malloc (plan->s1 * sizeof(residue_t));
    ASSERT (Xj != NULL);
    
    /* We treat the first two manually because those might correspond 
       to ap1_0 = V_1(X) and ap5_0 = V_5(X) */
    k = 0;
    if (plan->s1 > k && plan->S1[k] == 1)
      {
        mod_init_noset0 (Xj[k], m);
        mod_set (Xj[k++], ap1_0, m);
      }
    if (plan->s1 > k && plan->S1[k] == 5)
      {
        mod_init_noset0 (Xj[k], m);
        mod_set (Xj[k++], ap5_0, m);
      }
    
    i1 = 7;
    i5 = 11;
    while (k < plan->s1)
      {
        if (plan->S1[k] == i1)
          {
            mod_init_noset0 (Xj[k], m);
            mod_set (Xj[k], ap1_1, m);
#ifdef PARI
	    printf ("V(%u, X) == %lu /* = Xj[%d] */ /* PARI */\n",
                i1, mod_get_ul (Xj[k], m), k);
#endif
	    k++;
	    continue;
          }
        if (plan->S1[k] == i5)
          {
            mod_init_noset0 (Xj[k], m);
            mod_set (Xj[k], ap5_1, m);
#ifdef PARI
	    printf ("V(%u, X) == %lu /* = Xj[%d] */ /* PARI */\n",
		    i5, mod_get_ul (Xj[k], m), k);
#endif
	    k++;
	    continue;
          }
	
        mod_mul (t, ap1_1, X6, m);
        mod_sub (t, t, ap1_0, m);
        mod_set (ap1_0, ap1_1, m);
        mod_set (ap1_1, t, m);
        i1 += 6;
	
        mod_mul (t, ap5_1, X6, m);
        mod_sub (t, t, ap5_0, m);
        mod_set (ap5_0, ap5_1, m);
        mod_set (ap5_1, t, m);
        i5 += 6;
#ifdef PARI
	printf ("V(%u, X) == %lu /* new ap1_1 */ /* PARI */\n",
		i1, mod_get_ul (ap1_1, m));
	printf ("V(%u, X) == %lu /* new ap5_1 */ /* PARI */\n",
		i5, mod_get_ul (ap5_1, m));
#endif
      }
    
    /* Also compute Xd = V_d(X) while we've got V_6(X) */
    mod_V_ul (Xd, X6, two, plan->d / 6, m);

#ifdef PARI
    printf ("d = %u; Xd = V(d, X); Xd == %lu /* PARI */\n",
	    plan->d, mod_get_ul (Xd, m));
#endif
    
    mod_clear (ap1_0, m);
    mod_clear (ap1_1, m);
    mod_clear (ap5_0, m);
    mod_clear (ap5_1, m);
    mod_clear (X6, m);
    mod_clear (X2, m);
  }
#endif /* if FAST_XJ_INIT */
  
  mod_init_noset0 (Xid, m);
  mod_init_noset0 (Xid1, m);
  mod_init_noset0 (a, m);
  mod_set_ul_reduced (a, 1UL, m);
  l = 0;
  
  {
    /* Compute V_{i0 * d}(X) and V_{(i0 + 1) * d}(X) so we can
       compute the remaining V_{id}(X) via an arithmetic progression.
       TODO: init both with the same binary chain. */
    
    mod_V_ul (Xid, Xd, two, plan->i0, m);
    mod_V_ul (Xid1, Xd, two, plan->i0 + 1, m);
#ifdef PARI
    printf ("V(%u, X) == %lu /* PARI */\n",
	    plan->d * plan->i0, mod_get_ul (Xid, m));
    printf ("V(%u, X) == %lu /* PARI */\n",
	    plan->d * (plan->i0 + 1), mod_get_ul (Xid1, m));
#endif
    
    while (plan->pairs[l] != NEXT_PASS) 
      {
	while (plan->pairs[l] < NEXT_D && plan->pairs[l] < NEXT_PASS)
	  {
	    mod_sub (t, Xid, Xj[plan->pairs[l]], m);
	    mod_mul (a, a, t, m);
	    l++;
	  }
	
	/* Advance i by 1 */
	if (plan->pairs[l] == NEXT_D)
	  {
	    mod_mul (t, Xid1, Xd, m);
	    mod_sub (t, t, Xid, m);
	    mod_set (Xid, Xid1, m);
	    mod_set (Xid1, t, m);
	    l++; /* Skip over NEXT_D */
#ifdef PARI
	    id += plan->d;
	    printf ("V(%u, X) == %lu /* PARI */\n",
		    id, mod_get_ul (Xid, m));
	    printf ("V(%u, X) == %lu /* PARI */\n",
		    id + plan->d, mod_get_ul (Xid1, m));
#endif
	  }
      }
    l++; /* Skip over NEXT_PASS */
  }
  
  mod_set (r, a, m);

  for (k = 0; k < plan->s1; k++)
    mod_clear (Xj[k], m);

  free (Xj);

  for (k = 0; k < plan->s1; k++)
    mod_clear (Xj[k], m);
  mod_clear (Xid, m);
  mod_clear (Xid1, m);
  mod_clear (a, m);
  mod_clear (t, m);
  
  return 1UL;
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
  
  pm1_stage2 (t, X, &(plan->stage2), two, m);
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
