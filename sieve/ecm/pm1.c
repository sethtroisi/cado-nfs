#include <stdio.h>
#include "pp1.h"
#include "pm1.h"
#if defined(MODREDCUL)
#include "modredc_ul.h"
#include "modredc_ul_default.h"
#define pm1 pm1_ul
#define pm1_stage1 pm1_stage1_ul
#define pp1_stage2 pp1_stage2_ul
#elif defined(MODREDC15UL)
#include "modredc_15ul.h"
#include "modredc_15ul_default.h"
#define pm1 pm1_15ul
#define pm1_stage1 pm1_stage1_15ul
#define pp1_stage2 pp1_stage2_15ul
#else
#error Please define MODREDCUL or MODREDC15UL
#endif

/* #define PARI */

void
pm1_stage1 (residue_t x, const unsigned long *E, const int E_nrwords, 
            const modulus_t m)
{
  mod_pow_mp (x, x, E, E_nrwords, m);
}


/* Looks for a factor of the modulus m, using the P-1 algorithm.
   The parameteres of P-1 are given in plan.
   Returns 1 if backtracking was used, 0 otherwise. FIXME: Right now always 
   returns 0. */

int 
pm1 (modint_t f, const modulus_t m, const pm1_plan_t *plan)
{
  residue_t x, t, X, one, two;
  unsigned int i;
  
  mod_init_noset0 (x, m);
  mod_init_noset0 (one, m);
  mod_init_noset0 (two, m);
  mod_init_noset0 (t, m);
  mod_set1 (one, m);
  mod_add (two, one, one, m);
  
  /* Stage 1, a simple exponentiation ... */
  mod_2pow_mp (x, plan->E, plan->E_nrwords, m);
  /* ... except for the backtracking part for the 2's in the exponent */
  mod_set (t, x, m);
  for (i = 0; i < plan->exp2; i++)
    {
      mod_mul (x, x, x, m);
      if (mod_is1 (x, m))
        {
          mod_set (x, t, m);
          break;
        }
    }
  
#ifdef PARI
  printf ("E = B1_exponent (%u); x = Mod(2, %lu)^E; x == %lu /* PARI */\n", 
          plan->B1, mod_getmod_ul (m), mod_get_ul (x, m));
#endif
  
  mod_sub (t, x, one, m);
  mod_gcd (f, t, m);
  
  if (mod_intcmp_ul (f, 1UL) > 0 || plan->B1 >= plan->stage2.B2)
    {
      mod_clear (one, m);
      mod_clear (two, m);
      mod_clear (t, m);
      mod_clear (x, m);
      return 0;
    }
  
  /* Compute X = x + 1/x */
  mod_init_noset0 (X, m);
  mod_inv (X, x, m);
  mod_add (X, X, x, m);
  
#ifdef PARI
  printf ("X = x+1/x; X == %lu /* PARI */\n", mod_get_ul (X, m));
#endif
  
  pp1_stage2 (t, X, &(plan->stage2), two, m);
  mod_gcd (f, t, m);
  
  mod_clear (one, m);
  mod_clear (two, m);
  mod_clear (t, m);
  mod_clear (X, m);
  mod_clear (x, m);
  return 0;
}
