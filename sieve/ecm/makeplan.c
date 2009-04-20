#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "utils.h" /* for getprime() */
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "prac_bc.h"

#define PP1_DICT_NRENTRIES 6
static size_t pp1_dict_len[PP1_DICT_NRENTRIES] = {1, 1, 2, 2, 3, 4};
static literal_t *pp1_dict_entry[PP1_DICT_NRENTRIES] = 
  {"\xB", "\xA", "\xB\xA", "\x3\x0", "\x3\xB\xA", "\x3\x0\x3\x0"};
static code_t pp1_dict_code[PP1_DICT_NRENTRIES] = {0, 0, 10, 11, 13, 14};

static bc_dict_t pp1_dict = 
  {PP1_DICT_NRENTRIES, pp1_dict_len, pp1_dict_entry, pp1_dict_code};


#define ECM_DICT_NRENTRIES 4
static size_t ecm_dict_len[ECM_DICT_NRENTRIES] = {1, 1, 2, 2};
static literal_t *ecm_dict_entry[ECM_DICT_NRENTRIES] = 
{"\xB", "\xA", "\xB\xA", "\x3\x0"};
static code_t ecm_dict_code[ECM_DICT_NRENTRIES] = {0, 0, 10, 11};

static bc_dict_t ecm_dict = 
  {ECM_DICT_NRENTRIES, ecm_dict_len, ecm_dict_entry, ecm_dict_code};


void 
pm1_make_plan (pm1_plan_t *plan, const unsigned int B1, const unsigned int B2,
	       int verbose)
{
  mpz_t E;
  unsigned int p;
  size_t tmp_E_nrwords;
  
  /* Generate the exponent for stage 1 */
  plan->exp2 = 0;
  for (p = 1; p <= B1 / 2; p *= 2)
    plan->exp2++;

  plan->B1 = B1;
  mpz_init (E);
  mpz_set_ui (E, 1UL);
  p = (unsigned int) getprime (p);
  ASSERT (p == 3);
  for ( ; p <= B1; p = (unsigned int) getprime (p))
    {
      unsigned long q;
      /* Uses p^k s.t. (p-1)p^(k-1) <= B1, except for p=2 because our 
         base 2 is a QR for primes == 1 (mod 8) already */
      for (q = 1; q <= B1 / (p - 1); q *= p)
        mpz_mul_ui (E, E, p);
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


/* Make byte code for addition chain for stage 1, and the parameters for 
   stage 2 */

void 
pp1_make_plan (pp1_plan_t *plan, const unsigned int B1, const unsigned int B2,
	       int verbose)
{
  unsigned int p;
  const unsigned int addcost = 1, doublecost = 1, bytecost = 1, changecost = 1;
  /* FIXME: without compression, the implicit leading subchain init and 
     final subchain end are superfluous! */
  const unsigned int compress = 1;
  bc_state_t *bc_state;
  
  /* Make bytecode for stage 1 */
  plan->exp2 = 0;
  for (p = 1; p <= B1 / 2; p *= 2)
    plan->exp2++;

  plan->B1 = B1;
  bc_state = bytecoder_init (compress ? &pp1_dict : NULL);
  p = (unsigned int) getprime (p);
  ASSERT (p == 3);
  for ( ; p <= B1; p = (unsigned int) getprime (p))
    {
      unsigned long q;
      for (q = 1; q <= B1 / p; q *= p)
	prac_bytecode (p, addcost, doublecost, bytecost, changecost, bc_state);
    }
  bytecoder ((literal_t) 12, bc_state);
  bytecoder_flush (bc_state);
  plan->bc_len = bytecoder_size (bc_state);
  plan->bc = (char *) malloc (plan->bc_len);
  ASSERT (plan->bc);
  bytecoder_read (plan->bc, bc_state);
  bytecoder_clear (bc_state);

  if (verbose)
    {
      int changes = 0;
      printf ("Byte code for stage 1: ");
      for (p = 0; p < plan->bc_len; p++)
        {
	  printf ("%s%d", (p == 0) ? "" : ", ", (int) (plan->bc[p]));
	  changes += (p > 0 && plan->bc[p-1] != plan->bc[p]);
        }
      printf ("\n");
      printf ("Length %d, %d code changes\n", plan->bc_len, changes);
    }
    
  /* Make stage 2 plan */
  stage2_make_plan (&(plan->stage2), B1, B2, verbose);
}

void 
pp1_clear_plan (pp1_plan_t *plan)
{
  stage2_clear_plan (&(plan->stage2));
  free (plan->bc);
  plan->bc = NULL;
  plan->bc_len = 0;
  plan->B1 = 0;
}


/* Make byte code for addition chain for stage 1, and the parameters for 
   stage 2. Parameterization chooses Brent-Suyama curves with order divisible
   by 12 (BRENT12), Montgomery with torsion 12 over Q (MONTY12) or Montgomery
   with torsion 16 over Q (MONTY16), sigma is the associated parameter.
   Extra primes controls whether some primes should be added (or left out!)
   on top of the primes and prime powers <= B1, for example to take into 
   account the known factors in the group order. */

void 
ecm_make_plan (ecm_plan_t *plan, const unsigned int B1, const unsigned int B2,
	       const int parameterization, const unsigned long sigma, 
	       const int extra_primes, const int verbose)
{
  unsigned int p, q;
  const unsigned int addcost = 6, doublecost = 5, /* TODO: find good ratio */
    bytecost = 1, changecost = 1;
  /* FIXME: without compression, the implicit leading subchain init and 
     final subchain end are superfluous! */
  const unsigned int compress = 1;
  bc_state_t *bc_state;
  
  /* If group order is divisible by 12 or 16, add two or four 2s to stage 1 */
  if (extra_primes)
    plan->exp2 = (parameterization == MONTY16) ? 4 : 2;
  else
    plan->exp2 = 0;
  for (q = 1; q <= B1 / 2; q *= 2)
    plan->exp2++;
  
  /* Make bytecode for stage 1 */
  plan->B1 = B1;
  plan->parameterization = parameterization;
  plan->sigma = sigma;
  bc_state = bytecoder_init (compress ? &ecm_dict : NULL);
  p = (unsigned int) getprime (2UL);
  ASSERT (p == 3);
  /* If group order is divisible by 12, add another 3 to stage 1 primes */
  if (extra_primes && 
      (parameterization == BRENT12 || parameterization == MONTY12))
    prac_bytecode (3, addcost, doublecost, bytecost, changecost, bc_state);
  for ( ; p <= B1; p = (unsigned int) getprime (p))
    {
      for (q = 1; q <= B1 / p; q *= p)
	prac_bytecode (p, addcost, doublecost, bytecost, changecost, bc_state);
    }
  bytecoder ((literal_t) 12, bc_state);
  bytecoder_flush (bc_state);
  plan->bc_len = bytecoder_size (bc_state);
  plan->bc = (char *) malloc (plan->bc_len);
  ASSERT (plan->bc);
  bytecoder_read (plan->bc, bc_state);
  bytecoder_clear (bc_state);
  getprime (0);

  if (verbose)
    {
      int changes = 0;
      printf ("Byte code for stage 1: ");
      for (p = 0; p < plan->bc_len; p++)
        {
	  printf ("%s%d", (p == 0) ? "" : ", ", (int) (plan->bc[p]));
	  changes += (p > 0 && plan->bc[p-1] != plan->bc[p]);
        }
      printf ("\n");
      printf ("Length %d, %d code changes\n", plan->bc_len, changes);
    }
    
  /* Make stage 2 plan */
  stage2_make_plan (&(plan->stage2), B1, B2, verbose);
}

void 
ecm_clear_plan (ecm_plan_t *plan)
{
  stage2_clear_plan (&(plan->stage2));
  free (plan->bc);
  plan->bc = NULL;
  plan->bc_len = 0;
  plan->B1 = 0;
}
