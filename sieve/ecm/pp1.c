#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "prac_bc.h"
#include "pm1.h"
#include "pp1.h"


static inline void
pp1_add (residue_t r, const residue_t a, const residue_t b, 
	 const residue_t d, const modulus_t m)
{
  ASSERT (r != d);
  mod_mul (r, a, b, m);
  mod_sub (r, r, d, m);
}

static inline void
pp1_double (residue_t r, const residue_t a, const residue_t two, 
	    const modulus_t m)
{
  mod_mul (r, a, a, m);
  mod_sub (r, r, two, m);
}

void
pp1_stage1 (residue_t X, const char *code, const unsigned int l, 
	    const residue_t two, const modulus_t m)
{
  unsigned long i;
  residue_t A, B, C, t, t2;
  
  mod_init (A, m);
  mod_init (B, m);
  mod_init (C, m);
  mod_init (t, m);
  mod_init (t2, m);

  mod_set (A, X, m);

  for (i = 0; i < l; i++)
    {
      switch (code[i])
        {
	  case 0: /* Swap A, B */
            mod_swap (A, B, m);
            break;
          case 1:
            pp1_add (t, A, B, C, m);
            pp1_add (t2, t, A, B, m);
            pp1_add (B, B, t, A, m);
            mod_set (A, t2, m);
            break;
          case 2:
            pp1_add (B, A, B, C, m);
            pp1_double (A, A, two, m);
            break;
          case 3:
            pp1_add (t, B, A, C, m);
	    mod_set (C, B, m);
	    mod_set (B, t, m);
            break;
          case 4:
            pp1_add (B, B, A, C, m);
            pp1_double (A, A, two, m);
            break;
          case 5:
            pp1_add (C, C, A, B, m);
            pp1_double (A, A, two, m);
            break;
          case 6:
            pp1_double (t, A, two, m);
            pp1_add (t2, A, B, C, m);
            pp1_add (t2, t, t2, C, m);
	    mod_set (C, t2, m);
            pp1_add (t2, t, A, A, m);
	    mod_set (A, t2, m);
            mod_swap (B, C, m);
            break;
          case 7:
            pp1_add (t, A, B, C, m);
            pp1_add (t2, t, A, B, m);
	    mod_set (B, t2, m);
            pp1_double (t, A, two, m);
            pp1_add (A, A, t, A, m);
            break;
          case 8:
            pp1_add (t, A, B, C, m);
            pp1_add (C, C, A, B, m);
            mod_swap (B, t, m);
            pp1_double (t, A, two, m);
            pp1_add (t2, A, t, A, m);
	    mod_set (A, t2, m);
            break;
          case 9:
            pp1_add (C, C, B, A, m);
            pp1_double (B, B, two, m);
            break;
	  case 10: /* Init of subchain, B=A, C=A, A=2*A */
	    mod_set (B, A, m);
            mod_set (C, A, m);
            pp1_double (A, A, two, m);
            break;
          case 11:
            pp1_add (A, A, B, C, m); /* Final add */
            break;
          case 12:
            pp1_double (A, A, two, m); /* For p=2 */
            break;
	  case 13:
	    /* Rule 11, then rule 10 */
	    pp1_add (B, A, B, C, m);
            mod_set (C, B, m);
            pp1_double (A, B, two, m);
            break;
	  case 14:
	    /* Rule 3, then rule 0 */
	    mod_set (t, A, m);
	    pp1_add (A, B, A, C, m);
	    mod_set (C, B, m);
	    mod_set (B, t, m);
	    break;
	  case 15:
	    /* Rule 3, then rule 11, then rule 10 */
            pp1_add (t, B, A, C, m);
            pp1_add (C, A, t, B, m);
	    mod_set (B, C, m);
            pp1_double (A, C, two, m);
            break;
          case 16:
	    /* Rule 3, rule 0, rule 3 and rule 0, merged a bit */
	    mod_set (t, B, m);
	    pp1_add (B, B, A, C, m);
	    mod_set (C, A, m);
	    pp1_add (A, A, B, t, m);
            break;
#if 0
	  case 17:
	    /* Rule 3, then rule 11, then rule 10, then rule 3, rule 0, 
	       rule 3 and rule 0 */
            pp1_add (t, B, A, C, m);
            pp1_add (C, A, t, B, m);
	    mod_set (B, C, m);
            pp1_double (A, C, two, m);
	    mod_set (t, B, m);
	    pp1_add (B, B, A, C, m);
	    mod_set (C, A, m);
	    pp1_add (A, A, B, t, m);
            break;
#endif
	  default:
            abort ();
        }
    }

  mod_set (X, A, m);

  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (t, m);
  mod_clear (t2, m);
}

unsigned long 
pp1 (residue_t X, const modulus_t m, const pp1_plan_t *plan)
{
  residue_t b, two, save, t;
  unsigned long f;

  mod_init_noset0 (b, m);
  mod_init_noset0 (two, m);
  mod_init_noset0 (save, m);
  mod_init_noset0 (t, m);
  mod_set_ul_reduced (two, 2UL, m);
  
  /* Compute 2/7 (mod N) */
#if 0
  /* Faster method, but uses N of type unsigned long without REDC */
  if (N % 7 == 0)
    {
      f = 7;
      goto end;
    }
  
  {
    unsigned long a, n = N / 7UL, l = N % 7UL;
    /* inv7[i] stores 1/i (mod 7) */ 
    static const unsigned char inv7[7] = {0,1,4,5,2,3,6};
    /* kl1[l] stores (k*l-1)/7 for kl==1 (mod 7) and k,l < 7 */
    static const unsigned char kl1[7] = {0,0,1,2,1,2,5};
    a = inv7[l]*n + kl1[l];
    mod_set_ul_reduced (b, a, m);
    mod_neg (b, b, m);
    mod_add (b, b, b, m);
  }
#else
  /* Slow method, but works for any modulus_t type */
  mod_set_ul_reduced (b, 7UL, m);
  mod_inv (b, b, m);
  mod_add (b, b, b, m);
#endif

  pp1_stage1 (b, plan->bc, plan->bc_len, two, m);
  mod_sub (t, b, two, m);
  mod_gcd (&f, t, m);

  if (f == 1UL && plan->stage2.B2 > plan->B1)
    {
      pm1_stage2 (t, b, &(plan->stage2), two, m);
      mod_gcd (&f, t, m);
    }
  
  mod_set (X, b, m);

  mod_clear (b, m);
  mod_clear (save, m);
  mod_clear (two, m);
  mod_clear (t, m);

  return f;
}


/* Make byte code for addition chain for stage 1, and the parameters for 
   stage 2 */

void 
pp1_make_plan (pp1_plan_t *plan, const unsigned int B1, const unsigned int B2,
	       int verbose)
{
  unsigned int p;
  const unsigned int addcost = 1, doublecost = 1;
  const unsigned int compress = 1;
  
  /* Make bytecode for stage 1 */
  plan->B1 = B1;
  bytecoder_init (compress);
  for (p = 2; p <= B1; p = (unsigned int) getprime (p))
    {
      unsigned long q;
      for (q = p; q <= B1; q *= p)
	prac_bytecode (p, addcost, doublecost);
    }
  bytecoder_flush ();
  plan->bc_len = bytecoder_size ();
  plan->bc = (char *) malloc (plan->bc_len);
  ASSERT (plan->bc);
  bytecoder_read (plan->bc);
  bytecoder_clear ();

  if (verbose)
    {
      printf ("Byte code for stage 1 (length %d): ", plan->bc_len);
      for (p = 0; p < plan->bc_len; p++)
	printf ("%s%d", (p == 0) ? "" : ", ", (int) (plan->bc[p]));
      printf ("\n");
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
