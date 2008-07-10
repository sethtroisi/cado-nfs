#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "prac_bc.h"
#include "pm1.h"
#include "pp1.h"

#define dup(a,b) mod_mul(a,b,b,m);mod_sub(a,a,two,m);
#define add(a,b,c,d) mod_mul(a,b,c,m);mod_sub(a,a,d,m);
/* Safe version of add that can be used if "a" and "d" are identical */
#define adds(a,b,c,d) mod_set(t3,d,m);mod_mul(a,b,c,m);mod_sub(a,a,t3,m);
#define set(a,b) mod_set(a,b,m);
#define swap(a,b) mod_swap(a,b,m);

static void
pp1_range_50 (residue_t r, const residue_t b, const residue_t two, 
              const modulus_t m)
{
  residue_t A, B, C, t, t2, t3;

  mod_init_noset0 (A, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (C, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (t2, m);
  mod_init_noset0 (t3, m);

  mod_set (A, b, m);
 
#include "chain_50.c"
  
  mod_set (r, A, m);
  
  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (t, m);
  mod_clear (t2, m);
  mod_clear (t3, m);
}

static void
pp1_range_50_100 (residue_t r, const residue_t b, const residue_t two, 
                  const modulus_t m)
{
  residue_t A, B, C, t, t2, t3;

  mod_init_noset0 (A, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (C, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (t2, m);
  mod_init_noset0 (t3, m);

  mod_set (A, b, m);
 
#include "chain_50_100.c"
  
  mod_set (r, A, m);
  
  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (t, m);
  mod_clear (t2, m);
  mod_clear (t3, m);
}

static void
pp1_range_100_150 (residue_t r, const residue_t b, const residue_t two, 
                   const modulus_t m)
{
  residue_t A, B, C, t, t2, t3;

  mod_init_noset0 (A, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (C, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (t2, m);
  mod_init_noset0 (t3, m);

  mod_set (A, b, m);
 
#include "chain_100_150.c"
  
  mod_set (r, A, m);
  
  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (t, m);
  mod_clear (t2, m);
  mod_clear (t3, m);
}

static void
pp1_range_150_200 (residue_t r, const residue_t b, const residue_t two, 
                   const modulus_t m)
{
  residue_t A, B, C, t, t2, t3;

  mod_init_noset0 (A, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (C, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (t2, m);
  mod_init_noset0 (t3, m);

  mod_set (A, b, m);
 
#include "chain_150_200.c"
  
  mod_set (r, A, m);
  
  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (t, m);
  mod_clear (t2, m);
  mod_clear (t3, m);
}

static void
pp1_range_200_300 (residue_t r, const residue_t b, const residue_t two, 
                   const modulus_t m)
{
  residue_t A, B, C, t, t2, t3;

  mod_init_noset0 (A, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (C, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (t2, m);
  mod_init_noset0 (t3, m);

  mod_set (A, b, m);
 
#include "chain_200_300.c"
  
  mod_set (r, A, m);
  
  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (t, m);
  mod_clear (t2, m);
  mod_clear (t3, m);
}

static void
pp1_range_300_400 (residue_t r, const residue_t b, const residue_t two, 
                   const modulus_t m)
{
  residue_t A, B, C, t, t2, t3;

  mod_init_noset0 (A, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (C, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (t2, m);
  mod_init_noset0 (t3, m);

  mod_set (A, b, m);
 
#include "chain_300_400.c"
  
  mod_set (r, A, m);
  
  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (t, m);
  mod_clear (t2, m);
  mod_clear (t3, m);
}

static void
pp1_range_400_500 (residue_t r, const residue_t b, const residue_t two, 
                   const modulus_t m)
{
  residue_t A, B, C, t, t2, t3;

  mod_init_noset0 (A, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (C, m);
  mod_init_noset0 (t, m);
  mod_init_noset0 (t2, m);
  mod_init_noset0 (t3, m);

  mod_set (A, b, m);
 
#include "chain_400_500.c"
  
  mod_set (r, A, m);
  
  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (t, m);
  mod_clear (t2, m);
  mod_clear (t3, m);
}

void
pp1_stage1 (residue_t r, residue_t save, const residue_t b, const int B1, 
	    const residue_t two, const modulus_t m)
{
  if (B1 >= 50)
    {
      pp1_range_50 (r, b, two, m);
      mod_set (save, r, m);
    }
  if (B1 >= 100)
    {
      pp1_range_50_100 (r, r, two, m);
      if (mod_equal (r, two, m))
	goto end;
      mod_set (save, r, m);
    }
  if (B1 >= 150)
    {
      pp1_range_100_150 (r, r, two, m);
      if (mod_equal (r, two, m))
	  goto end;
      mod_set (save, r, m);
    }
  if (B1 >= 200)
    {
      pp1_range_150_200 (r, r, two, m);
      if (mod_equal (r, two, m))
	goto end;
      mod_set (save, r, m);
    }
  if (B1 >= 300)
    {
      pp1_range_200_300 (r, r, two, m);
      if (mod_equal (r, two, m))
	goto end;
      mod_set (save, r, m);
    }
  if (B1 >= 400)
    {
      pp1_range_300_400 (r, r, two, m);
      if (mod_equal (r, two, m))
	goto end;
      mod_set (save, r, m);
    }
  if (B1 >= 500)
    {
      pp1_range_400_500 (r, r, two, m);
      if (mod_equal (r, two, m))
	goto end;
      mod_set (save, r, m);
    }

 end:
 ;
}

unsigned long 
pp1 (const modulus_t m, const pp1_plan_t *plan)
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
  /* Faster method, but uses N of type unsigned long */
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

  pp1_stage1 (b, save, b, plan->B1, two, m);
  mod_sub (t, b, two, m);
  mod_gcd (&f, t, m);

  if (f == 1UL)
    {
      pm1_stage2 (t, b, &(plan->stage2), two, m);
      mod_gcd (&f, t, m);
    }

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
  
  /* Make bytecode for stage 1 */
  plan->B1 = B1;
  bytecoder_init (0);
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
      printf ("Byte code for stage 1: ");
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
