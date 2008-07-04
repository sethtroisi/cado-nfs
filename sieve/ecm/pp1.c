#include "mod_ul.h"
#include "pp1.h"

#define dup(a,b) mod_mulredc(a,b,b,invm,m);mod_sub(a,a,two,m);
#define add(a,b,c,d) mod_mulredc(a,b,c,invm,m);mod_sub(a,a,d,m);
/* Safe version of add that can be used if "a" and "d" are identical */
#define adds(a,b,c,d) mod_set(t3,d,m);mod_mulredc(a,b,c,invm,m);mod_sub(a,a,t3,m);
#define set(a,b) mod_set(a,b,m);
#define swap(a,b) mod_swap(a,b,m);

static void
pp1_range_50 (residue_t r, const residue_t b, const residue_t two, 
              const unsigned long invm, const modulus_t m)
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
                  const unsigned long invm, const residue_t m)
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
                   const unsigned long invm, const residue_t m)
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
                   const unsigned long invm, const residue_t m)
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
                   const unsigned long invm, const residue_t m)
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
                   const unsigned long invm, const residue_t m)
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
                   const unsigned long invm, const residue_t m)
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
	    const residue_t two, const unsigned long invm, const residue_t m)
{
  if (B1 >= 50)
    {
      pp1_range_50 (r, b, two, invm, m);
      mod_set (save, r, m);
    }
  if (B1 >= 100)
    {
      pp1_range_50_100 (r, r, two, invm, m);
      if (mod_equal (r, two, m))
	goto end;
      mod_set (save, r, m);
    }
  if (B1 >= 150)
    {
      pp1_range_100_150 (r, r, two, invm, m);
      if (mod_equal (r, two, m))
	  goto end;
      mod_set (save, r, m);
    }
  if (B1 >= 200)
    {
      pp1_range_150_200 (r, r, two, invm, m);
      if (mod_equal (r, two, m))
	goto end;
      mod_set (save, r, m);
    }
  if (B1 >= 300)
    {
      pp1_range_200_300 (r, r, two, invm, m);
      if (mod_equal (r, two, m))
	goto end;
      mod_set (save, r, m);
    }
  if (B1 >= 400)
    {
      pp1_range_300_400 (r, r, two, invm, m);
      if (mod_equal (r, two, m))
	goto end;
      mod_set (save, r, m);
    }
  if (B1 >= 500)
    {
      pp1_range_400_500 (r, r, two, invm, m);
      if (mod_equal (r, two, m))
	goto end;
      mod_set (save, r, m);
    }

 end:
 ;
}

void
pp1 (residue_t f, const modulus_t m, const unsigned long invm, const int B1, 
     const int B2)
{
  residue_t b, two, save;

  mod_init_noset0 (b, m);
  mod_init_noset0 (two, m);
  mod_init_noset0 (save, m);
  mod_set_ul_reduced (two, 2UL, m);
  mod_tomontgomery (two, two, m);
  
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
  mod_tomontgomery (b, b, m);
  mod_add (b, b, b, m);
#endif

  pp1_stage1 (b, save, b, B1, two, invm, m);
  mod_sub (b, b, two, m);
  mod_gcd (f, b, m);
  mod_add (b, b, two, m);

  mod_clear (b, m);
  mod_clear (save, m);
  mod_clear (two, m);
}

