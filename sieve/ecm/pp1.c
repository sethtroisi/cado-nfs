#include <stdio.h>
#include <stdlib.h>
#include "pp1.h"
#if defined(MODREDCUL)
#include "modredc_ul.h"
#include "modredc_ul_default.h"
#define pp1 pp1_ul
#define pp1_stage2 pp1_stage2_ul
#elif defined(MODREDC15UL)
#include "modredc_15ul.h"
#include "modredc_15ul_default.h"
#define pp1 pp1_15ul
#define pp1_stage2 pp1_stage2_15ul
#else
#error Please define MODREDCUL or MODREDC15UL
#endif

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
  ASSERT (r != two);
  mod_mul (r, a, a, m);
  mod_sub (r, r, two, m);
}

static void
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
            mod_set (t2, A, m);
            pp1_add (A, A, t, t2, m);
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

void 
pp1_stage2 (residue_t r, const residue_t X, const stage2_plan_t *plan, 
	    const residue_t two, const modulus_t m)
{
  residue_t Xd, Xid, Xid1, a, a_bk, t;
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
    mod_set (ap1_0, X, m);         /* ap1_0 = V_1(X) = X */
    pp1_double (X2, X, two, m);    /* X2 = V_2(X) = X^2 - 2 */
    pp1_add (X6, X2, X, X, m);     /* V_3(X) = V_2(X) * V_1(X) - V_1(X) */
    pp1_add (ap5_0, X6, X2, X, m); /* V_5(X) = V_3(X) * V_2(X) - V_1(X) */
    pp1_double (X6, X6, two, m);   /* V_6(X) = V_3(X)*V_3(X) - 2 */
    pp1_add (ap1_1, X6, X, ap5_0, m); /* V_7(X) = V_6(X) * V_1(X) - V_5(X) */
    pp1_add (ap5_1, X6, ap5_0, X, m); /* V_11(X) = V_6(X) * V_5(X) - V_1(X) */
    
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
	
        pp1_add (t, ap1_1, X6, ap1_0, m);
        mod_set (ap1_0, ap1_1, m);
        mod_set (ap1_1, t, m);
        i1 += 6;
	
        pp1_add (t, ap5_1, X6, ap5_0, m);
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
  mod_init_noset0 (a_bk, m);
  mod_set1 (a, m);
  mod_set (a_bk, a, m);
  l = 0;
  
  {
    /* Compute V_{i0 * d}(X) and V_{(i0 + 1) * d}(X) so we can
       compute the remaining V_{id}(X) via an arithmetic progression.
       TODO: init both with the same binary chain. */
    
    /* Todo: do both with the same addition chain */
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
	
	/* See if we got a == 0. If yes, restore previous a value and
	   end stage 2 */
	if (mod_is0 (a, m))
	  {
	    mod_set (a, a_bk, m);
	    break;
	  }
	mod_set (a_bk, a, m); /* Save new a value */
	
	/* Advance i by 1 */
	if (plan->pairs[l] == NEXT_D)
	  {
	    pp1_add (t, Xid1, Xd, Xid, m);
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
  mod_clear (a_bk, m);
  mod_clear (t, m);
}


int  
pp1 (modint_t f, const modulus_t m, const pp1_plan_t *plan)
{
  residue_t b, two, save, t;

  mod_init_noset0 (b, m);
  mod_init_noset0 (two, m);
  mod_init_noset0 (save, m);
  mod_init_noset0 (t, m);
  mod_set1 (two, m);
  mod_add (two, two, two, m);
  
  /* Compute 2/7 (mod N) */
  mod_set (b, two, m);
  mod_div7 (b, b, m);
  
  pp1_stage1 (b, plan->bc, plan->bc_len, two, m);
  mod_sub (t, b, two, m);
  mod_gcd (f, t, m);

  if (mod_intcmp_ul(f, 1UL) == 0 && plan->stage2.B2 > plan->B1)
    {
      pp1_stage2 (t, b, &(plan->stage2), two, m);
      mod_gcd (f, t, m);
    }
  
  mod_clear (b, m);
  mod_clear (save, m);
  mod_clear (two, m);
  mod_clear (t, m);

  return 0;
}
