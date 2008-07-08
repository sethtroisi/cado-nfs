#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "utils.h"
#include "../basicnt.h"
#include "pm1.h"

/* #define PARI */

void
pm1_stage1 (residue_t x, const unsigned long *E, const int E_nrwords, 
            const modulus_t m)
{
  mod_pow_mp (x, x, E, E_nrwords, m);
}


unsigned long 
pm1_stage2 (residue_t r, const residue_t X, const pm1_plan_t *plan, 
	    const residue_t two, const modulus_t m)
{
  residue_t Xd, Xid, Xid1, a, t;
  residue_t *Xj;
  unsigned int k, l;
#ifdef PARI
  unsigned int id;
#endif
  
  if (plan->s2 != 1)
    {
      fprintf (stderr, "Stage 2 with more than 1 pass not implemented yet\n");
      return 1UL;
    }
  
  mod_init_noset0 (t, m);

  /* Compute Xd = V_d(X) */
  mod_V_ul (Xd, X, two, plan->d, m);

#ifdef PARI
      printf ("d = %u; Xd = V(d, X); Xd == %lu /* PARI */\n",
	      plan->d, mod_get_ul (Xd, m));
#endif
  
#define FAST_XJ_INIT
#ifndef FAST_XJ_INIT
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
    
    mod_clear (ap1_0, m);
    mod_clear (ap1_1, m);
    mod_clear (ap5_0, m);
    mod_clear (ap5_1, m);
    mod_clear (X6, m);
    mod_clear (X2, m);
  }
#endif /* if FAST_XJ_INIT */
  
  /* Do all the passes */
  mod_init_noset0 (Xid, m);
  mod_init_noset0 (Xid1, m);
  mod_init_noset0 (a, m);
  mod_set_ul (a, 1UL, m);
  l = 0;
  for (k = 0; k < plan->s2; k++)
    {
      /* For each pass, compute V_{S[k]}(X) and V_{S[k]+d}(X) so we can
	 compute the remaining V_{S[k]+id}(X) via an arithmetic progression.
	 TODO: init both with the same binary chain. */
      
      mod_V_ul (Xid, X, two, plan->S2[k], m);
      mod_V_ul (Xid1, X, two, plan->S2[k] + plan->d, m);
#ifdef PARI
      id = plan->S2[k];
      printf ("V(%u, X) == %lu /* PARI */\n",
	      plan->S2[k], mod_get_ul (Xid, m));
      printf ("V(%u, X) == %lu /* PARI */\n",
	      plan->S2[k] + plan->d, mod_get_ul (Xid1, m));
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

  if (f > 1UL)
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
  
  pm1_stage2 (t, X, plan, two, m);
  mod_gcd (&f, t, m);
  
  mod_clear (one, m);
  mod_clear (two, m);
  mod_clear (t, m);
  mod_clear (X, m);
  return f;
}



/* Looks for x in the sorted array a which has length l. Requires that
   x actually appears in a[]. */
unsigned int 
binsearch (const int *a, const unsigned int l, const int x)
{
  unsigned int low = 0, high = l - 1, mid;
  
  while (low <= high)
  {
    mid = low + (high - low) / 2;
    if (a[mid] > x)
      high = mid - 1;
    else if (a[mid] < x)
      low = mid + 1;
    else
      return mid; /* a[mid] == x */
  }
  abort(); /* Not found! */
}


static void
bitset (unsigned char *a, unsigned int i)
{
  a[i/8] |= 1 << (i%8);
}


static void
bitclear (unsigned char *a, unsigned int i)
{
  a[i/8] &= ~(1 << (i%8));
}


static int
bittest (unsigned char *a, unsigned int i)
{
  return (a[i/8] & (1 << (i%8))) != 0;
}


void 
pm1_make_plan (pm1_plan_t *plan, const unsigned int B1, const unsigned int B2,
	       int verbose)
{
  mpz_t E;
  unsigned int p, nr_primes, nr_pairs;
  unsigned int i, min_i, max_i, j, m, n;
  unsigned char *primes;
  size_t tmp_E_nrwords;
  int need_NEXT_D;
  
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

  
  /* Choose stage 2 parameters */
  
  /* Choose a value for d. Should depend on B2-B1, for a start we fix d=210 */
  plan->d = 210;
  
  /* We could do stage 2 in several passes, to reduce the cost of the 
     initialisaton of V_j(x + 1/x), 1 <= j <= d / 2, gcd(j,d) = 1.
     To do so, we need to factor this set of j into a sum of two sets,
     but this isn't possible with the condition 1 <= j <= d / 2. We may
     have to allow larger j. For now, we do only one pass. */
  plan->s2 = 1;
  plan->S2 = malloc (plan->s2 * sizeof (int));
  ASSERT (plan->S2 != NULL);
  
  /* List of the j values for which we need to precompute V_j(x + 1/x) */
  ASSERT (eulerphi_ul ((unsigned long) plan->d) % (2 * plan->s2) == 0);
  plan->s1 = (unsigned int) 
    eulerphi_ul ((unsigned long) plan->d) / (2 * plan->s2);
  plan->S1 = malloc (plan->s1 * sizeof (int));
  ASSERT (plan->S1 != NULL);
  for (i = 0, j = 1; j < plan->d / 2; j += 2 /* Assumes 2|d */)
    if (gcd_ul ((unsigned long) j, (unsigned long) plan->d) == 1)
      plan->S1[i++] = j;
  ASSERT (i == plan->s1);

  if (verbose)
    {
      printf ("S_1 = {");
      for (i = 0; i < plan->s1; i++)
	printf ("%d%s", plan->S1[i], (i + 1 < plan->s1) ? ", " : "");
      printf ("}\n");
    }
  
  /* Generate the list of pairs for stage 2 */
  /* For each prime B1 < p <= B2, we write p = id-(k1+k2), k1 in S1, s2 in S2,
     and write into *pairs the index of k1 within the S_2 array. 
     The first pass has k2 = S_2[0], the second one has k2 = S_2[1], etc 
     until k2 = S_2[s2 - 1].
     When it's time to increase i, NEXT_D is written to *pairs. 
     When it's time to start a new pass, NEXT_PASS is written. This is also
     the signal to end stage 2, when s2 passes have been processed. */
  
  /* If we use s1*s2 = eulerphi(d)/2, there's exactly one way to write each
     prime as i*d +- (k1 + k2). If we allow s1*s2 >= eulerphi(d)/2, there may
     be several ways which could help pairing up primes. For now we do the
     simple way */
  
  /* Make bit array where bit at index p, prime p with B1 < p <= B2, are 
     set to 1. */
  primes = malloc (B2 / 8 + 1);
  ASSERT (primes != NULL);
  memset (primes, 0, B2 / 8 + 1);
  nr_primes = 0;
  for ( ; p <= B2; p = getprime (p))
    {
      nr_primes++;
      bitset (primes, p);
    }
  
  /* We need at most one pair per prime, plus the number of NEXT_D and 
     NEXT_PASS codes */
  plan->pairs = 
    malloc ((nr_primes + plan->s2 * ((B2 - B1) / plan->d + 1)) * sizeof (char));
  ASSERT (plan->pairs != NULL);
  
  /* For now we have s2 = 1, k2 = 0 */
  
  /* p > B1, so i*d > B1 - max(S_1), or i >= ceil((B1 - max(S_1)) / d).
     For now we have max(S_1) < d/2, so we can write
     i >= floor((B1 + d/2) / d)*/
  /* p <= B2, so i*d +- k1 <= B2, so i*d - k1 <= B2, so i*d <= B2 + max(S_1),
     or i <= floor((B2 + max(S_1)) / d). */

  n = 0;
  nr_pairs = 0;
  min_i = (B1 + plan->d/2) / plan->d;
  max_i = (B2 + plan->d/2) / plan->d;
  need_NEXT_D = 0;
  for (i = min_i; i <= max_i; i++)
    {
      plan->S2[0] = min_i * plan->d;
      for (m = 0; m < plan->s1; m++)
        {
          unsigned int k1 = plan->S1[m];
          /* See if this is a i*d +- j we need to include */
	  /* TODO: take composite id+-j into account for pairing! */
          if ((i*plan->d >= k1 && i*plan->d - k1 <= B2 && 
	       bittest (primes, i*plan->d - k1)) || 
              (i*plan->d + k1 <= B2 && bittest (primes, i*plan->d + k1)))
            {
	      while (need_NEXT_D)
		{
		  plan->pairs[n++] = NEXT_D;
		  if (verbose)
		    printf ("Adding NEXT_D to list\n");
		  need_NEXT_D--;
		}
	      
	      if (verbose)
		{
		  printf ("Adding %d*d +- %d (=S1[%d]) to list, includes "
			  "primes ", i, k1, m);
		  if (i*plan->d >= k1 && bittest (primes, i*plan->d - k1))
		    printf ("%d ", i*plan->d - k1);
		  if (i*plan->d + k1 <= B2 && bittest (primes, i*plan->d + k1))
		    printf ("%d", i*plan->d + k1);
		  printf ("\n");
		}
              ASSERT (m < 255);
              plan->pairs[n++] = (unsigned char) m;
	      nr_pairs++;
              if (i*plan->d >= k1)
                bitclear (primes, i*plan->d - k1);

              if (i*plan->d + k1 <= B2)
                bitclear (primes, i*plan->d + k1);
            }
        }
      
      need_NEXT_D++;
    }
  plan->pairs[n++] = NEXT_PASS;
  if (verbose)
    printf ("Adding NEXT_PASS to list\n");
  
  if (verbose)
    {
      printf ("S_2 = {");
      for (i = 0; i < plan->s2; i++)
	printf ("%d%s", plan->S2[i], (i + 1 < plan->s2) ? ", " : "");
      printf ("}\n");
    }

  if (verbose)
    {
      printf ("pairs = ");
      for (i = 0; i < n; i++)
	{
	  if (plan->pairs[i] == NEXT_D)
	    printf ("NEXT_D ");
	  else if (plan->pairs[i] == NEXT_PASS)
	    printf ("NEXT_PASS ");
	  else
	    printf ("%d ", plan->pairs[i]);
	}
      printf ("\n");
    }
  
  if (verbose)
    printf ("Used %u pairs to include %u primes, avg. %.2f primes/pair\n", 
	    nr_pairs, nr_primes, (double)nr_primes / (double) nr_pairs);

  for (i = B1 + 1; i <= B2; i++)
    {
      if (bittest (primes, i))
	{
	  fprintf (stderr, "Error, prime %d is still set in bitfield\n", i);
	  abort ();
	}
    }

  free (primes);
  getprime (0);
}


void 
pm1_clear_plan (pm1_plan_t *plan)
{
  free (plan->pairs);
  plan->pairs = NULL;
  
  free (plan->S1);
  plan->S1 = NULL;
  plan->s1 = 0;
  free (plan->S2);
  plan->S2 = NULL;
  plan->s2 = 0;

  free (plan->E);
  plan->E = NULL;
  plan->E_nrwords = 0;
  plan->B1 = 0;
  plan->d = 0;
}
