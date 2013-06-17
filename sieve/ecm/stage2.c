/* Code to produce a stage 2 plan for P-1, P+1, and ECM.
   For given B2min, B2 determines which value of d and pairs (i,j) to use
   so that the id+-j values cover the primes in ]B2min, B2] */

#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "stage2.h"
#include "portability.h"

#if 0
/* Looks for x in the sorted array a which has length l. Requires that
   x actually appears in a[]. */
static unsigned int 
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
#endif


static inline unsigned long 
eulerphi_ul (unsigned long n)
{
  unsigned long p, r = 1UL;
  
  if (n == 0UL) /* Undefined, we return 0 */
    return 0UL;

  if (n % 2UL == 0UL)
    {
      n /= 2UL;
      while (n % 2UL == 0UL)
        {
          n /= 2UL;
          r *= 2UL;
        }
    }

  for (p = 3UL; p*p <= n; p += 2UL)
    {
      if (n % p == 0UL)
        {
          n /= p;
          r *= p - 1UL;
          while (n % p == 0UL)
            {
              n /= p;
              r *= p;
            }
        }
    }
  /* Now n is either 1 or a prime */
  if (n > 1UL)
    r *= n - 1UL;
  
  return r;
}

/* Tries to find a factor p of n such that primes[p] != 0.
   We assume p is the largest prime factor of n.
   If successful, return p, otherwise return 0. */
unsigned int 
find_prime_div (const unsigned int n, const unsigned int min_f, 
		const unsigned char *primes)
{
    unsigned int p = n, f = min_f;

    ASSERT (n % 2 != 0);
    while (p >= f*f)
    {
	while (p % f == 0)
	{
	    p /= f;
 	    if (primes[p])
		return p;
	}
	f += 2; /* Assumes n even */
    }
    
    if (primes[p])
	return p;
    
    return 0;
}

void 
stage2_make_plan (stage2_plan_t *plan, const unsigned int B2min, 
                  const unsigned int B2, const int verbose)
{
  unsigned int p, nr_primes, nr_pairs;
  unsigned int i, min_i, max_i, j, d, m, n, not_in_d;
  unsigned char *primes;
  int need_NEXT_D;
  const int composite_pairing = 1;

  plan->B2 = B2;
  if (B2 <= B2min)
    {
      plan->d = 0;
      plan->s1 = 0;
      plan->S1 = NULL;
      plan->pairs = NULL;
      return;
    }

  /* Choose stage 2 parameters */


  /* Choose a value for d. Should depend on B2-B2min, for a start we 
     fix d=210 */
  d = 210;
  plan->d = d;
  not_in_d = 11; /* the smallest prime not in d */
  
  /* List of the j values for which we need to precompute V_j(x + 1/x) */
  plan->s1 = (unsigned int) 
    eulerphi_ul ((unsigned long) d) / 2;
  plan->S1 = malloc (plan->s1 * sizeof (int));
  ASSERT (plan->S1 != NULL);
  for (i = 0, j = 1; j < d / 2; j += 2 /* Assumes 2|d */)
    if (gcd_ul ((unsigned long) j, (unsigned long) d) == 1UL)
      plan->S1[i++] = j;
  ASSERT (i == plan->s1);

  if (verbose)
    {
      printf ("S_1 = {");
      for (i = 0; i < plan->s1; i++)
	printf ("%d%s", plan->S1[i], (i + 1 < plan->s1) ? ", " : "");
      printf ("}\n");
    }
  
  /* Preliminary choice for the smallest and largest i value we might need.
     The smallest may increase yet due to pairing with composite values. */
  /* p > B2min, so i*d > B2min - max(S_1), or i >= ceil((B2min - max(S_1)) / d).
     For now we have max(S_1) < d/2, so we can write
     i >= floor((B2min + d/2) / d)*/
  /* p <= B2, so i*d +- j <= B2, so i*d - j <= B2, so i*d <= B2 + max(S_1),
     or i <= floor((B2 + max(S_1)) / d). */
  min_i = (B2min + d/2) / d;
  max_i = (B2 + d/2) / d;
  if (verbose)
      printf ("Initial choice for min_i = %u, max_i = %u\n", min_i, max_i);

  /* Generate the list of pairs for stage 2 */
  /* For each prime B2min < p <= B2, we write p = id-j, j in S1, 
     and write into *pairs the index of j within the S_2 array. 
     When it's time to increase i, NEXT_D is written to *pairs. 
     NEXT_PASS is the signal to end stage 2. */
  
  /* Make array where entry at index p, prime p with B2min < p <= B2, are 
     set to 1. The largest value we can write as i*d+j with i < max_i and 
     j < d/2 is max_i * d + d/2 - 1 */
  primes = malloc (max_i * d + d/2);
  ASSERT (primes != NULL);
  memset (primes, 0, max_i * d + d/2);
  nr_primes = 0;

  for (p = 2; p <= B2min; p = (unsigned int) getprime (p));
  /* Now p is smallest prime > B2min */
  for ( ; p <= B2; p = getprime (p))
    {
      nr_primes++;
      primes[p] = 1;
    }
  
  /* We need at most one pair per prime, plus the number of NEXT_D and 
     NEXT_PASS codes */
  plan->pairs = 
    malloc ((nr_primes + (B2 - B2min) / d + 1) * sizeof (char));
  ASSERT (plan->pairs != NULL);
  

  /* Lower max_i so that max_i * d +- j, 0 < j < d/2, actually includes 
     any primes */
  for (i = max_i; i >= min_i; i--)
    {
      for (m = 0; m < plan->s1; m++)
        {
	  j = plan->S1[m];
          if ((i*d >= j && primes[i*d - j]) || primes[i*d + j])
	    {
	      if (verbose)
	        {
		  printf("Final max_i = %u, found prime %u or %u\n",
			 i, i*d - j, i*d + j);
	        }
	      break;
	    }
	}
      if (m < plan->s1)
	break;
    }
  max_i = i;

  n = 0;
  nr_pairs = 0;
  need_NEXT_D = 0;
  
  if (composite_pairing)
  {
      /* Do a pass over the primes in reverse, flagging off those primes
	 that are included as composite values id+-j, where id-+j is prime.
	 The smallest prime not in d is not_in_d, so any proper divisor of
	 i*d+-j is at most (i*d+d/2+1) / not_in_d */
      for (p = B2; p >= (max_i * d + d/2 - 1) / not_in_d; p--)
      {
	  if (primes[p])
	  {
	      unsigned int q, r;
	      int sj;
	      /* p = id + j or id - j with 0 < j < d/2 */
	      j = p % d;
	      sj = (int) j - ((j > d/2) ? (int) d : 0);
	      /* p = i*d + sj, q = i*d - sj = p - 2*sj */
	      if ((int)p < 2*sj)
		  continue;
	      q = (int) p - 2*sj;	
	      
	      if (!primes[q])      /* If q is not a prime, see if a stage 2  
				      prime divides q */
	      {
		  r = find_prime_div (q, not_in_d, primes);
		  if (r)
		  {
		      if (verbose)
			  printf ("Including %u as factor of %u = %u * %u "
				  "+ %d which pairs with prime %u = %u"
				  " * %u + %d\n", 
				  r, q, (p - sj) / d, d, -sj, p, 
				  (p-sj)/d, d, sj);
		      primes[r] = 0;
		  }
	      }
	  }
      }
     
      /* Some small primes may remain. Go through all possible (i,j)
	 values in order of decreasing id+j, and choose those that
	 cover two small primes */
      for (i = max_i; i >= min_i && i > 0; i--)
      {
	  for (m = 0; m < plan->s1; m++)
	  {
	      unsigned int q;
	      j = plan->S1[m];
	      p = i * d - j;
	      q = i * d + j;
	      if (!primes[q] && !primes[q])
	      {
		  unsigned int r1, r2;
		  r1 = find_prime_div (p, not_in_d, primes);
		  r2 = find_prime_div (q, not_in_d, primes);
		  if (r1 && r2)
		  {
		      primes[p] = 1; /* Flag on one of the two "primes" that
					include these two composite values */
		      primes[r1] = 0;
		      primes[r2] = 0;
		      if (verbose)
			  printf ("Including %u * %u +- %u which includes "
				  "%u and %u as factors\n", i, d, j, r1, r2);
		  }
	      }
	  }
      }

      /* Some small primes may still remain. Go through all possible (i,j)
	 values in order of decreasing id+j again, and choose those that
	 cover at least one prime */
     for (i = max_i; i >= min_i && i > 0; i--)
      {
	  for (m = 0; m < plan->s1; m++)
	  {
	      unsigned int q;
	      j = plan->S1[m];
	      p = i * d - j;
	      q = i * d + j;
	      if (!primes[p] && !primes[q])
	      {
		  unsigned int r1, r2;
		  r1 = find_prime_div (p, not_in_d, primes);
		  r2 = find_prime_div (q, not_in_d, primes);
		  if (r1 || r2)
		  {
		      primes[p] = 1; /* Flag on one of the two "primes" that
					include this composite value */
		      primes[r1] = 0;
		      primes[r2] = 0;
		      if (verbose)
			  printf ("Including %u * %u +- %u which includes "
				  "%u as factors\n", i, d, j, 
				  (r1 == 0) ? r2 : r1);
		  }
	      }
	  }
      }
  }
  
  /* Increase min_i so that min_i * d +- j, 0 < j < d/2, actually includes 
     any primes */
  for (i = min_i; i <= max_i; i++)
    {
      for (m = 0; m < plan->s1; m++)
        {
	  j = plan->S1[m];
          if ((i*d >= j && primes[i*d - j]) || primes[i*d + j])
	    {
	      if (verbose)
	        {
		  printf("Final min_i = %u, found prime %u or %u\n",
			 i, i*d - j, i*d + j);
	        }
	      break;
	    }
	}
      if (m < plan->s1)
	break;
    }
  min_i = i;

  
  /* For the remaining primes, write the required (i,j)-pairs to a list */

  for (i = min_i; i <= max_i; i++)
    {
      for (m = 0; m < plan->s1; m++)
        {
          j = plan->S1[m];
          /* See if this is a i*d +- j we need to include */
          if ((i*d >= j && primes[i*d - j]) || primes[i*d + j])
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
			  "primes ", i, j, m);
		  if (i*d >= j && primes[i*d - j])
		    printf ("%d ", i*d - j);
		  if (i*d + j <= B2 && primes[i*d + j])
		    printf ("%d", i*d + j);
		  printf ("\n");
		}
              ASSERT (m < 255);
              plan->pairs[n++] = (unsigned char) m;
	      nr_pairs++;
              if (i*d >= j)
                primes[i*d - j] = 0;

              if (i*d + j <= B2)
                primes[i*d + j] = 0;
            }
        }
      
      need_NEXT_D++;
    }
  plan->pairs[n++] = NEXT_PASS;
  if (verbose)
    printf ("Adding NEXT_PASS to list\n");
  plan->i0 = min_i;
  plan->i1 = max_i + 1;
  
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

  for (i = B2min + 1; i <= B2; i++)
    {
      if (primes[i])
	{
	  fprintf (stderr, "Error, prime %d is still set\n", i);
	  abort ();
	}
    }

  free (primes);
  getprime (0);
}

void stage2_clear_plan (stage2_plan_t *plan)
{
  free (plan->pairs);
  plan->pairs = NULL;
  free (plan->S1);
  plan->S1 = NULL;
  plan->s1 = 0;
  plan->i0 = 0;
  plan->i1 = 0;
  plan->d = 0;
}
