#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../basicnt.h"
#include "utils.h"
#include "stage2.h"

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
stage2_make_plan (stage2_plan_t *plan, unsigned int B2min, unsigned int B2,
		     int verbose)
{
  unsigned int p, nr_primes, nr_pairs;
  unsigned int i, min_i, max_i, j, m, n;
  unsigned char *primes;
  int need_NEXT_D;

  plan->B2 = B2;
  if (B2 <= B2min)
    {
      plan->d = 0;
      plan->s1 = 0;
      plan->S1 = NULL;
      plan->s2 = 0;
      plan->S2 = NULL;
      plan->pairs = NULL;
      return;
    }

  /* Choose stage 2 parameters */
  for (p = 2; p <= B2min; p = (unsigned int) getprime (p));
  /* Now p is smallest prime > B2min */
  
  /* Choose a value for d. Should depend on B2-B2min, for a start we 
     fix d=210 */
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
  /* For each prime B2min < p <= B2, we write p = id-(k1+k2), k1 in S1, 
     s2 in S2, and write into *pairs the index of k1 within the S_2 array. 
     The first pass has k2 = S_2[0], the second one has k2 = S_2[1], etc 
     until k2 = S_2[s2 - 1].
     When it's time to increase i, NEXT_D is written to *pairs. 
     When it's time to start a new pass, NEXT_PASS is written. This is also
     the signal to end stage 2, when s2 passes have been processed. */
  
  /* If we use s1*s2 = eulerphi(d)/2, there's exactly one way to write each
     prime as i*d +- (k1 + k2). If we allow s1*s2 >= eulerphi(d)/2, there may
     be several ways which could help pairing up primes. For now we do the
     simple way */
  
  /* Make bit array where bit at index p, prime p with B2min < p <= B2, are 
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
    malloc ((nr_primes + plan->s2 * ((B2 - B2min) / plan->d + 1)) * sizeof (char));
  ASSERT (plan->pairs != NULL);
  
  /* For now we have s2 = 1, k2 = 0 */
  
  /* p > B2min, so i*d > B2min - max(S_1), or i >= ceil((B2min - max(S_1)) / d).
     For now we have max(S_1) < d/2, so we can write
     i >= floor((B2min + d/2) / d)*/
  /* p <= B2, so i*d +- k1 <= B2, so i*d - k1 <= B2, so i*d <= B2 + max(S_1),
     or i <= floor((B2 + max(S_1)) / d). */

  n = 0;
  nr_pairs = 0;
  min_i = (B2min + plan->d/2) / plan->d;
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

  for (i = B2min + 1; i <= B2; i++)
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

void stage2_clear_plan (stage2_plan_t *plan)
{
  free (plan->pairs);
  plan->pairs = NULL;
  free (plan->S1);
  plan->S1 = NULL;
  plan->s1 = 0;
  free (plan->S2);
  plan->S2 = NULL;
  plan->s2 = 0;
  plan->d = 0;
}
