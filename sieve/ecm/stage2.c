#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../basicnt.h"
#include "utils.h"
#include "stage2.h"

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
      plan->pairs = NULL;
      return;
    }

  /* Choose stage 2 parameters */
  for (p = 2; p <= B2min; p = (unsigned int) getprime (p));
  /* Now p is smallest prime > B2min */
  
  /* Choose a value for d. Should depend on B2-B2min, for a start we 
     fix d=210 */
  plan->d = 210;
  
  /* List of the j values for which we need to precompute V_j(x + 1/x) */
  plan->s1 = (unsigned int) 
    eulerphi_ul ((unsigned long) plan->d) / 2;
  plan->S1 = malloc (plan->s1 * sizeof (int));
  ASSERT (plan->S1 != NULL);
  for (i = 0, j = 1; j < plan->d / 2; j += 2 /* Assumes 2|d */)
    if (gcd_ul ((unsigned long) j, (unsigned long) plan->d) == 1UL)
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
  /* For each prime B2min < p <= B2, we write p = id-j, j in S1, 
     and write into *pairs the index of j within the S_2 array. 
     When it's time to increase i, NEXT_D is written to *pairs. 
     NEXT_PASS is the signal to end stage 2. */
  
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
    malloc ((nr_primes + (B2 - B2min) / plan->d + 1) * sizeof (char));
  ASSERT (plan->pairs != NULL);
  
  /* p > B2min, so i*d > B2min - max(S_1), or i >= ceil((B2min - max(S_1)) / d).
     For now we have max(S_1) < d/2, so we can write
     i >= floor((B2min + d/2) / d)*/
  /* p <= B2, so i*d +- j <= B2, so i*d - j <= B2, so i*d <= B2 + max(S_1),
     or i <= floor((B2 + max(S_1)) / d). */

  n = 0;
  nr_pairs = 0;
  min_i = (B2min + plan->d/2) / plan->d;
  max_i = (B2 + plan->d/2) / plan->d;
  need_NEXT_D = 0;
  plan->i0 = min_i;
  plan->i1 = max_i + 1;
  for (i = min_i; i <= max_i; i++)
    {
      for (m = 0; m < plan->s1; m++)
        {
          j = plan->S1[m];
          /* See if this is a i*d +- j we need to include */
	  /* TODO: take composite id +- j into account for pairing! */
          if ((i*plan->d >= j && i*plan->d - j <= B2 && 
	       bittest (primes, i*plan->d - j)) || 
              (i*plan->d + j <= B2 && bittest (primes, i*plan->d + j)))
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
		  if (i*plan->d >= j && bittest (primes, i*plan->d - j))
		    printf ("%d ", i*plan->d - j);
		  if (i*plan->d + j <= B2 && bittest (primes, i*plan->d + j))
		    printf ("%d", i*plan->d + j);
		  printf ("\n");
		}
              ASSERT (m < 255);
              plan->pairs[n++] = (unsigned char) m;
	      nr_pairs++;
              if (i*plan->d >= j)
                bitclear (primes, i*plan->d - j);

              if (i*plan->d + j <= B2)
                bitclear (primes, i*plan->d + j);
            }
        }
      
      need_NEXT_D++;
    }
  plan->pairs[n++] = NEXT_PASS;
  if (verbose)
    printf ("Adding NEXT_PASS to list\n");
  
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
  plan->i0 = 0;
  plan->i1 = 0;
  plan->d = 0;
}
