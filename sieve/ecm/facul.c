/* Factors integers with P-1, P+1 and ECM. Input is in an mpz_t, 
   factors are unsigned long. Returns number of factors found, 
   or -1 in case of error. */

#include <stdlib.h>
#include <stdio.h>
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "modredc_ul.h"
#include "facul.h"

#define PM1_METHOD 1
#define PP1_METHOD 2
#define EC_METHOD 3


#define STATS_LEN 128

static unsigned long stats_called[STATS_LEN] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};
static unsigned long stats_found_n[STATS_LEN] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};


/* Make a simple minded strategy for factoring. We start with P-1 and
   P+1 (with x0=2/7), then an ECM curve with low bounds, then a bunch of
   ECM curves with larger bounds. How many methods to do in total is
   controlled by the n parameter: P-1, P+1 and the first ECM curve
   (with small bounds) are always done, then n ECM curves (with larger bounds)
*/

facul_strategy_t *
facul_make_strategy (const int n)
{
  facul_strategy_t *strategy;
  int i;
  strategy = malloc ((n + 4) * sizeof (facul_strategy_t));
  
  strategy[0].method = PM1_METHOD;
  strategy[0].plan = malloc (sizeof (pm1_plan_t));
  pm1_make_plan (strategy[0].plan, 315, 2205, 0);
  
  strategy[1].method = PP1_METHOD;
  strategy[1].plan = malloc (sizeof (pp1_plan_t));
  pp1_make_plan (strategy[1].plan, 525, 3255, 0);
  
  strategy[2].method = EC_METHOD;
  strategy[2].plan = malloc (sizeof (ecm_plan_t));
  ecm_make_plan (strategy[2].plan, 105, 3255, BRENT12, 10, 0);
  
  for (i = 3; i < n + 3; i++)
    {
      strategy[i].method = EC_METHOD;
      strategy[i].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (strategy[i].plan, 315, 5355, BRENT12, i + 8, 0);
    }

  strategy[n + 3].method = 0;
  strategy[n + 3].plan = NULL;

  return strategy;
}

void 
facul_clear_strategy (facul_strategy_t *strategy)
{
  int i = 0;
  for (i = 0; strategy[i].method != 0; i++)
    {
      if (strategy[i].method == PM1_METHOD)
	pm1_clear_plan (strategy[i].plan);
      else if (strategy[i].method == PP1_METHOD)
	pp1_clear_plan (strategy[i].plan);
      else if (strategy[i].method == EC_METHOD)
	ecm_clear_plan (strategy[i].plan);
      strategy[i].method = 0;
      strategy[i].plan = NULL;
    }
  free (strategy);
}

static int
cmp_ul (const unsigned long *a, const unsigned long *b)
{
  if (*a < *b) return -1;
  if (*a == *b) return 0;
  return 1;
}


void facul_print_stats (FILE *stream)
{
  int i, notfirst;
  unsigned long sum;

  fprintf (stream, "# facul statistics.\n# histogram of methods called: ");
  notfirst = 0;
  sum = 0;
  for (i = 0; i < STATS_LEN; i++)
    {
      sum += stats_called[i];
      if (stats_called[i] > 0UL)
	fprintf (stream, "%s %d: %lu", 
		 (notfirst++) ? ", " : "", i, stats_called[i]);
    }
  fprintf (stream, ". Total: %lu\n", sum);

  fprintf (stream, "# histogram of input numbers found: ");
  notfirst = 0;
  sum = 0;
  for (i = 0; i < STATS_LEN; i++)
    {
      sum += stats_found_n[i];
      if (stats_found_n[i] > 0UL)
	fprintf (stream, "%s %d: %lu", 
		 (notfirst++) ? ", " : "", i, stats_found_n[i]);
    }
  fprintf (stream, ". Total: %lu\n", sum);
}


int
facul (unsigned long *factors, const mpz_t N, facul_strategy_t *strategy)
{
  unsigned long n, f;
  modulusredcul_t m;
  residueredcul_t r;
  int i, found = 0;

  if (mpz_sgn (N) <= 0)
    return -1;
  if (mpz_cmp_ui (N, 1UL) == 0)
    return 0;
  
  /* Right now we can only deal with moduli that are 1 unsigned long in size.
     Need to write arithmetic for larger moduli. */
  if (!mpz_fits_ulong_p (N))
    return 0;

  n = mpz_get_ui (N);

  modredcul_initmod_ul (m, n);
  modredcul_init (r, m);
  
  i = 0;
  while (strategy[i].method != 0)
    {
      if (i < STATS_LEN)
	  stats_called[i]++;
      
      if (strategy[i].method == PM1_METHOD)
	f = pm1 (r, m, (pm1_plan_t *) (strategy[i].plan));
      else if (strategy[i].method == PP1_METHOD)
	f = pp1 (r, m, (pp1_plan_t *) (strategy[i].plan));
      else if (strategy[i].method == EC_METHOD)
	f = ecm (r, m, (ecm_plan_t *) (strategy[i].plan));
      else 
	{
	  /* A method value we don't know about. Something's wrong, bail out */
	  modredcul_clear (r, m);
	  modredcul_clearmod (m);
	  return -1;
	}
      
      if (f == n)
	{
	  if (i < STATS_LEN)
	    stats_found_n[i]++;
	}
      else if (f > 1)
	{
	  factors[found++] = f;
	  modredcul_clear (r, m);
	  modredcul_clearmod (m);
	  n /= f;
	  modredcul_initmod_ul (m, n);
	  modredcul_init (r, m);
	  modredcul_set_ul (r, 2UL, m);
	  if (mod_sprp (r, m))
	    {
	      modredcul_set_ul (r, 3UL, m);
	      if (mod_sprp (r, m))
		{
		  // fprintf (stderr, "facul(): %lu is a prp, exiting\n", n);
		  factors[found++] = n;
		  break;
		}
	    }
	}
      i++;
    }
  modredcul_clear (r, m);
  modredcul_clearmod (m);

  if (found > 1)
    {
      /* Sort the factors we found */
      qsort (factors, found, sizeof (unsigned long), 
	     (int (*)(const void *, const void *)) &cmp_ul);
      /* Typecasts like this are what C is all about */
    }

  return found;
}
