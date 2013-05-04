/* Factors integers with P-1, P+1 and ECM. Input is in an mpz_t, 
   factors are unsigned long. Returns number of factors found, 
   or -1 in case of error. */

#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"      /* Happens to define ULONG_BITS, which we use */
#include "portability.h"
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "facul.h"
#include "facul_doit.h"

/* These global variables are only for statistics. In case of
 * multithreaded sieving, the stats might be wrong...
 */

unsigned long stats_called[STATS_LEN] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};
unsigned long stats_found_n[STATS_LEN] = {
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
facul_make_strategy (const int n, const unsigned long fbb, 
		     const unsigned int lpb)
{
  facul_strategy_t *strategy;
  facul_method_t *methods;
  int i;
  
  strategy = malloc (sizeof (facul_strategy_t));
  strategy->lpb = 1UL << lpb;
  if (lpb == ULONG_BITS)
      strategy->lpb = ~0UL;
  /* Store fbb^2 in assume_prime_thresh */
  if (fbb > UINT32_MAX)
    strategy->assume_prime_thresh = UINT64_MAX;
  else
    strategy->assume_prime_thresh = (uint64_t) fbb * (uint64_t) fbb;

  methods = malloc ((n + 4) * sizeof (facul_method_t));
  strategy->methods = methods;

  /* run one P-1 curve with B1=315 and B2=2205 */
  methods[0].method = PM1_METHOD;
  methods[0].plan = malloc (sizeof (pm1_plan_t));
  pm1_make_plan (methods[0].plan, 315, 2205, 0);

  /* run one P+1 curve with B1=525 and B2=3255 */
  methods[1].method = PP1_METHOD;
  methods[1].plan = malloc (sizeof (pp1_plan_t));
  pp1_make_plan (methods[1].plan, 525, 3255, 0);

  /* run one ECM curve with Montgomery parametrization, B1=105, B2=3255 */
  methods[2].method = EC_METHOD;
  methods[2].plan = malloc (sizeof (ecm_plan_t));
  ecm_make_plan (methods[2].plan, 105, 3255, MONTY12, 2, 1, 0);
  
  if (n > 0)
    {
      methods[3].method = EC_METHOD;
      methods[3].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (methods[3].plan, 315, 5355, BRENT12, 11, 1, 0);
    }

  for (i = 4; i < n + 3; i++)
    {
      methods[i].method = EC_METHOD;
      methods[i].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (methods[i].plan, 315, 5355, MONTY12, i - 1, 1, 0);
    }

  methods[n + 3].method = 0;
  methods[n + 3].plan = NULL;

  return strategy;
}


void 
facul_clear_strategy (facul_strategy_t *strategy)
{
  facul_method_t *methods = strategy->methods;
  int i = 0;

  for (i = 0; methods[i].method != 0; i++)
    {
      if (methods[i].method == PM1_METHOD)
        pm1_clear_plan (methods[i].plan);
      else if (methods[i].method == PP1_METHOD)
	pp1_clear_plan (methods[i].plan);
      else if (methods[i].method == EC_METHOD)
	ecm_clear_plan (methods[i].plan);
      methods[i].method = 0;
      free (methods[i].plan);
      methods[i].plan = NULL;
    }
  free (methods);
  methods = NULL;
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
facul (unsigned long *factors, const mpz_t N, const facul_strategy_t *strategy)
{
  int found = 0;
  size_t bits;
  
#ifdef PARI
  gmp_fprintf (stderr, "%Zd", N);
#endif

  if (mpz_sgn (N) <= 0)
    return -1;
  if (mpz_cmp_ui (N, 1UL) == 0)
    return 0;
  
  /* If the composite does not fit into our modular arithmetic, return
     no factor */
  bits = mpz_sizeinbase (N, 2);
  if (bits > MODREDC2UL2_MAXBITS)
    return 0;
  
  /* Use the fastest modular arithmetic that's large enough for this input */
  if (bits <= MODREDCUL_MAXBITS)
    {
      modulusredcul_t m;
      ASSERT(mpz_fits_ulong_p(N));
      modredcul_initmod_ul (m, mpz_get_ui(N));
      found = facul_doit_ul (factors, m, strategy, 0);
      modredcul_clearmod (m);
    }
  else if (bits <= MODREDC15UL_MAXBITS)
    {
      modulusredc15ul_t m;
      unsigned long t[2];
      modintredc15ul_t n;
      size_t written;
      mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, N);
      ASSERT_ALWAYS(written <= 2);
      modredc15ul_intset_uls (n, t, written);
      modredc15ul_initmod_int (m, n);
      found = facul_doit_15ul (factors, m, strategy, 0);
      modredc15ul_clearmod (m);
    }
  else if (bits <= MODREDC2UL2_MAXBITS)
    {
      modulusredc2ul2_t m;
      unsigned long t[2];
      modintredc2ul2_t n;
      size_t written;
      mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, N);
      ASSERT_ALWAYS(written <= 2);
      modredc2ul2_intset_uls (n, t, written);
      modredc2ul2_initmod_int (m, n);
      found = facul_doit_2ul2 (factors, m, strategy, 0);
      modredc2ul2_clearmod (m);
    } else {
      modulusmpz_t m;
      modmpz_initmod_int (m, N);
      found = facul_doit_mpz (factors, m, strategy, 0);
      modmpz_clearmod (m);
    }
  
  if (found > 1)
    {
      /* Sort the factors we found */
      qsort (factors, found, sizeof (unsigned long), 
	     (int (*)(const void *, const void *)) &cmp_ul);
    }

#ifdef PARI
  if (found > 1)
    {
      fprintf (stderr, " == ");
      for (i = 0; i < found; i++)
	fprintf (stderr, "%lu%s", factors[i], 
		 (i+1 < found) ? " * " : " /* PARI */\n");
    }
  else
    fprintf (stderr, "; /* PARI */\n");
#endif

  return found;
}
