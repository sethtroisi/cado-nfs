/* Factors integers with P-1, P+1 and ECM. Input is in an mpz_t, 
   factors are unsigned long. Returns number of factors found, 
   or -1 in case of error. */

#include "cado.h"
#include <stdint.h>	/* AIX wants it first (it's a bug) */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <regex.h>

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

int
nb_curves (const unsigned int lpb)
{
  /* the following table, computed with the proba_cofactor() function in the
     facul.sage file, ensures a probability of at least about 90% to find a
     factor below 2^lpb with n = T[lpb] */
#define LPB_MAX 33
  int T[LPB_MAX+1] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 0-9 */
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 10-19 */
                      0, 0, 1 /*22:0.9074*/, 2 /*23:0.9059*/, 3 /*24:0.8990*/,
                      5 /*25:0.9194*/, 6 /*26:0.9065*/, 8 /*27:0.9053*/,
                      10 /*28:0.9010*/, 13 /*29:0.9091*/, 16 /*30:0.9134*/,
                      18 /*31:0.9039*/, 21 /*32:0.9076*/, 24/*33:0.8963*/};
  return (lpb <= LPB_MAX) ? T[lpb] : T[LPB_MAX];
#undef LPB_MAX
}

// ncurves = 100 corresponds to looking for 50-bit factors.
// ncurves = 200 corresponds to looking for 64-bit factors.


/* Make a simple minded strategy for factoring. We start with P-1 and
   P+1 (with x0=2/7), then an ECM curve with low bounds, then a bunch of
   ECM curves with larger bounds. How many methods to do in total is
   controlled by the n parameter: P-1, P+1 and the first ECM curve
   (with small bounds) are always done, then n ECM curves (with larger bounds)
*/

facul_strategy_t *
facul_make_strategy (const unsigned long fbb, const unsigned int lpb,
                     int n, const int verbose)
{
  facul_strategy_t *strategy;
  facul_method_t *methods;
  int i;

  if (n == 0)
    n = nb_curves (lpb);
  ASSERT_ALWAYS(n < STATS_LEN);
  strategy = malloc (sizeof (facul_strategy_t));
  strategy->lpb = lpb;
  /* Store fbb^2 in assume_prime_thresh */
  strategy->assume_prime_thresh = (double) fbb * (double) fbb;
  strategy->BBB = (double) fbb * strategy->assume_prime_thresh;

  methods = malloc ((n + 4) * sizeof (facul_method_t));
  strategy->methods = methods;

  /* run one P-1 curve with B1=315 and B2=2205 */
  methods[0].method = PM1_METHOD;
  methods[0].plan = malloc (sizeof (pm1_plan_t));
  pm1_make_plan (methods[0].plan, 315, 2205, verbose);

  /* run one P+1 curve with B1=525 and B2=3255 */
  methods[1].method = PP1_27_METHOD;
  methods[1].plan = malloc (sizeof (pp1_plan_t));
  pp1_make_plan (methods[1].plan, 525, 3255, verbose);

  /* run one ECM curve with Montgomery parametrization, B1=105, B2=3255 */
  methods[2].method = EC_METHOD;
  methods[2].plan = malloc (sizeof (ecm_plan_t));
  ecm_make_plan (methods[2].plan, 105, 3255, MONTY12, 2, 1, verbose);

  if (n > 0)
    {
      methods[3].method = EC_METHOD;
      methods[3].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (methods[3].plan, 315, 5355, BRENT12, 11, 1, verbose);
    }

  /* heuristic strategy where B1 is increased by sqrt(B1) at each curve */
  double B1 = 105.0;
  for (i = 4; i < n + 3; i++)
    {
      double B2;
      unsigned int k;

      B1 += sqrt (B1);
      B2 = 17.0 * B1;
      /* we round B2 to (2k+1)*105, thus k is the integer nearest to
         B2/210-0.5 */
      k = B2 / 210.0;
      methods[i].method = EC_METHOD;
      methods[i].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (methods[i].plan, (unsigned int) B1, (2 * k + 1) * 105,
                     MONTY12, i - 1, 1, 0);
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
      else if (methods[i].method == PP1_27_METHOD)
	pp1_clear_plan (methods[i].plan);
      else if (methods[i].method == PP1_65_METHOD)
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
  if (bits > MODMPZ_MAXBITS)
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
    } 
  else 
    {
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


/*****************************************************************************/
/*                       STRATEGY BOOK                                       */
/*****************************************************************************/

static int
get_index_plan (precompute_plan_t*tab, int len, int B1, int B2, int method)
{
  for (int i = 0; i < len; i++)
    if (tab[i].B1 == B1 && tab[i].B2 == B2 && tab[i].method == method)
      return i;
  return -1;
}


static void
return_data_ex (char** res, regmatch_t *pmatch, size_t nmatch, const char * str_process)
{
  //(re)init res
  for (size_t i = 0; i < nmatch;i++)
    res[i] = NULL;

  if ( pmatch[0].rm_so != pmatch[0].rm_eo)
    {
      for (size_t i = 1; i < nmatch; i++)
	{
	  int start = pmatch[i].rm_so;
	  int end = pmatch[i].rm_eo;
	  if (start == -1)
	    {
	      break;
	    }
	  else
	    {
	      int size = end-start;
	      char* el = malloc (sizeof(char)*(size+1));
	      assert (el != NULL);
	      strncpy (el, &str_process[start], size);
	      el[size] = '\0';
	      res[i-1] = el;
	    }
	}
    }
}


static int
process_line (facul_strategies_t* strategies, unsigned int* index_st,
	      int* length_plan, const unsigned int rfbb, const unsigned int afbb,
	      const unsigned long rlpb, const unsigned long alpb, const char *str,
	      const int verbose)
{
  int index_method = 0; /* this is the index of the current factoring
			   methods */

  regex_t preg_index, preg_fm;
  //regular expression for the sides
  const char *str_preg_index = "r0=([[:digit:]]+),[[:space:]]r1=([[:digit:]]+)";
  //regular expression for the strategy
  const char *str_preg_fm = "([[:alnum:]]{2,3}):[[:space:]](...-?[[:alnum:]]{2,3}?),[[:space:]]([[:digit:]]+),[[:space:]]([[:digit:]]+)";
  regcomp (&preg_index, str_preg_index, REG_ICASE|REG_EXTENDED);
  regcomp (&preg_fm, str_preg_fm, REG_ICASE|REG_EXTENDED);
  
  //process the ligne
  const char * str_process = &str[0];

  while (str_process[0] != '\0' )
    {
      //init
      size_t nmatch = 5;
      regmatch_t *pmatch= calloc (sizeof(*pmatch), nmatch);
      char **res = malloc (sizeof(char*) * nmatch);
      /*TEST REGULAR EXPRESSION  'preg_index*/
      regexec (&preg_index, str_process, nmatch, pmatch, 0);
      return_data_ex (res, pmatch, nmatch, str_process);
      if (res[0] != NULL)
	{
	  /* changes the current strategy and allocates it! */
	  facul_both_strategy_t strategy;

	  strategy.lpb[0] = rlpb;
	  strategy.lpb[1] = alpb;
	  /* Store fbb^2 in assume_prime_thresh */
	  strategy.assume_prime_thresh[0] = (double) rfbb * (double) rfbb;
	  strategy.assume_prime_thresh[1] = (double) afbb * (double) afbb;

	  strategy.BBB[0] = (double) rfbb * strategy.assume_prime_thresh[0];
	  strategy.BBB[1] = (double) afbb * strategy.assume_prime_thresh[1];

	  facul_method_t *methods;
	  methods = malloc (NB_MAX_METHODS * sizeof (facul_method_t));
	  ASSERT_ALWAYS (methods != NULL);
	  strategy.methods = methods;
	  strategy.methods[0].method = 0;
	  strategy.methods[0].plan = NULL;

	  index_st[0] = atoi(res[0]);
	  index_st[1] = atoi(res[1]);
	  strategies->strategy[ atoi(res[0]) ][ atoi(res[1]) ] = strategy;
	}

      /*else TEST REGULAR EXPRESSION  'preg_alg'*/
      else 
	{
	  regexec (&preg_fm, str_process, nmatch, pmatch, 0);
	  return_data_ex (res, pmatch, nmatch, str_process);
	  if (res[0] != NULL)
	    {
	      /*add the new factoring method to the current strategy! */
	      facul_both_strategy_t* strategy =
		&strategies->strategy[index_st[0]][index_st[1]];
	      facul_method_t* methods= strategy->methods;
	      ASSERT_ALWAYS (methods!=NULL);

	      int side;
	      //todo: remove the words ALG and RAT!!!//S0 or S1
	      if (strcmp (res[0], "ALG") == 0 || strcmp (res[0], "S1") == 0)
		side = 1;
	      else if (strcmp (res[0], "RAT") == 0 || strcmp(res[0], "S0") == 0)
		side = 0;
	      else 
		{
		  side = atoi(res[0]);
		}
	      int B1 = atoi (res[2]);
	      int B2 = atoi (res[3]);

	      if ( B1 == 0 || B2 == 0)//This is the zero strategy
	      	{
		  methods[index_method].method = PM1_METHOD;
	      	  methods[index_method].side = side;
		  methods[index_method].plan = NULL;
	      	}
	      else
		{
		  int curve = 0;
		  if (strcmp (res[1], "PM1") == 0)
		    methods[index_method].method = PM1_METHOD;
		  else if (strcmp (res[1], "PP1-27") == 0)
		    methods[index_method].method = PP1_27_METHOD;
		  else if (strcmp (res[1], "PP1-65") == 0)
		    methods[index_method].method = PP1_65_METHOD;
		  else 
		    {
		      methods[index_method].method = EC_METHOD;
		      if (strcmp (res[1], "ECM-B12") == 0)
			curve = BRENT12;
		      else if (strcmp (res[1], "ECM-M12") == 0)
			curve = MONTY12;
		      else if (strcmp (res[1], "ECM-M16") == 0)
			curve = MONTY16;
		      else
			{
			  fprintf (stderr,
				   "error : the method '%s' is unknown!\n",
				   res[1]);
			  return -1;
			}
		    }
		  //search if the plan isn't already computed!
		  methods[index_method].side = side;
		  int index_plan = get_index_plan(strategies->plan,
						  *length_plan, B1, B2,
						  methods[index_method].method);
		  //index_plan = -1;
		  if (index_plan == -1)
		    {
		      void* plan = NULL;
		      if (methods[index_method].method == PM1_METHOD)
		      	{
			  plan = malloc (sizeof (pm1_plan_t));
			  pm1_make_plan (plan, B1, B2, verbose);
			}
		      else if (methods[index_method].method == PP1_27_METHOD ||
			       methods[index_method].method == PP1_65_METHOD)
		      	{
			  plan = malloc (sizeof (pp1_plan_t));
			  pp1_make_plan (plan, B1, B2, verbose);
			}
		      else {
		      	//int param =
			//(curve==MONTY16)?1:3+index_method; todo:
			//change -1 by +3 (-1 is to do the same thing
			//that make_strategy)
			//todo: choose better values for param.
			int param = (curve==MONTY16)?1:index_method-2;
			if (param <= 1)
			  param = index_method +10;//for example
			if (index_method == 3)
			  param = 2;//first  M12
			if (curve == BRENT12)
			  param = 11;//B12
			plan = malloc (sizeof (ecm_plan_t));

		      	ecm_make_plan (plan,
		      		       B1,B2, curve, param, 1, verbose);
		      }
		      strategies->plan[*length_plan].method = methods[index_method].method;
		      strategies->plan[*length_plan].B1 = B1;
		      strategies->plan[*length_plan].B2 = B2;
		      strategies->plan[*length_plan].plan = plan;

		      index_plan = *length_plan;
		      (*length_plan)++;
		      ASSERT (*length_plan < NB_MAX_METHODS);
		      //to show the end of plan!
		      strategies->plan[*length_plan].plan = NULL;
		      strategies->plan[*length_plan].B1 = 0;
		      strategies->plan[*length_plan].B2 = 0;
		      strategies->plan[*length_plan].method = 0;
		    }

		  methods[index_method].plan = strategies->plan[index_plan].plan;

		}
	      index_method++;
	      ASSERT_ALWAYS (index_method < NB_MAX_METHODS);

	      //to show the end of methods!
	      methods[index_method].plan = NULL;
	      methods[index_method].method = 0;
	    }
	  else//to end the while
	    {
	      pmatch[0].rm_eo = strlen(str_process);
	    }
	}
      str_process = &str_process[pmatch[0].rm_eo];
      //free
      for (size_t i = 0; i < nmatch; i++)
	free(res[i]);
      free(res);
      free (pmatch);
    }
  
  regfree(&preg_index);
  regfree(&preg_fm);
  return 0;
}


facul_strategies_t*
facul_make_strategies(const unsigned long rfbb, const unsigned int rlpb,
		      const unsigned int rmfb, const unsigned long afbb,
		      const unsigned int alpb, const unsigned int amfb,
		      FILE* file, const int verbose)
{
  if (file == NULL)
    return NULL;
  
  //printf ("create strategies!\n");
  facul_strategies_t* strategies = malloc (sizeof(facul_strategies_t));
  ASSERT (strategies != NULL);
  strategies->mfb[0] = rmfb;
  strategies->mfb[1] = amfb;
  facul_both_strategy_t** strategy = malloc (sizeof (*strategy) * (rmfb+1));
  ASSERT (strategy != NULL);
  unsigned int r, a;
  for (r = 0; r <= rmfb; r++)
    {
      strategy[r] = malloc (sizeof (*strategy[r]) * (amfb+1));
      ASSERT (strategy[r] != NULL);
    } 

  /* to precompute our method */
  precompute_plan_t* plan = malloc (sizeof(precompute_plan_t*)*NB_MAX_METHODS);
  plan[0].method = 0;
  plan[0].plan = NULL;

  strategies->plan = plan;
  strategies->strategy = strategy;
  
  unsigned int index_strategies[2] = {0,0};
  int length_plan = 0;

  char line[10000];

  fseek (file, 0, SEEK_SET);
  while (fgets (line, sizeof(line), file) != NULL)
    {
      //process each line of 'file'
      int err = process_line (strategies, index_strategies,
			      &length_plan, rfbb, afbb,
			      rlpb, alpb, line, verbose);
      ASSERT (index_strategies[0] <= rmfb && index_strategies[1] <= amfb);
      if (err == -1)
	return NULL;
    }

  /*
    For each strategy, one finds what is the last method used on each
    side.
  */
  for (r = 0; r <= rmfb; r++)
    for (a = 0; a <= amfb; a++){
      facul_both_strategy_t* st = &strategies->strategy[r][a];
      ASSERT (st != NULL);
      facul_method_t* methods = st->methods;
      int index_last_method[2] = {-1, -1};
      ASSERT (methods != NULL);
      unsigned int i;
      for (i = 0; methods[i].method != 0; i++) {
	methods[i].is_the_last = 0;
	index_last_method[methods[i].side] = i;
      }
      /*if index_last_method[*] that imply that it doesn't exist method
	for one side!!!*/
      ASSERT (index_last_method[0] != -1 &&
	      index_last_method[1] != -1);
      methods[index_last_method[0]].is_the_last = 1;
      methods[index_last_method[1]].is_the_last = 1;
    }
  return strategies;
}


void 
facul_clear_strategies (facul_strategies_t *strategies)
{
  if (strategies == NULL)
    return ;
  
  //free plan
  precompute_plan_t* plan = strategies->plan;
  for (int j = 0; plan[j].method!=0; j++)
    {
      if (plan[j].method == PM1_METHOD)
	pm1_clear_plan (plan[j].plan);
      else if (plan[j].method == PP1_27_METHOD ||
	       plan[j].method == PP1_65_METHOD)
	pp1_clear_plan (plan[j].plan);
      else if (plan[j].method == EC_METHOD)
	ecm_clear_plan (plan[j].plan);
      free (plan[j].plan);
      plan[j].method = 0;
      plan[j].plan = NULL;
    } 
  free (plan);
  plan = NULL;

  //free strategies
  unsigned int r;
  for (r = 0; r <= strategies->mfb[0]; r++) {
    unsigned int a;
    for (a = 0; a <= strategies->mfb[1]; a++)
      free (strategies->strategy[r][a].methods);
    free (strategies->strategy[r]);
  }
  free (strategies->strategy);
  free (strategies);
  
}

int
facul_fprint_strategies (FILE* file, facul_strategies_t* strategies)
{
  if (file == NULL)
    return -1;
  unsigned int r;
  for (r = 0; r <= strategies->mfb[0]; r++) {
    unsigned int a;
    for (a = 0; a <= strategies->mfb[1]; a++) {
      printf ("[r = %d, a = %d]", r, a);
      facul_both_strategy_t* st = &strategies->strategy[r][a];
      //print info lpb ...
      if (st == NULL)
	continue;
      printf ("(lpb = [%ld,%ld], as...=[%lf, %lf], BBB = [%lf, %lf])\n",
	      st->lpb[0], st->lpb[1], st->assume_prime_thresh[0],
	      st->assume_prime_thresh[1], st->BBB[0], st->BBB[1]);
      facul_method_t* methods = st->methods;
      if (methods == NULL)
      	continue;
      for (int i = 0; methods[i].method != 0; i++)
	{
	  if (methods[i].plan == NULL)//zero method!!!
	    printf ("[side=%d, FM=%ld, B1=0, B2=0] ", methods[i].side,
		    methods[i].method);
	  else {
	    if (methods[i].method == PM1_METHOD)
	      {
		pm1_plan_t* plan = (pm1_plan_t*) methods[i].plan;
		printf ("[side=%d, FM=%ld, B1=%d, B2=%d] ", methods[i].side,
			methods[i].method, plan->B1, plan->stage2.B2);
	      }
	    else if (methods[i].method == PP1_27_METHOD ||
		     methods[i].method == PP1_65_METHOD)
	      {
		pp1_plan_t* plan = (pp1_plan_t*) methods[i].plan;
		printf ("[side=%d, FM=%ld, B1=%d, B2=%d] ", methods[i].side,
			methods[i].method, plan->B1, plan->stage2.B2);
	      }
	    else if (methods[i].method == EC_METHOD)
	      {
		ecm_plan_t* plan = (ecm_plan_t*) methods[i].plan;
		printf ("[side=%d, FM=%ld, B1=%d, B2=%d] ", methods[i].side,
			methods[i].method, plan->B1,plan->stage2.B2);
	      }
	    else
	      return -1;
	  }
	}
      printf ("\n");
    }
  }
  return 1;
}


/********************************************/
/*            modset_t                      */
/********************************************/


void 
modset_clear (modset_t *modset)
{
  //ASSERT_ALWAYS(modset->arith != CHOOSE_NONE);//why!!!
  switch (modset->arith) {
  case CHOOSE_NONE: /* already clear */
    break;
  case CHOOSE_UL:
    modredcul_clearmod (modset->m_ul);
    break;
//#if     MOD_MAXBITS > MODREDCUL_MAXBITS
  case CHOOSE_15UL:
    modredc15ul_clearmod (modset->m_15ul);
    break;
//#endif
//#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
  case CHOOSE_2UL2:
    modredc2ul2_clearmod (modset->m_2ul2);
    break;
//#endif
//#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
  case CHOOSE_MPZ:
    modmpz_clearmod (modset->m_mpz);
    break;
//#endif
  default:
    abort();
  }
  modset->arith = CHOOSE_NONE;
}


/*
  This function tries to factor a pair of cofactors (m[0], m[1]) from
  strategies. It returns the number of factors found on each side, or
  -1 if the factor is not smooth.
  Remarks: - the values of factors found are stored in 'factors'.
           - the variable 'is_smooth' allows to know if a cofactor
             is already factored.
 */

static int*
facul_both_src (unsigned long **factors, const modset_t* m,
		 const facul_strategies_t *strategies, int* cof,
		 int* is_smooth)
{
  int* found = calloc(2, sizeof(int));
  
  facul_both_strategy_t* strategy = &strategies->strategy[cof[0]][cof[1]];
  
  if (strategy->methods == NULL)
    return found;

  for (int i = 0; strategy->methods[i].method != 0; i++)
    {
      modset_t fm, cfm;
      unsigned int len_fm = 0;
      unsigned int len_cfm = 0;
      int side = strategy->methods[i].side;
      if (is_smooth[side])
	continue; 

      int res_fac = 0;
      switch (m[side].arith) {
      case CHOOSE_UL:
	res_fac = facul_doit_onefm_ul(factors[side], m[side].m_ul,
				      strategy->methods[i], &fm, &cfm,
				      &len_fm, &len_cfm, strategy->lpb[side],
				      strategy->assume_prime_thresh[side],
				      strategy->BBB[side]);
	break;
//#if     MOD_MAXBITS > MODREDCUL_MAXBITS
      case CHOOSE_15UL:
	res_fac = facul_doit_onefm_15ul(factors[side], m[side].m_15ul,
					strategy->methods[i], &fm, &cfm,
					&len_fm, &len_cfm, strategy->lpb[side],
					strategy->assume_prime_thresh[side],
					strategy->BBB[side]);
	break;
//#endif
//#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
      case CHOOSE_2UL2:
	res_fac = facul_doit_onefm_2ul2 (factors[side], m[side].m_2ul2,
					 strategy->methods[i], &fm, &cfm,
					 &len_fm, &len_cfm, strategy->lpb[side],
					 strategy->assume_prime_thresh[side],
					 strategy->BBB[side]);
	break;
//#endif
//#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
      case CHOOSE_MPZ:
	res_fac = facul_doit_onefm_mpz (factors[side], m[side].m_mpz,
					strategy->methods[i], &fm, &cfm,
					&len_fm, &len_cfm, strategy->lpb[side],
					strategy->assume_prime_thresh[side],
					strategy->BBB[side]);
	break;
//#endif
      default: abort();
      }
      //check our result!
      //res_fac contains the number of factors found!

      if (res_fac == -1)
	{
	  /*
	    The cofactor m[side] is not smooth. So, one stops the
	    cofactorisation.
	  */
	  found[side] = -1;
	  break;
	}
      if (res_fac == 0)
	{
	  /* Zero factor found. If it was the last method for this
	     side, then one stops the cofactorisation. Otherwise, one
	     tries with an other method! */
	  if (strategy->methods[i].is_the_last)
	      break;
	  else
	    continue;
	}
      found[side] += res_fac;

      if (res_fac == 2)
	/*
	  Indeed, using only one factoring method, you found two
	  primes factors of m (f, m/f) then m is factored and your
	  work is finished for this cofactor!!
	*/
	is_smooth[side] = true;
      
      /*
	res_fac == 1!  Only one factor has been found. Hence, our
	factorisation is not finished.
       */
      if (fm.arith != CHOOSE_NONE)
	{
	  unsigned long** factors2 = malloc(2 * sizeof(factors2));
	  factors2[side] = factors[side]+res_fac;
	  factors2[(1+side)%2] = factors[(1+side)%2];
	  modset_t m2[2];
	  m2[side] = fm;
	  m2[(1+side)%2] = m[(1+side)%2];
	  cof[side] =  len_fm;
	  int* found2 = facul_both_src (factors2, m2, strategies, cof,
					is_smooth);
	  free (factors2);

	  found[0] += found2[0];
	  found[1] += found2[1];
	  free (found2);
	  modset_clear (&fm);
	  break;
	}
      if (cfm.arith != CHOOSE_NONE)
	{
	  unsigned long** factors2 = malloc(2 * sizeof(factors2));
	  factors2[side] = factors[side]+res_fac;
	  factors2[(1+side)%2] = factors[(1+side)%2];
	  modset_t m2[2];
	  m2[side] = cfm;
	  m2[(1+side)%2] = m[(1+side)%2];
	  cof[side] =  len_cfm;
	  int* found2 = facul_both_src (factors2, m2, strategies, cof,
					is_smooth);

	  free (factors2);
	  found[0] += found2[0];
	  found[1] += found2[1];
	  free (found2);
	  //free modset
	  modset_clear (&cfm);
	  break;
	}
    }
  
  return found;
}


/*
  This function is like facul, but we will work with the both norms
  together.
  returns the number of factors for each side!
*/
int*
facul_both (unsigned long **factors, mpz_t* N,
	    const facul_strategies_t *strategies, int* is_smooth)
{
  int cof[2];
  size_t bits;
  int* found = NULL;

  modset_t n[2];
  n[0].arith = CHOOSE_NONE;
  n[1].arith = CHOOSE_NONE;
  
#ifdef PARI
  gmp_fprintf (stderr, "(%Zd %Zd)", N[0], N[1]);
#endif

  if (mpz_sgn (N[0]) <= 0 || mpz_sgn (N[1]) <= 0)
      return found;
  
  if (mpz_cmp_ui (N[0], 1UL) == 0)
    is_smooth[0] = 1;
  if (mpz_cmp_ui (N[1], 1UL) == 0)
    is_smooth[1] = 1;

  for (int side = 0; side < 2; side++)
    {
      /* If the composite does not fit into our modular arithmetic, return
	 no factor */
      bits = mpz_sizeinbase (N[side], 2);
      cof[side] = bits;
      if (bits > MODMPZ_MAXBITS)
	return 0;
  
      /* Use the fastest modular arithmetic that's large enough for this input */
      if (bits <= MODREDCUL_MAXBITS)
	{
	  ASSERT(mpz_fits_ulong_p(N[side]));
	  modredcul_initmod_ul (n[side].m_ul, mpz_get_ui(N[side]));
	  n[side].arith = CHOOSE_UL;
	}
      else if (bits <= MODREDC15UL_MAXBITS)
	{
	  unsigned long t[2];
	  modintredc15ul_t m;
	  size_t written;
	  mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, N[side]);
	  ASSERT_ALWAYS(written <= 2);
	  modredc15ul_intset_uls (m, t, written);
	  modredc15ul_initmod_int (n[side].m_15ul, m);
	  n[side].arith = CHOOSE_15UL;
	}
      else if (bits <= MODREDC2UL2_MAXBITS)
	{
	  unsigned long t[2];
	  modintredc2ul2_t m;
	  size_t written;
	  mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, N[side]);
	  ASSERT_ALWAYS(written <= 2);
	  modredc2ul2_intset_uls (m, t, written);
	  modredc2ul2_initmod_int (n[side].m_2ul2, m);
	  n[side].arith = CHOOSE_2UL2;
	}
      else
	{
	  modmpz_initmod_int (n[side].m_mpz, N[side]);
	  n[side].arith = CHOOSE_MPZ;
	}
      ASSERT (n[side].arith != CHOOSE_NONE);
    }

  found = facul_both_src (factors, n, strategies, cof, is_smooth);

  for (int side = 0; side < 2; side++)
    {
      if (found[side] > 1)
	{
	  /* Sort the factors we found */
	  qsort (factors[side], found[side], sizeof (unsigned long),
		 (int (*)(const void *, const void *)) &cmp_ul);
	}

#ifdef PARI
      if (found[side] > 1)
	{
	  fprintf (stderr, " == ");
	  for (i = 0; i < found[side]; i++)
	    fprintf (stderr, "%lu%s", factors[side][i],
		     (i+1 < found[side]) ? " * " : " /* PARI */\n");
	}
      else
	fprintf (stderr, "; /* PARI */\n");
#endif
    }

  //Free
  modset_clear (&n[0]);
  modset_clear (&n[1]);
  
  return found;
}


