#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <gmp.h>
#include "facul.h"
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "modredc_ul.h"
#include "modredc_ul_default.h"
#include "getprime.h"

#define MAX_METHODS 20

const char *method_name[] = {"P-1", "P+1", "ECM"};

void print_help (char *programname)
{
  printf ("%s [options] <start> <stop> <B1> <B2>\n", programname);
  printf ("Run factoring method on numbers in interval [<start>, <stop>]\n");
  printf ("-pm1 <B1> <B2>      Run P-1 with B1, B2\n");
  printf ("-pp1 <B1> <B2>      Run P+1 with B1, B2\n");
  printf ("-ecm <B1> <B2> <s>  Run ECM with B1, B2 and sigma s\n");
  printf ("-strat   Use the facul default strategy. Don't use with -pm1, -pp1, -ecm\n");
  printf ("-fbb <n> Use <n> as factor base bound, e.g. for primality checks\n");
  printf ("-lpb <n> Use <n> as large prime bound, e.g. for early abort\n");
  printf ("-p       Try only primes in [<start>, <stop>] (default: all odd "
	  "numbers)\n");
  printf ("-v       More verbose output. Once: print parameters. Twice: print "
	  "residues\n");
  printf ("-vf      Print factors that are found\n");
  printf ("-cof <n> Multiply each number to test by <num> before calling facul.\n"
	  "         NOTE: facul does not report input numbers as factors,\n"
	  "         with -p you MUST use -cof or no factors will be found\n");
  printf ("-m <n>   Choose modulus, do separate counts of p %% modulus\n");
  printf ("-inp <f> Read decimal numbers to factor from file <f>\n");
  printf ("-inpraw <f>  Read numbers in GMP raw format from file <f>\n");
  printf ("-inpstop <n> Stop after reading <n> numbers\n");
}

int main (int argc, char **argv)
{
  unsigned long start, stop, i, mod = 0UL, inpstop = ULONG_MAX;
  unsigned long hits = 0, total = 0;
  unsigned long fbb = 0, lpb = ~(0UL);
  char *inp_fn = NULL;
  FILE *inp;
  struct rusage usage;
  mpz_t N, cof;
  facul_strategy_t *strategy;
  int nr_methods = 0;
  int only_primes = 0, verbose = 0;
  int printfactors = 0;
  int inp_raw = 0;
  int got_usage;
  int strat = 0;
  unsigned long *primmod = NULL, *hitsmod = NULL;
  unsigned long f[16];
  int facul_code;

  strategy = malloc (sizeof(facul_strategy_t));
  strategy->methods = malloc ((MAX_METHODS + 1) * sizeof (facul_method_t));
  strategy->lpb = ~(0UL);
  strategy->fbb2[0] = 0UL;
  strategy->fbb2[1] = 0UL;

  /* Parse options */
  mpz_init (N);
  mpz_init (cof);
  mpz_set_ui (cof, 1UL);
  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 1 && strcmp (argv[1], "-h") == 0)
	{
	  print_help (argv[0]);
	  return 0;
	}
      else if (argc > 3 && strcmp (argv[1], "-pm1") == 0 && 
	       nr_methods < MAX_METHODS)
	{
	  unsigned long B1, B2;
	  B1 = strtoul (argv[2], NULL, 10);
	  B2 = strtoul (argv[3], NULL, 10);
	  strategy->methods[nr_methods].method = PM1_METHOD;
	  strategy->methods[nr_methods].plan = malloc (sizeof (pm1_plan_t));
	  ASSERT (strategy->methods[nr_methods].plan != NULL);
	  pm1_make_plan (strategy->methods[nr_methods].plan, B1, B2, 
			 (verbose >= 2));
	  nr_methods++;
	  argc -= 3;
	  argv += 3;
	}
      else if (argc > 3 && strcmp (argv[1], "-pp1") == 0 && 
	       nr_methods < MAX_METHODS)
	{
	  unsigned long B1, B2;
	  B1 = strtoul (argv[2], NULL, 10);
	  B2 = strtoul (argv[3], NULL, 10);
	  strategy->methods[nr_methods].method = PP1_METHOD;
	  strategy->methods[nr_methods].plan = malloc (sizeof (pp1_plan_t));
	  ASSERT (strategy->methods[nr_methods].plan != NULL);
	  pp1_make_plan (strategy->methods[nr_methods].plan, B1, B2, 
			 (verbose >= 2));
	  nr_methods++;
	  argc -= 3;
	  argv += 3;
	}
      else if (argc > 4 && strcmp (argv[1], "-ecm") == 0 && 
	       nr_methods < MAX_METHODS)
	{
	  unsigned long B1, B2;
	  long sigma;
	  B1 = strtoul (argv[2], NULL, 10);
	  B2 = strtoul (argv[3], NULL, 10);
	  sigma = strtol (argv[4], NULL, 10);
	  strategy->methods[nr_methods].method = EC_METHOD;
	  strategy->methods[nr_methods].plan = malloc (sizeof (ecm_plan_t));
	  ASSERT (strategy->methods[nr_methods].plan != NULL);
	  ecm_make_plan (strategy->methods[nr_methods].plan, B1, B2, 
                         (sigma > 0) ? BRENT12 : MONTY12, labs (sigma), 
                         (verbose >= 2));
	  nr_methods++;
	  argc -= 4;
	  argv += 4;
	}
      else if (argc > 1 && strcmp (argv[1], "-strat") == 0 && 
	       nr_methods == 0)
        {
	  strat = 1;
	  argc -= 1;
	  argv += 1;
	}
      else if (argc > 2 && strcmp (argv[1], "-fbb") == 0)
	{
	  fbb = strtoul (argv[2], NULL, 10);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-lpb") == 0)
	{
	  lpb = strtoul (argv[2], NULL, 10);
	  argc -= 2;
	  argv += 2;
	}
     else if (argc > 1 && strcmp (argv[1], "-p") == 0)
	{
	  only_primes = 1;
	  argc -= 1;
	  argv += 1;
	}
      else if (argc > 1 && strcmp (argv[1], "-v") == 0)
	{
	  verbose++;
	  argc -= 1;
	  argv += 1;
	}
      else if (argc > 1 && strcmp (argv[1], "-vf") == 0)
	{
	  printfactors = 1;
	  argc -= 1;
	  argv += 1;
	}
      else if (argc > 2 && strcmp (argv[1], "-inpstop") == 0)
	{
	  inpstop = strtoul (argv[2], NULL, 10);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-m") == 0)
	{
	  mod= strtoul (argv[2], NULL, 10);
	  hitsmod = (unsigned long *) malloc (mod * sizeof (unsigned long));
	  primmod = (unsigned long *) malloc (mod * sizeof (unsigned long));
	  for (i = 0; i < mod; i++)
	    hitsmod[i] = primmod[i]= 0;
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-cof") == 0)
	{
	  mpz_set_str (cof, argv[2], 10);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-inp") == 0)
	{
	  inp_fn = argv[2];
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-inpraw") == 0)
	{
	  inp_fn = argv[2];
	  inp_raw = 1;
	  argc -= 2;
	  argv += 2;
	}
      else
        {
	  printf ("Unrecoglized option: %s\n", argv[1]);
	  exit (EXIT_FAILURE);
        }
    }
  
  if (strat && nr_methods != 0)
    {
      printf ("Don't use -strat with -pm1, -pp1 or -ecm\n");
      exit (EXIT_FAILURE);
    }

  if (strat)
    {
      facul_clear_strategy (strategy);
      strategy = facul_make_strategy (15, fbb, (lpb == 0) ? 0 : 1UL << lpb);
    }
  else
    {
      printf ("Strategy has %d methods\n", nr_methods);
      strategy->methods[nr_methods].method = 0;
    }

  if (inp_fn == NULL)
    {
      if (argc < 3)
	{
	  print_help (argv[0]);
	  return 1;
	}
      
      /* Get range to test */
      start = strtoul (argv[1], NULL, 10);
      if (start % 2UL == 0UL)
	start++;
      stop = strtoul (argv[2], NULL, 10);
      
      if (only_primes)
	{
	  do {
            i = getprime (1);
          } while (i < start);
	}
      else
	i = start;
      
      /* The main loop */
      while (i <= stop)
	{
	  if (mod > 0)
	    primmod[i % mod]++;
	  
	  total++;
	  mpz_mul_ui (N, cof, i);
          if (verbose >= 2)
            gmp_printf ("Trying to factor %Zd\n", N);
          facul_code = facul (f, N, strategy);
	  
	  if (facul_code > 0)
	    {
	      hits++;
	      if (mod > 0)
		hitsmod[i % mod]++;
	    }
	  
	  if (printfactors && facul_code > 0)
	    {
	      int j;
	      for (j = 0; j < facul_code; j++)
		printf ("%lu ", f[j]);
	      printf ("\n");
	    }
	  
	  if (only_primes)
	    i = getprime (1);
	  else
	    i += 2;
	}
  } else {
    inp = fopen (inp_fn, "r");
    if (inp == NULL)
      {
	printf ("Could not open %s\n", inp_fn);
	exit (EXIT_FAILURE);
      }

    /* Read lines from inp */
    while (!feof(inp) && inpstop-- > 0)
      {
	if (inp_raw)
	  mpz_inp_raw (N, inp);
	else
	  mpz_inp_str (N, inp, 0);
	if (mpz_sgn (N) <= 0)
	  continue;
	mpz_mul (N, N, cof);
	total++;
        if (verbose >= 2)
          gmp_printf ("Trying to factor %Zd\n", N);
	facul_code = facul (f, N, strategy);
	
	if (facul_code > 0)
	  hits++;
	
	if (printfactors && facul_code > 0)
	  {
	    int j;
	    for (j = 0; j < facul_code; j++)
	      printf ("%lu ", f[j]);
	    printf ("\n");
	  }
      }
    fclose (inp);
  }
  

  facul_clear_strategy (strategy);

  got_usage = getrusage(RUSAGE_SELF, &usage);

  printf ("Of %lu tries there were %lu with a factor found\n", 
           total, hits);
  printf ("Ratio: %f\n", (double)hits / (double)total);
  
  if (got_usage == 0)
    {
      double usrtime;
      usrtime = (double) usage.ru_utime.tv_sec * 1000000. +
                (double) usage.ru_utime.tv_usec;
      printf ("Total time: %.2f s, per call: %f usec, per factor: %f usec\n",
              usrtime / 1000000., usrtime / (double) total, 
              usrtime / (double) hits);
    }
    
  if (mod > 0)
    {
      printf ("Distribution of factors found over residue classes mod %lu:\n",
              mod);
      for (i = 0; i < mod; i++)
        if (hitsmod[i] > 0)
          printf ("%lu:%lu (%f)%s", 
                  i, hitsmod[i], (double)hitsmod[i] / (double) primmod[i],
                  (i < mod - 1) ? ", " : "\n");
    }

  mpz_clear (N);
  mpz_clear (cof);

  if (mod > 0)
    {
      free (hitsmod);
      free (primmod);
    }

  if (verbose)
    facul_print_stats (stdout);

  return 0;
}
