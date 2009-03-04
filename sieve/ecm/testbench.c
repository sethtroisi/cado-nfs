#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <gmp.h>
#include <primegen.h>
#include "facul.h"
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "modredc_ul.h"
#include "modredc_ul_default.h"

#define MAX_METHODS 20

const char *method_name[] = {"P-1", "P+1", "ECM"};

void print_help (char *programname)
{
  printf ("%s [options] <start> <stop> <B1> <B2>\n", programname);
  printf ("Run factoring method on numbers in interval [<start>, <stop>]\n");
  printf ("-pm1 <B1> <B2>      Run P-1 with B1, B2\n");
  printf ("-pp1 <B1> <B2>      Run P+1 with B1, B2\n");
  printf ("-ecm <B1> <B2> <s>  Run ECM with B1, B2 and sigma s\n");
  printf ("-m12     Use Montgomery parameterization with 12 torsion for ECM\n");
  printf ("-p       Try only primes in [<start>, <stop>] (default: all odd "
	  "numbers)\n");
  printf ("-v       More verbose output. Once: print parameters. Twice: print "
	  "residues\n");
  printf ("-vf      Print factors that are found\n");
  printf ("-cof <n> Multiply each number to test by <num> before calling\n"
	  "         factoring routine\n");
  printf ("-m <n>   Choose modulus, do separate counts of p %% modulus\n");
}

int main (int argc, char **argv)
{
  unsigned long start, stop, i, mod = 0UL;
  unsigned long hits = 0, hits_input = 0, total = 0;
  char *inp_fn = NULL;
  FILE *inp;
  struct rusage usage;
  mpz_t N, cof;
  facul_strategy_t *strategy;
  primegen pg[1];
  int nr_methods = 0;
  int parameterization = 0;
  int only_primes = 0, verbose = 0;
  int printfactors = 0;
  int got_usage;
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
	  pp1_make_plan (strategy->methods[nr_methods].plan, B1, B2, 
			 (verbose >= 2));
	  nr_methods++;
	  argc -= 3;
	  argv += 3;
	}
      else if (argc > 4 && strcmp (argv[1], "-ecm") == 0 && 
	       nr_methods < MAX_METHODS)
	{
	  unsigned long B1, B2, sigma;
	  B1 = strtoul (argv[2], NULL, 10);
	  B2 = strtoul (argv[3], NULL, 10);
	  sigma = strtoul (argv[4], NULL, 10);
	  strategy->methods[nr_methods].method = EC_METHOD;
	  strategy->methods[nr_methods].plan = malloc (sizeof (ecm_plan_t));
	  ecm_make_plan (strategy->methods[nr_methods].plan, B1, B2, BRENT12, 
			 sigma, (verbose >= 2));
	  nr_methods++;
	  argc -= 4;
	  argv += 4;
	}
      else if (argc > 1 && strcmp (argv[1], "-m12") == 0)
	{
	  parameterization = 1;
	  argc--;
	  argv++;
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
      else
        {
	  printf ("Unrecoglized option: %s\n", argv[1]);
	  exit (EXIT_FAILURE);
        }
    }
  
  printf ("Strategy has %d methods\n", nr_methods);
  strategy->methods[nr_methods].method = 0;

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
	  primegen_init (pg);
	  primegen_skipto (pg, start);
	  i = primegen_next (pg);
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
	      i = primegen_next (pg);
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
      while (!feof(inp))
      {
	  mpz_inp_str (N, inp, 0);
	  if (mpz_sgn (N) <= 0)
	      continue;
	  mpz_mul (N, N, cof);
	  total++;
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
  if (mpz_cmp_ui (cof, 1UL) != 0)
    {
      printf ("Input number was found %lu times\n", hits_input);
    }
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
