#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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
  printf ("-pm1     Run P-1 (default)\n");
  printf ("-pp1     Run P+1\n");
  printf ("-ecm     Run ECM\n");
  printf ("-c       Run both bytecode and binary chain P+1 and compare results\n");
  printf ("-m12     Use Montgomery parameterization with 12 torsion for ECM\n");
  printf ("-p       Try only primes in [<start>, <stop>] (default: all odd "
	  "numbers)\n");
  printf ("-v       More verbose output. Once: print parameters. Twice: print "
	  "residues\n");
  printf ("-vf      Print factors that are found\n");
  printf ("-x0 <n>  Use <n> as starting value for P-1/P+1, or as sigma "
           "value for ECM\n");
  printf ("-cof <n> Multiply numbers from search interval by <num> before "
	  "calling\n"
	  "         factoring routine\n");
}

int main (int argc, char **argv)
{
  unsigned long start, stop, x0 = 2UL, i, B1, B2;
  unsigned long hits = 0, hits_input = 0, total = 0;
  mpz_t N, cof;
  facul_strategy_t *strategy;
  primegen pg[1];
  int nr_methods = 0;
  int compare = 0;
  int method = 0; /* 0 = P-1, 1 = P+1, 2 = ECM */
  int parameterization = 0;
  int only_primes = 0, verbose = 0;
  int printfactors = 0;

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
      else if (argc > 1 && strcmp (argv[1], "-pm1") == 0)
	{
	  method = 0;
	  argc--;
	  argv++;
	}
      else if (argc > 1 && strcmp (argv[1], "-pp1") == 0)
	{
	  method = 1;
	  x0 = 7;
	  argc--;
	  argv++;
	}
      else if (argc > 1 && strcmp (argv[1], "-ecm") == 0)
	{
	  method = 2;
	  x0 = 7;
	  argc--;
	  argv++;
	}
      else if (argc > 1 && strcmp (argv[1], "-m12") == 0)
	{
	  parameterization = 1;
	  argc--;
	  argv++;
	}
      else if (argc > 1 && strcmp (argv[1], "-c") == 0)
	{
	  compare = 1;
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
      else if (argc > 2 && strcmp (argv[1], "-x0") == 0)
	{
	  x0 = strtoul (argv[2], NULL, 10);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-cof") == 0)
	{
	  mpz_set_str (cof, argv[2], 10);
	  argc -= 2;
	  argv += 2;
	}
    }
  
  if (argc < 5)
    {
      print_help (argv[0]);
      return 1;
    }

  /* Get parameters */
  start = strtoul (argv[1], NULL, 10);
  if (start % 2UL == 0UL)
    start++;
  if (start <= x0)
    start = (x0 + 1UL) | 1UL;
  stop = strtoul (argv[2], NULL, 10);
  B1 = strtoul (argv[3], NULL, 10);
  B2 = strtoul (argv[4], NULL, 10);

  if (verbose)
    {
      printf ("Running %s on %s in [%lu, %lu] with B1 = %lu, B2 = %lu, "
	      "%s=%lu\n",
              method_name[method], 
              (only_primes) ? "primes" : "odd numbers", 
              start, stop, B1, B2, 
              (method == 2) ? "sigma" : "x0", x0);
    }


  strategy = malloc (sizeof(facul_strategy_t));
  strategy->methods = malloc (MAX_METHODS * sizeof (facul_method_t));
  strategy->lpb = ~(0UL);
  strategy->fbb2[0] = 0UL;
  strategy->fbb2[1] = 0UL;
  if (method == 0)
    {
      strategy->methods[nr_methods].method = PM1_METHOD;
      strategy->methods[nr_methods].plan = malloc (sizeof (pm1_plan_t));
      pm1_make_plan (strategy->methods[nr_methods].plan, B1, B2, (verbose >= 2));
      nr_methods++;
    }
  else if (method == 1)
    {
      strategy->methods[nr_methods].method = PP1_METHOD;
      strategy->methods[nr_methods].plan = malloc (sizeof (pp1_plan_t));
      pp1_make_plan (strategy->methods[nr_methods].plan, B1, B2, (verbose >= 2));
      nr_methods++;
    }
  else if (method == 2)
    {
      strategy->methods[nr_methods].method = EC_METHOD;
      strategy->methods[nr_methods].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (strategy->methods[nr_methods].plan, B1, B2, BRENT12, x0, (verbose >= 2));
      nr_methods++;
    }
  strategy->methods[nr_methods].method = 0;


  /* The main loop */
  if (only_primes)
    {
      primegen_init (pg);
      primegen_skipto (pg, start);
      i = primegen_next (pg);
    }
  else
    i = start;
    
  while (i <= stop)
    {
      unsigned long f[16];
      int facul_code;

      total++;
      mpz_mul_ui (N, cof, i);
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
      
      if (only_primes)
	i = primegen_next (pg);
      else
	i += 2;
    }

  facul_clear_strategy (strategy);

  printf ("Of %lu %s in [%lu, %lu] there were %lu with smooth order\n", 
	  total, (only_primes) ? "primes" : "odd numbers", start, stop, hits);
  printf ("Ratio: %f\n", (double)hits / (double)total);
  if (mpz_cmp_ui (cof, 1UL) != 0)
    {
      printf ("Input number was found %lu times\n", hits_input);
    }
  mpz_clear (N);
  mpz_clear (cof);

  facul_print_stats (stdout);

  return 0;
}
