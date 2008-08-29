#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <primegen.h>
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "modredc_ul.h"
#include "modredc_ul_default.h"


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
  unsigned long start, stop, x0 = 2UL, i, B1, B2, p, cof = 1UL;
  unsigned long hits = 0, hits_input = 0, total = 0;
  mpz_t E;
  int compare = 0;
  int method = 0; /* 0 = P-1, 1 = P+1, 2 = ECM */
  int parameterization = 0;
  int only_primes = 0, verbose = 0;
  int printfactors = 0;
  primegen pg[1];
  pm1_plan_t pm1plan;
  pp1_plan_t pp1plan;
  ecm_plan_t ecmplan;

  /* Parse options */
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
	  cof = strtoul (argv[2], NULL, 10);
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


  if (method == 0)
    pm1_make_plan (&pm1plan, B1, B2, (verbose >= 2));
  if (method == 1)
    pp1_make_plan (&pp1plan, B1, B2, (verbose >= 2));
  if (method == 2)
    ecm_make_plan (&ecmplan, B1, B2, BRENT12, x0, (verbose >= 2));

  /* Compute exponent (a.k.a. multiplier for P+1/ECM) for stage 1.
     This is currently only used for comparing speed/results of
     P+1 byte code and binary chain. */
  mpz_init (E);
  mpz_set_ui (E, 1UL);
  primegen_init (pg);
  for (p = primegen_next (pg); p <= B1; p = primegen_next (pg))
    {
      unsigned long q;
      for (q = p; q * p < B1; q *= p);
      mpz_mul_ui (E, E, q);
    }
  if (verbose)
    gmp_printf ("E = %Zd;\n", E);


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
      modulus_t m;
      residue_t b, r, r2;
      unsigned long f, ic = i * cof;

      total++;

      /* Init the modulus */
      mod_initmod_ul (m, ic);
      mod_init_noset0 (r, m);
      mod_init_noset0 (b, m);

      if (method == 0) /* P-1 */
	{
	  f = pm1_15 (r, m, &pm1plan);

	  if (verbose >= 2)
	    printf ("Mod(%lu, %lu)^E == %lu\n", x0, ic, mod_get_ul (r, m));
	}
      else if (method == 1) /* P+1 */
        {
	  if (1) /* Switch between byte code interpreter and mod_V_mp() */
	    f = pp1 (r, m, &pp1plan);
	  else
	    {
	      mod_set_ul (b, 7UL, m);
	      mod_inv (b, b, m);
	      mod_add (b, b, b, m); /* b = 2/7 (mod m) */
	      mod_V_mp (r, b, E->_mp_d, E->_mp_size, m);
	    }
	  
          if (compare)
            {
	      mod_init_noset0 (r2, m);
	      mod_set_ul (b, 7UL, m);
	      mod_inv (b, b, m);
	      mod_add (b, b, b, m); /* b = 2/7 (mod m) */
              mod_V_mp (r2, b, E->_mp_d, E->_mp_size, m);
              if (!mod_equal (r, r2, m))
                printf ("Error, pp1() and mod_V_mp() differ for "
                        "modulus %lu\n", i);
	      else if (verbose >= 2)
		printf ("pp1() and mod_V_mp() agree for modulus %lu\n", i);

	      mod_clear (r2, m);
            }

	  if (verbose >= 2)
	    printf ("V(E, Mod(%lu, %lu) == %lu\n", x0, ic, mod_get_ul (r, m));
	}
      else
	{
	  f = ecm (r, m, &ecmplan);
	  if (verbose >= 2)
	    printf ("x-coordinate after stage 1 is Mod(%lu, %lu)\n",
		    mod_get_ul (r, m), ic);
	}
      
      if (printfactors && f > 1UL)
	printf ("%lu\n", f);

      if (f > 1UL && (cof == 1 || f != ic))
	hits++;
      if (cof != 1 && f == ic)
	hits_input++;

      mod_clear (r, m);
      mod_clear (b, m);
      mod_clearmod (m);

      if (only_primes)
	i = primegen_next (pg);
      else
	i += 2;
    }

  if (method == 0)
    pm1_clear_plan (&pm1plan);
  if (method == 1)
    pp1_clear_plan (&pp1plan);
  if (method == 2)
    ecm_clear_plan (&ecmplan);

  mpz_clear (E);
  printf ("Of %lu %s in [%lu, %lu] there were %lu with smooth order\n", 
	  total, (only_primes) ? "primes" : "odd numbers", start, stop, hits);
  printf ("Ratio: %f\n", (double)hits / (double)total);
  if (cof > 1)
    {
      printf ("Input number was found %lu times\n", hits_input);
    }

  return 0;
}
