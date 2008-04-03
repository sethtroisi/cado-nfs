#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <primegen.h>
#include "mod_ul.h"
#include "pp1.h"
#include "ecm.h"

const char *method_name[] = {"P-1", "P+1", "ECM"};

void print_help (char *programname)
{
  printf ("%s [options] <start> <stop> <B1>\n", programname);
  printf ("Run factoring method on numbers in interval [<start>, <stop>]\n");
  printf ("-pm1     Run P-1 (default)\n");
  printf ("-pp1     Run P+1\n");
  printf ("-ecm     Run ECM (not implemented yet)\n");
  printf ("-f       Run fast P+1 (code generated, unrolled code)\n");
  printf ("-c       Run both fast and slow P+1 and compare results\n");
  printf ("-m12     Use Montgomery parameterization with 12 torsion for ECM\n");
  printf ("-p       Try only primes in [<start>, <stop>] (default: all odd "
	  "numbers)\n");
  printf ("-v       More verbose output. Once: print parameters. Twice: print "
	  "residues\n");
  printf ("-x0 <n>  Use <n> as starting value for P-1/P+1, or as sigma "
           "value for ECM\n");
}

int main (int argc, char **argv)
{
  unsigned long start = 999, stop = 1001, x0 = 7UL, i, B1, p;
  unsigned long hits = 0, total = 0;
  mpz_t E;
  int fast = 0, compare = 0;
  int method = 0; /* 0 = P-1, 1 = P+1, 2 = ECM */
  int parameterization = 0;
  int only_primes = 0, verbose = 0;
  primegen pg[1];

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
	  argc--;
	  argv++;
	}
      else if (argc > 1 && strcmp (argv[1], "-ecm") == 0)
	{
	  method = 2;
	  argc--;
	  argv++;
	}
      else if (argc > 1 && strcmp (argv[1], "-f") == 0)
	{
	  fast = 1;
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
      else if (argc > 2 && strcmp (argv[1], "-x0") == 0)
	{
	  x0 = strtoul (argv[2], NULL, 10);;
	  argc -= 2;
	  argv += 2;
	}
    }
  
  if (argc < 4)
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
  if (verbose)
    {
      printf ("Running %s on %s in [%lu, %lu] with B1 = %lu, %s=%lu\n",
              method_name[method], 
              (only_primes) ? "primes" : "odd numbers", 
              start, stop, B1, 
              (method == 2) ? "sigma" : "x0", x0);
    }


  /* Compute exponent (a.k.a. multiplier for P+1/ECM) for stage 1 */
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
      unsigned long invm;

      total++;

      /* Init the modulus */
      mod_initmod_ul (m, i);
      mod_init_noset0 (r, m);
      mod_set_ul_reduced (b, x0, m);
      invm = -mod_invmodlong (m);
      mod_init_noset0 (b, m);
      mod_set_ul_reduced (b, x0, m);
      mod_tomontgomery (b, b, m);

      if (method == 0) /* P-1 */
	{
	  mod_powredc_mp (r, b, E->_mp_d, E->_mp_size, invm, m);
	  mod_frommontgomery (r, r, invm, m);
	  if (verbose >= 2)
	    printf ("Mod(%lu, %lu)^E == %lu\n", x0, i, mod_get_ul (r, m));
	  if (mod_get_ul (r, m) == 1UL)
	    hits++;
	}
      else if (method == 1) /* P+1 */
        {
          if (fast)
            pp1_stage1 (r, r, b, (int) B1, invm, m);
          else if (compare)
            {
	      mod_init_noset0 (r2, m);
              pp1_stage1 (r, r, b, (int) B1, invm, m);
              mod_Vredc_mp (r2, b, E->_mp_d, E->_mp_size, invm, m);
              if (!mod_equal (r, r2, m))
                printf ("Error, pp1_stage1_100() and mod_Vredc_mp() differ for "
                        "modulus %lu\n", i);
	      mod_clear (r2, m);
            }
          else
            mod_Vredc_mp (r, b, E->_mp_d, E->_mp_size, invm, m);
          mod_frommontgomery (r, r, invm, m);

	  if (verbose >= 2)
	    printf ("V(E, Mod(%lu, %lu) == %lu\n", x0, i, mod_get_ul (r, m));

          if (mod_get_ul (r, m) == 2UL)
            hits++;
	}
      else
	{
	  mod_set_ul_reduced (b, x0, m);
	  if (ecm_stage1 (r, (int) B1, b, parameterization, m))
	    {
	      if (verbose >= 2)
		printf ("Found %lu\n", mod_gcd (r, m));
	      hits++;
	    }
	  else if (verbose >= 2)
	    printf ("Sigma = %lu, x1 = %lu\n", 
		    mod_get_ul (b, m), mod_get_ul (r, m));
	}
      
      mod_clear (r, m);
      mod_clear (b, m);
      mod_clearmod (m);

      if (only_primes)
	i = primegen_next (pg);
      else
	i += 2;
    }

  mpz_clear (E);
  printf ("Of %lu %s in [%lu, %lu] there were %lu with smooth order\n", 
	  total, (only_primes) ? "primes" : "odd numbers", start, stop, hits);
  printf ("Ratio: %f\n", (double)hits / (double)total);
  return 0;
}
