/* Usage: find_primes_hit P B1 sigma - find all primes <= P that are found with
   ECM(B1,sigma) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "getprime.c"

#undef MODTRACE
#include "mod_ul.h"

#define BRENT12 0
#define MONTY12 1

int verbose = 0;


static void
usage (const char *s)
{
  printf ("Usage: %s Pmin Pmax B1 sigma\n", s); 
  printf ("Find all primes >= Pmin, <= Pmax hit by ECM(B1,sigma)\n");
  printf ("Options:\n");
  printf ("-v        Increase verbosity\n");
  printf ("-gnuplot  Make output suitable for gnuplot\n");
  printf ("-c        Print average exponent of 2, 3, 5, 7 in group order\n");
  printf ("-m12      Make curves of rational 12 torsion from Montgomery's thesis\n");
  exit (1);
}

unsigned int
prime_exponent (unsigned long n, unsigned long p)
{
  unsigned int i = 0;
  while (n % p == 0)
  {
    i++;
    n /= p;
  }
  return i;
}

int
main (int argc, char *argv[])
{
  unsigned long Pmin, Pmax, B1, sigma;
  unsigned long primes = 0, found = 0;
  double p;
  unsigned long nr_points, pow2 = 0, pow3 = 0, pow5 = 0, pow7 = 0;
  int do_count = 0, gnuplot = 0;
  int parameterization = BRENT12;
  const char *programname = argv[0];

  while (argc > 1 && argv[1][0] == '-')
    {
	if (strcmp(argv[1], "-v") == 0)
	{
	    verbose++;
	    argv++;
	    argc--;
	}
	else if (strcmp(argv[1], "-gnuplot") == 0)
	{
	    gnuplot = 1;
	    argv++;
	    argc--;
	}
	else if (strcmp(argv[1], "-c") == 0)
	{
	    do_count = 1;
	    argv++;
	    argc--;
	}
	else if (strcmp(argv[1], "-m12") == 0)
	{
	    parameterization = MONTY12;
	    argv++;
	    argc--;
	}
	else 
	  usage (programname);
    }

  if (argc < 5)
    usage (programname);

  Pmin = atol (argv[1]);
  Pmax = atol (argv[2]);
  B1 = atol (argv[3]);
  sigma = atol (argv[4]);

  for (p = 2.0; p < Pmin; p = getprime (p));

  for ( ; p <= Pmax; p = getprime (p))
    {
      primes ++;
      if (verbose)
	{
	  printf ("Testing p = %ld", (unsigned long) p);
	}

      if (ecm ((unsigned long) p, (double) B1, sigma, parameterization))
	{
	  found ++;
	  if (verbose)
	    printf (": smooth");
	}

      if (do_count)
	{
	  nr_points = ell_curveorder ((unsigned long) p, sigma, 
	                              parameterization);
	  if (verbose)
	      printf (", order = %ld", nr_points);
	  if (nr_points > 0)
	    {
	      pow2 += prime_exponent (nr_points, 2UL);
	      pow3 += prime_exponent (nr_points, 3UL);
	      pow5 += prime_exponent (nr_points, 5UL);
	      pow7 += prime_exponent (nr_points, 7UL);
	    }
	}
      if (verbose)
	  printf ("\n");
    }
  getprime (0.0);

  if (do_count)
    {
      printf ("Avg. exponent of 2: %f, 3: %f, 5: %f, 7: %f\n",
	      (double)pow2/(double)primes, (double)pow3/(double)primes, 
	      (double)pow5/(double)primes, (double)pow7/(double)primes);
    }

  if (gnuplot)
      printf ("%lu %f\n", Pmin, (double) found / (double) primes);
  else
      printf ("primes=%lu found=%lu (%1.0f%%)\n", primes, found,
	      100.0 * (double) found / (double) primes);
  
  return 0;
}
