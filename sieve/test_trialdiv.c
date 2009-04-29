#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "ularith.h"
#include "modredc_ul.h"
#include "trialdiv.h"
#include "getprime.h"
#include <sys/time.h>
#include <sys/resource.h>

int main (int argc, char **argv)
{
  fbprime_t *primes;
  trialdiv_divisor_t *d;
  const size_t max_div = 16;
  unsigned long *factors;
  int i, r = 0, expect = 0, len = 1, nr_primes = 1000, nr_N = 1000000;
  mpz_t M, N;
  struct rusage usage;
  double usrtime;
  
  if (argc > 1)
    len = atoi (argv[1]);
  
  if (argc > 2)
    nr_N = atoi (argv[2]);

  if (argc > 3)
    nr_primes = atoi (argv[3]);
  
  if (len > TRIALDIV_MAXLEN)
    {
      printf ("Error, trialdiv not implemented for input sizes greater than "
	      "%d words\n", TRIALDIV_MAXLEN);
      exit (EXIT_FAILURE);
    }

  mpz_init (M);
  mpz_init (N);
  mpz_set_ui (N, 1UL);
  mpz_mul_2exp (N, N, 8 * sizeof(unsigned long) * len);
  mpz_tdiv_q_ui (N, N, 3UL);
  mpz_nextprime (N, N);

  for (i = 0; i < 0; i++) /* To skip some primes */
    getprime(1UL);
  primes = malloc (nr_primes * sizeof (fbprime_t));
  for (i = 0; i < nr_primes; i++)
    {
      unsigned long ppow;
      primes[i] = getprime(1UL);
      ppow = 1UL;
      /* FIXME powers greater than ULONG_MAX may appear in N */
      while (ULONG_MAX / primes[i] >= ppow)
	{
	  ppow *= primes[i];
	  expect += nr_N / ppow;
	  if (mpz_tdiv_ui (N, ppow) + nr_N % ppow > ppow)
	    expect++;
	}
    }
  printf ("Trial dividing by primes in [%lu, %lu]\n", 
	  (unsigned long) primes[0], (unsigned long) primes[nr_primes - 1]);
  getprime (0);

  d = trialdiv_init (primes, nr_primes);
  free (primes);
  factors = (unsigned long *) malloc (max_div * sizeof (unsigned long));
  
  for (i = 0; i < nr_N; i++)
    {
      size_t t;
      mpz_set (M, N);
      do {
	t = trialdiv (factors, M, d, max_div);
	r += (t > max_div) ? max_div : t;
      } while (t > max_div);
      mpz_add_ui (N, N, 1UL);
    }
  
  getrusage(RUSAGE_SELF, &usage);
  usrtime = (double) usage.ru_utime.tv_sec * 1000000. +
      (double) usage.ru_utime.tv_usec;

  mpz_clear (M);
  mpz_clear (N);
  trialdiv_clear (d);
  free (factors);
  printf ("Found %d divisors, expected %d\n", r, expect);
  printf ("Time: %f s, per N: %f mus, per prime: %f ns\n", 
	  usrtime / 1e6, usrtime / nr_N, usrtime / nr_N / nr_primes * 1e3);
  if (r != expect)
    {
      printf ("Error: did not find the expected number of divisors!\n");
      exit (EXIT_FAILURE);
    }
  return 0;
}
