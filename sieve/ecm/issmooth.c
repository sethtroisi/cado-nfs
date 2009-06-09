/* Reads integers pairs "p N" from stdin and checks them for smoothness, 
   where smooth is defined as
     For given B1, B2, k, d,
     N is smooth if N/gcd(N,E) is 1 or is a prime B1 < p <= B2,
     or is coprime to d and <d/2, or divides d,
     where E = k * lcm(1, 2, 3, ..., B1).
   By default, k=1 and d=1, so that they have no effect.
   With -pmin and -pmax, test only those input where pmin <= p <= pmax
   Prints the number of smooth numbers, and with -v option each "p N" 
   wher N is smooth
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <gmp.h>

static int 
isprime (const unsigned long N, mpz_t t)
{
  mpz_set_ui (t, N);
  return mpz_probab_prime_p (t, 1);
}


static unsigned long
gcd (unsigned long a, unsigned long b)
{
  unsigned long t;

  if (a >= b)
    a %= b;

  while (a > 0)
    {
      t = b % a;
      b = a;
      a = t;
    }

  return b;
}


int 
main (int argc, char **argv)
{
  mpz_t E, m_c;
  unsigned long N, c, p, B1, B2, k = 1, d = 1, i, j, smooth = 0, 
                pmin = 0, pmax = ULONG_MAX;
  int quiet = 0, verbose = 0;
  char *d_coprime;
  
  /* Parse arguments */
  if (argc < 3)
    {
      printf ("Please provide B1 B2 (and optionally k) values\n");
      exit (EXIT_FAILURE);
    }

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argv[1][1] == 'v')
        {
          verbose = 1;
          argc--;
          argv++;
        }
      else if (argv[1][1] == 'q')
        {
          quiet = 1;
          argc--;
          argv++;
        }
      else if (argc > 2 && argv[1][1] == 'd')
        {
          d = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-pmin") == 0)
        {
          pmin = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && strcmp (argv[1], "-pmax") == 0)
        {
          pmax = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else
        {
          fprintf (stderr, "Unknown option %s\n", argv[1]);
          exit (EXIT_FAILURE);
        }
    }

  B1 = strtoul (argv[1], NULL, 10);
  B2 = strtoul (argv[2], NULL, 10);
  if (argc > 3)
    k = strtoul (argv[3], NULL, 10);
  d_coprime = (char *) malloc (d);
  for (i = 1; i < d; i++)
    d_coprime[i] = (gcd(i,d) == 1) ? 1 : 0;
  
  /* Compute E */
  mpz_init (E);
  mpz_set_ui (E, 1UL);
  for (i = 2UL; i <= B1; i = (i + 1UL) | 1UL)
    if (mpz_gcd_ui (NULL, E, i) == 1UL) /* A prime? */
      {
        for (j = 1UL; j <= B1 / i; j *= i);
        mpz_mul_ui (E, E, j);
      }
  mpz_mul_ui (E, E, k);
  if (!quiet)
    gmp_fprintf (stderr, "B1 = %lu, B2 = %lu, k = %lu, pmin = %lu, "
                 "pmax = %lu, E = %Zd\n", B1, B2, k, pmin, pmax, E);
  
  /* Test input numbers */
  mpz_init (m_c);
  while (scanf("%lu %lu\n", &p, &N) == 2)
    {
      if (pmin <= p && p <= pmax)
        {
          c = N / mpz_gcd_ui (NULL, E, N);
          if (c == 1UL || 
              d % c == 0UL || (c < d / 2 && d_coprime[c] == 1) ||
              (B1 < c && c <= B2 && isprime(c, m_c)))
            {
              smooth++;
              if (verbose)
                printf ("%lu %lu\n", p, N);
            }
        }
    }
  mpz_clear (E);
  mpz_clear (m_c);

  if (!quiet)
    printf ("%lu smooth numbers\n", smooth);
  
  return 0;
}
