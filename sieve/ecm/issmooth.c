/* Reads integers pairs "p O N" from stdin and checks N for smoothness, 
   where smooth is defined as
     For given B1, B2, k, d,
     N is smooth if N/gcd(N,E) is 1 or is a prime B1 < p <= B2,
     or is coprime to d and <d/2, or divides d,
     where E = k * lcm(1, 2, 3, ..., B1).
   By default, k=1 and d=1, so that they have no effect.
   With -pmin and -pmax, test only those input where pmin <= p <= pmax
   Prints the number of smooth numbers, and with -v option each "p O N" 
   where N is smooth. 
   Also prints the average exponent of primes up to 19 in O/N 
   (which requires N|O).
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <gmp.h>

#define NR_EXPONENTS 8
#define EXP_PRIMES {2,3,5,7,11,13,17,19}

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

static unsigned int
valuation (const unsigned long n, const unsigned long p)
{
  unsigned long c = n;
  unsigned int i = 0;
  while (c % p == 0)
    {
      c /= p;
      i++;
    }
  return i;
}


int 
main (int argc, char **argv)
{
  mpz_t E, m_c;
  unsigned long N, O, c, p, i, j;
  /* Parameters */
  unsigned long B1, *B2, maxB2, k = 1, d = 1, pmin = 0, pmax = ULONG_MAX, m = 1;
  int quiet = 0, verbose = 0;
  unsigned int nr_B2 = 0;
  /* Counters */
  unsigned long input = 0, *smooth = 0; 
  int scanf_code;
  char *d_coprime;
  /* Input and smoothness count for each residue class modulo m */
  unsigned long *res_input, *res_smooth; 
  char buf[256];
  unsigned int exp_order[NR_EXPONENTS], exp_index[NR_EXPONENTS];
  const unsigned long exp_primes[NR_EXPONENTS] = EXP_PRIMES;
  
  /* Parse arguments */
  if (argc < 3)
    {
      printf ("Please provide B1 and B2 values\n");
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
      else if (argc > 2 && argv[1][1] == 'k')
        {
          k = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && argv[1][1] == 'd')
        {
          d = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
        }
      else if (argc > 2 && argv[1][1] == 'm')
        {
          m = strtoul (argv[2], NULL, 10);
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
  nr_B2 = argc - 2;
  B2 = malloc (nr_B2 * sizeof (unsigned long));
  smooth = malloc (nr_B2 * sizeof (unsigned long));
  maxB2 = 0;
  for (i = 0; i < nr_B2; i++)
    {
      smooth[i] = 0;
      B2[i] = strtoul (argv[2 + i], NULL, 10);
      maxB2 = (B2[i] > maxB2) ? B2[i] : maxB2;
    }

  d_coprime = (char *) malloc (d);
  for (i = 1; i < d; i++)
    d_coprime[i] = (gcd(i,d) == 1) ? 1 : 0;
  for (i = 0; i < NR_EXPONENTS; i++)
    exp_order[i] = exp_index[i] = 0;
  
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
    {
      fprintf (stderr, "B1 = %lu, ", B1);
      for (i = 0; i < nr_B2; i++)
        fprintf (stderr, "B2[%lu] = %lu, ", i, B2[i]);
      gmp_fprintf (stderr, "k = %lu, pmin = %lu, pmax = %lu, E = %Zd\n", 
                   k, pmin, pmax, E);
    }
  
  res_input = (unsigned long *) malloc (m * sizeof(unsigned long));
  res_smooth = (unsigned long *) malloc (m * sizeof(unsigned long));
  for (i = 0; i < m; i++)
    res_input[i] = res_smooth[i] = 0;
  
  /* Test input numbers */
  mpz_init (m_c);
  while (!feof(stdin))
    {
      if (fgets (buf, 256, stdin) == NULL)
        break;
      scanf_code = sscanf (buf, "%lu %lu %lu\n", &p, &O, &N);
      if (pmin <= p && p <= pmax)
        {
          /* input line may contain only "N", or "p N", or "p O N" */
          input++;
          if (scanf_code == 1)
            N = p;
          else if (scanf_code == 2)
            N = O;
          res_input[p % m]++;
          c = N / mpz_gcd_ui (NULL, E, N);
          if (c == 1UL || 
              d % c == 0UL || (c < d / 2 && d_coprime[c] == 1) ||
              (B1 < c && c <= maxB2 && isprime(c, m_c)))
            {
              for (i = 0; i < nr_B2; i++)
                smooth[i] += (c <= B2[i]) ? 1UL : 0UL;
              res_smooth[p % m]++;
              if (verbose)
                printf ("%lu %lu\n", p, N);
            }
          if (scanf_code == 3)
            {
              assert (O % N == 0);
              for (i = 0; i < NR_EXPONENTS; i++)
                {
                  exp_order[i] += valuation (O, exp_primes[i]);
                  exp_index[i] += valuation (O/N, exp_primes[i]);
                }
            }
        }
    }
  mpz_clear (E);
  mpz_clear (m_c);

  if (!quiet)
    {
      printf ("%lu input numbers in [%lu, %lu]\n", input, pmin, pmax);
      for (i = 0; i < nr_B2; i++)
        printf ("%lu %lu,%lu-smooth numbers, ratio %f\n", 
              smooth[i], B1, B2[i], (double) smooth[i] / (double) input);
      if (m > 1UL)
        {
          printf ("Number of smooth/input numbers and ratio where p == r (mod %lu):\n", m);
          for (i = 0; i < m; i++)
            if (res_input[i] != 0)
              printf ("r = %lu: %lu/%lu = %f  ", 
                      i, res_smooth[i], res_input[i], 
                      (double) res_smooth[i] / (double) res_input[i]);
          printf ("\n");
        }
      printf ("Sum of valuation of q in O:\n");
      for (i = 0; i < NR_EXPONENTS; i++)
        if (exp_order[i])
          printf ("q = %lu: %u  ", exp_primes[i], exp_order[i]);
      printf ("\n");
      printf ("Sum of valuation of q in N/O (index of starting element):\n");
      for (i = 0; i < NR_EXPONENTS; i++)
        if (exp_index[i])
          printf ("q = %lu: %u  ", exp_primes[i], exp_index[i]);
      printf ("\n");
    }
  
  return 0;
}
