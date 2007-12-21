/* check the relations produces by the siever, and factors the residues.
   Usage: checknorms -poly c80.poly c80.rels1 c80.rels2 ...
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "cado.h"
#include "ecm.h"
#include "utils/utils.h"

/* number of Miller-Rabin tests in mpz_probab_prime_p */
#define REPS 1

#define MAXREL_LENGTH 1024 /* max length in characters of a relation */

/* norm <- |g(a,b)| = |g[1]*a + g[0]*b| */
void
eval_rational (mpz_t norm, mpz_t *g, long a, unsigned long b)
{
  mpz_mul_si (norm, g[1], a);
  mpz_addmul_ui (norm, g[0], b);
  mpz_abs (norm, norm);
}

/* norm <- |f(a,b)| = |f[d]*a^d + ... + f[0]*b^d| */
void
eval_algebraic (mpz_t norm, mpz_t *f, int d, long a, unsigned long b)
{
  mpz_t B; /* powers of b */
  mpz_init_set_ui (B, 1);
  mpz_set (norm, f[d]);
  while (d-- > 0)
    {
      mpz_mul_si (norm, norm, a);
      mpz_mul_ui (B, B, b);
      mpz_addmul (norm, f[d], B);
    }
  mpz_clear (B);
  mpz_abs (norm, norm);
}

/* check sb < p */
void
check_prime (mpz_t p, unsigned long sb, long a, unsigned long b)
{
  static unsigned long reports = 0;
  if (mpz_cmp_ui (p, sb) <= 0)
    {
      if (++reports <= 10)
	gmp_fprintf (stderr, "Warning, (a=%ld,b=%lu): prime %Zd smaller than factor base bound %lu\n", a, b, p, sb);
    }
}

// Get a prime factor of norm and put it in p1.
void
find_prime_factor(mpz_t p1, mpz_t norm)
{
    mpz_t r;
    int res;

    /* stupid? */
    mpz_init(r);
    mpz_sqrtrem(p1, r, norm);
    res = mpz_cmp_ui(r, 0) == 0;
    mpz_clear(r);
    if(res){
	fprintf(stderr, "Perfect square!!!\n");
	return;
    }
    /* the following seems to work well for a large prime bound of 2^24,
       i.e., residues less than 2^48 */
    double B1 = 4.0 * (double) mpz_sizeinbase (norm, 2);
    do
	{
	    //	    gmp_fprintf (stderr, "ECM on %Zd\n", norm);
	    res = ecm_factor (p1, norm, B1, NULL);
	    //	    fprintf(stderr, "Used ecm_factor: %d\n", res);
	    B1 += sqrt (B1);
	}
    while (res == 0 || mpz_cmp (p1, norm) == 0 
	   || mpz_probab_prime_p (p1, REPS) == 0);
}

/* factor norm into at most 2 primes smaller than 2^lp:
   0) if norm = 1, return 0
   1) if norm is prime, norm < 2^lp, put it into p1 and return 1
   2) if norm = p1*p2 with p1,p2 primes < 2^lp, return 2
   3) otherwise if norm has 3 prime factors of more,
      or if one prime factor is >= 2^lp, return 3.

   TODO: check that the case p1 == p2 (rare) doesn't destroy anything

*/
unsigned long
factor (mpz_t p1, mpz_t p2, mpz_t norm, size_t lp)
{
  static unsigned long count_too_large_factor = 0;

  if (mpz_cmp_ui (norm, 1) == 0)
    return 0;
  else if (mpz_probab_prime_p (norm, REPS))
    {
      if (mpz_sizeinbase (norm, 2) > lp) /* this happens very frequently */
        return 3;
      mpz_set (p1, norm);
      return 1;
    }
  else /* two or more factors */
    {
      find_prime_factor(p1, norm);
      /* p1 is a non trivial prime factor of norm */
      /* we want p1 < 2^lp */
      if (mpz_sizeinbase (p1, 2) > lp) /* this should be rare */
        {
          if (count_too_large_factor++ < 10)
            gmp_fprintf (stderr, "Warning, composite norm %Zd has too large prime factor %Zd\n", norm, p1);
          return 3;
        }
      mpz_divexact (p2, norm, p1);
      /* we also want p2 < 2^lp */
      if (mpz_sizeinbase (p2, 2) > lp)
        {
          if (count_too_large_factor++ < 10)
            gmp_fprintf (stderr, "Warning, composite norm %Zd has too large factor %Zd\n", norm, p2);
          return 3;
        }
      if (mpz_probab_prime_p (p1, REPS) && mpz_probab_prime_p (p2, REPS))
	return 2;
      else
	{
          static unsigned long count = 0;
          if (count++ < 10)
            gmp_fprintf (stderr, "Warning, norm with 3 primes or more: %Zd\n", norm);
	  return 3;
	}
    }
}

int
coprime (long a, unsigned long b, mpz_t t, mpz_t u)
{
  mpz_set_si (t, a);
  mpz_set_ui (u, b);
  mpz_gcd (t, t, u);
  return mpz_cmp_ui (t, 1) == 0;
}

/* primes[0]..primes[nprimes-1] are the small primes omitted in relations */
unsigned long
checkrels (char *f, cado_poly cpoly, int verbose, size_t mfbr, size_t mfba,
           unsigned long *primes, unsigned long nprimes)
{
  FILE *fp;
  unsigned long b, nb, mfba_reports = 0, mfbr_reports = 0;
  unsigned long rels = 0; /* number of relations read */
  unsigned long out_rels = 0; /* number of relations output */
  fbprime_t p;
  char s[MAXREL_LENGTH];
  unsigned long len; /* length of current relation */
  long a;
  int c, ok;
  mpz_t norm, p1, p2;
  unsigned long i;

  mpz_init (norm);
  mpz_init (p1);
  mpz_init (p2);
  if (verbose)
    fprintf (stderr, "Checking relations from file %s\n", f);
  fp = fopen (f, "r");
  if (fp == NULL)
    {
      fprintf (stderr, "Error: unable to open file %s\n", f);
      exit (EXIT_FAILURE);
    }
  while (!feof (fp))
    {
      len = 0;
      ok = 0; /* ok=0 for a comment line, ok>=1 for a relation,
		 ok>=3 for a valid relation */
      nb = 0;
      /* skip comment lines */
      c = getc (fp);
      if (c == '#') /* read entire line */
	goto end;
      else
	ungetc (c, fp);
      c = fscanf (fp, "%ld,%lu:", &a, &b);
      if (c == 0 || c == EOF)
	break;

      ok = 1;
      len += sprintf (s, "%ld,%lu:", a, b);

      /* check a and b are coprime */
      if (coprime (a, b, p1, p2) == 0)
        {
          static unsigned long report = 0;
          if (report ++ < 10)
            fprintf (stderr, "Warning, discard relation (%ld,%lu) since a,b are not coprime\n", a, b);
          //          exit (EXIT_FAILURE);
          goto end;
        }

      /* evaluate norm on rational side */
      eval_rational (norm, cpoly->g, a, b);
      /* read primes on rational side */
      while (fscanf (fp, "%x", &p) != 0)
	{
	  /* check p divides the norm */
	  if (mpz_fdiv_q_ui (norm, norm, p) != 0)
	    {
	      fprintf (stderr, "Error, prime %x=%u does not divide norm on rational side for a=%ld, b=%lu\n", p, p, a, b);
	      exit (EXIT_FAILURE);
	    }
          /* check p is smaller than the large prime bound */
          if (p > 1UL << cpoly->lpbr)
            {
              fprintf (stderr, "Warning, rational prime %x exceeds large prime bound 2^%d\n", p, cpoly->lpbr);
              nb = 3; /* discard relation */
            }
	  len += sprintf (s + len, "%x", p);
	  c = getc (fp);
	  if (c == ':')
	    break;
	  len += sprintf (s + len, "%c", c);
	}
      /* divide by small primes */
      for (i = 0; i < nprimes; i++)
        while (mpz_divisible_ui_p (norm, primes[i]))
          {
            len += sprintf (s + len, ",%lx", primes[i]);
            mpz_divexact_ui (norm, norm, primes[i]);
          }
      /* check the residue is smaller than 2^mfbr */
      if (mpz_sizeinbase (norm, 2) > mfbr)
	{
	  if (++mfbr_reports <= 10)
	    gmp_fprintf (stderr, "Warning, rat. residue %Zd exceeds bound for a=%ld, b=%lu\n", norm, a, b);
	  nb = 3; /* discard relation */
	}
      else if (nb != 3) /* factor residue on rational side */
	nb = factor (p1, p2, norm, cpoly->lpbr);
      /* write additional primes */
      if (nb <= 2)
	{
	  ok ++;
	  if (nb >= 1)
	    {
	      check_prime (p1, cpoly->rlim, a, b);
	      len += gmp_sprintf (s + len, ",%Zx", p1);
	    }
	  if (nb == 2)
	    {
	      check_prime (p2, cpoly->rlim, a, b);
	      len += gmp_sprintf (s + len, ",%Zx", p2);
	    }
	}
      len += sprintf (s + len, ":");

      /* warning: there might be no prime on the algebraic side */
      c = getc (fp);
      if (c != '\n')
        ungetc (c, fp);

      /* evaluate norm on algebraic side */
      eval_algebraic (norm, cpoly->f, cpoly->degree, a, b);
      /* read primes on algebraic side */
      while (c != '\n' && fscanf (fp, "%x", &p) != 0)
	{
	  /* check p divides the norm */
	  if (mpz_fdiv_q_ui (norm, norm, p) != 0)
	    {
	      fprintf (stderr, "Error, prime %x=%u does not divide norm on algebraic side for a=%ld, b=%lu\n", p, p, a, b);
	      exit (EXIT_FAILURE);
	    }
          /* check p is smaller than the large prime bound */
          if (p > 1UL << cpoly->lpba)
            {
              fprintf (stderr, "Warning, algebraic prime %x exceeds large prime bound 2^%d\n", p, cpoly->lpba);
              nb = 3; /* discard relation */
            }
	  len += sprintf (s + len, "%x", p);
	  c = getc (fp);
	  if (c != ',')
	    break;
          len += sprintf (s + len, "%c", c);
	}
      /* divide by small primes */
      for (i = 0; i < nprimes; i++)
        while (mpz_divisible_ui_p (norm, primes[i]))
          {
            len += sprintf (s + len, ",%lx", primes[i]);
            mpz_divexact_ui (norm, norm, primes[i]);
          }
      /* check the residue is smaller than 2^mfba */
      if (nb == 3)
	nb = 3; /* relations already discarded by rational side */
      else if (mpz_sizeinbase (norm, 2) > (unsigned) mfba)
	{
	  if (++mfba_reports <= 10)
	    gmp_fprintf (stderr, "Warning, alg. residue %Zd exceeds bound for a=%ld, b=%lu\n", norm, a, b);
	  nb = 3; /* discard relation */
	}
      else /* factor residue on algebraic side */
	nb = factor (p1, p2, norm, cpoly->lpba);
      /* write additional primes */
      if (nb <= 2)
	{
	  ok ++;
	  if (nb >= 1)
	    {
	      check_prime (p1, cpoly->alim, a, b);
	      len += gmp_sprintf (s + len, ",%Zx", p1);
	    }
	  if (nb == 2)
	    {
	      check_prime (p2, cpoly->alim, a, b);
	      len += gmp_sprintf (s + len, ",%Zx", p2);
	    }
	}

      rels ++;
      if (ok >= 3){
	out_rels ++;
	if(!(out_rels % 100000))
	    fprintf(stderr, "%lu-th relation found // %lu read\n",
		    out_rels, rels);
      }
    end:
      /* read to end of line */
      while (c != '\n' && c != EOF)
	{
	  len += sprintf (s + len, "%c", c);
	  c = getc (fp);
	}
      s[len] = '\0';
      if (ok == 0 || ok >= 3)
	printf ("%s\n", s);
    }
  fclose (fp);
  mpz_clear (norm);
  mpz_clear (p1);
  mpz_clear (p2);
  fprintf (stderr, "File %s: output %lu out of %lu relations\n", f, out_rels,
	   rels);
  return out_rels;
}

int
main (int argc, char *argv[])
{
  int verbose = 0;
  char *polyfilename = NULL;
  cado_poly cpoly;
  unsigned long tot_rels = 0;
  int mfbr = 0, mfba = 0, nb_files, i;
  unsigned long tdmax = 1; /* bound for trial division */
  unsigned long *primes; /* small primes omitted in relations */
  int nprimes;
  unsigned long p;

  fprintf (stderr, "%s.r%s", argv[0], REV);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", argv[i]);
  fprintf (stderr, "\n");

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 1 && strcmp (argv[1], "-v") == 0)
	{
	  verbose++;
	  argc--;
	  argv++;
	}
      else if (argc > 2 && strcmp (argv[1], "-poly") == 0)
	{
	  polyfilename = argv[2];
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-mfbr") == 0)
	{
	  mfbr = atoi (argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-mfba") == 0)
	{
	  mfba = atoi (argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-t") == 0)
	{
	  tdmax = atoi (argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else 
	{
	  fprintf (stderr, "Usage: %s [-v] [-mfbr <n>] [-mfba <n>] [-t <n>] -poly <file> rels1 ... relsn\n", argv[0]);
	  exit (EXIT_FAILURE);
	}
    }

  if (polyfilename == NULL)
    {
      fprintf (stderr, "Please specify a polynomial file with -poly\n");
      exit (EXIT_FAILURE);
    }

  if (!read_polynomial (cpoly, polyfilename))
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }

  /* default residue bounds are those from file, if not overridden by
     command line parameters */
  if (mfbr == 0)
    mfbr = cpoly->mfbr;
  if (mfba == 0)
    mfba = cpoly->mfba;

  /* compute small primes */
  primes = (unsigned long*) malloc (tdmax * sizeof(unsigned long));
  nprimes = 0;
  for (p = 2; p < tdmax; p += 1 + (p & 1))
    {
      for (i = 0; (i < nprimes) && (primes[i] * primes[i] <= p); i++)
        if (p % primes[i] == 0)
          {
            i = nprimes;
            break;
          }
      if (nprimes == 0 || i < nprimes)
        primes[nprimes++] = p;
    }
  primes = realloc (primes, nprimes * sizeof(unsigned long));

  nb_files = argc - 1;
  while (argc > 1)
    {
      tot_rels += checkrels (argv[1], cpoly, verbose, mfbr, mfba, primes,
                             nprimes);
      argc --;
      argv ++;
    }

  free (primes);

  if (nb_files > 1)
    fprintf (stderr, "All input files: output %lu relations\n", tot_rels);

  return 0;
}
