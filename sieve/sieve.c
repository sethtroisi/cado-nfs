#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#ifdef HAVE_MSRH
#include <asm/msr.h>
#else
#define rdtscll(x)
#endif
#include "cado.h"
#include "fb.h"
#include "mod_ul.c"
#include "sieve_aux.h"
#include "basicnt.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define REFAC_PRP_SIZE_THRES 50
#define REFAC_PRP_THRES 500
#define SIEVE_PERMISSIBLE_ERROR 5
#define PRP_REPS 1 /* number of tests in mpz_probab_prime_p */

/* Some multiple precision functions we'll need */

void
mp_poly_eval (mpz_t r, mpz_t *poly, int deg, long a)
{
  int i;

  mpz_set (r, poly[deg]);
  for (i = deg - 1; i >= 0; i--)
    {
      mpz_mul_si (r, r, a);
      mpz_add (r, r, poly[i]);
    }
}


/* Scale coefficient f_i by a^i (inv == 1) or by c^(deg-i) (inv == -1) */
void
mp_poly_scale (mpz_t *r, mpz_t *poly, const int deg, const long c, 
	       const int inv)
{
  int i;
  mpz_t t;

  ASSERT (inv == 1 || inv == -1);

  mpz_init (t);
  mpz_set_ui (t, 1);
  for (i = 0; i <= deg; i++)
    {
      const int j = (inv == 1) ? i : deg - i;
      mpz_mul (r[j], poly[j], t);
      mpz_mul_si (t, t, c);
    }
  
  mpz_clear (t);
}


void 
mp_poly_print (mpz_t *poly, int deg, const char *name, int homogeneous)
{
  int i;
  printf ("%s", name);

  if (!homogeneous)
    {
      for (i = deg; i >= 2; i--)
	{
	  if (mpz_sgn (poly[i]) != 0)
	    gmp_printf ("%s%Zd * x^%d ", 
			mpz_sgn(poly[i]) > 0 ? "+" : "", poly[i], i);
	}
      if (deg >= 1 && mpz_sgn (poly[1]) != 0)
	gmp_printf ("%s%Zd * x ", mpz_sgn(poly[1]) > 0 ? "+" : "", poly[1]);
      if (deg >= 0 && mpz_sgn (poly[0]) != 0)
	gmp_printf ("%s%Zd", mpz_sgn(poly[0]) > 0 ? "+" : "", poly[0]);
    } 
  else 
    {
      for (i = deg; i >= 0; i--)
	{
	  if (mpz_sgn (poly[i]) != 0)
	    {
	      gmp_printf (" %s%Zd", 
			  mpz_sgn(poly[i]) > 0 && i < deg ? "+" : "", 
			  poly[i]);
	      if (i > 1)
		printf (" *a^%d", i);
	      if (i == 1)
		printf ("*a");
	      if (i + 1 < deg)
		printf ("*b^%d", deg - i);
	      if (i + 1 == deg)
		printf ("*b");
	    }
	}
    }
}


/* Used only in compute_norms() */
static unsigned char
log_norm (const double *f, const int deg, const double x, 
	  const double log_scale, const double log_proj_roots)
{
  double r;
  int i;
  unsigned char n;

  r = f[deg];
  for (i = deg - 1; i >= 0; i--)
    r = r * x + f[i];

  n = fb_log (fabs (r), log_scale, - log_proj_roots);
  /* printf ("Norm at x = %.0f is %.0f, rounded log is %d\n", x, r, (int) n); */
  return n;
}


/* Very slow but thorough way of computing norms. Returns the maximum rounded
   log norm in sievearray. */

unsigned char
compute_norms (unsigned char *sievearray, const long amin, const long amax, 
	       const unsigned long b, const double *poly, const int deg, 
	       const double proj_roots, const double log_scale, const int odd, 
	       const int verbose)
{
  double f[MAXDEGREE + 1]; /* Polynomial in $a$ for a given fixed $b$ */
  double bpow;
  unsigned long long tsc1, tsc2;
  const double log_proj_roots = log(proj_roots) * log_scale;
  long a, a2;
  int i;
  unsigned char n1, n2, nmax;
  const int stride = 128;
  
  ASSERT (odd == 0 || odd == 1);
  ASSERT (!odd || (amin & 1) == 1);
  ASSERT (!odd || (amax & 1) == 1);
  ASSERT (amin <= amax);

  rdtscll (tsc1);
  
  bpow = 1.;
  for (i = deg; i >= 0; i--)
    {
      f[i] = poly[i] * bpow;
      bpow *= (double) b;
    }
  
  n1 = log_norm (f, deg, (double) amin, log_scale, log_proj_roots);
  nmax = n1;

  a = amin;
  a2 = amin;
  /* a == amin + (d << odd), d = (a - amin) >> odd */

  while (a <= amax)
    {
      /* We'll cover [a ... a2[ now */
      ASSERT (a == a2);
      ASSERT (n1 == log_norm (f, deg, (double) a, log_scale, log_proj_roots));
      a2 = MIN(a + stride, amax + (1 << odd));
      ASSERT (!odd || (a2 - a) % 2 == 0);
      
      n2 = log_norm (f, deg, (double) a2, log_scale, log_proj_roots);
      
      if (n1 == n2) 
	{
	  /* Let's assume the log norm is n1 everywhere in this interval */
	  memset (sievearray + ((a - amin) >> odd), n1, (a2 - a) >> odd);
	}
      else
	{
	  /* n1 and n2 are different. Do each a up to (exclusive) a2 
	     individually */
	  /* printf ("log_c(F(%ld, %lu)) == %d != log_c(F(%ld, %lu)) == %d\n",
	     a, b, n1, a2, b, n2); */
	  sievearray[(a - amin) >> odd] = n1;
	  a += 1 << odd;
	  for ( ; a < a2; a += 1 << odd)
	    {
	      unsigned char n = log_norm (f, deg, (double) a, log_scale, 
					  log_proj_roots);
	      sievearray[(a - amin) >> odd] = n;
	      if (n > nmax)
		nmax = n;
	    }
	  
	}
      a = a2;
      n1 = n2;
    }
  
  rdtscll (tsc2);
#ifdef HAVE_MSRH
  if (verbose)
    {
      printf ("# Computing norms took %lld clocks\n", tsc2 - tsc1);
      printf ("# Maximum rounded log norm is %u\n", (unsigned int) nmax);
    }
#endif
  
#ifdef PARI
  for (a = amin; a <= amax; a += 1 << odd)
    printf ("if(log_c(%c(%ld, %lu), %f) != %d,"
	    "print (\"log_c(%c(\", %ld,\", \", %lu\"), %f) \", %d)) /* PARI */\n", 
	    deg > 1 ? 'F' : 'G', a, b, log_proj_roots, sievearray[(a - amin) >> odd],
	    deg > 1 ? 'F' : 'G', a, b, log_proj_roots, sievearray[(a - amin) >> odd]); 
#endif

  return nmax;
}


static void 
sieve_small_slow (unsigned char *sievearray, factorbase_small_inited_t *fb,
		  const unsigned int arraylen)
{
  fbprime_t p;
  unsigned int d;
  unsigned char l;

  for (; fb->p > 0; fb++)
    {
      p = fb->p;
      d = fb->loc_and_log & 0xffffff; /* Mask low 24 bits */
      l = fb->loc_and_log >> 24;
      while (d < arraylen)
	{
	  sievearray[d] -= l;
	  d += p;
	}
      d -= arraylen;
      fb->loc_and_log = d + (l << 24);
    }
}


/* sievearray must be long enough to hold amax-amin+1 chars.
   proj_roots is the rounded log of the projective roots */

/* odd: if 1, the sieve array contains only locations for odd a, 
        amin <= a <= amax */

/* Eventually we'll do:
   1. Generate bucket-sorted sieve updates for large fb primes. One bucket
      should hold the updates to positions in the same L2 chunk
   2. For each L2 chunk in this line:
     2.1 For each L1 chunk in this L2 chunk:
       2.1.1 Compute norms, dividing out tiny primes, i.e. p < 10
       2.1.2 Sieve by small primes, p < L1SIZE
     2.2 Do the sieve updates from this L2 chunk's bucket, writing
         sieve reports with useful primes and approx. log
*/


void 
sieve (unsigned char *sievearray, factorbase_degn_t *fb, 
       const long amin, const long amax, const unsigned long b, 
       const unsigned char threshold, fbprime_t *useful_primes_par,
       const unsigned int useful_length_par,
       const int odd)
{
  /* The sievearray[0] entry corresponds to (amin, b), and
     sievearray[d] to (amin + d * (1 + odd), b) */
  fbprime_t *useful_primes = useful_primes_par;
  const uint32_t l = (amax - amin) / (1 + odd) + 1;
  uint32_t i, amin_p, p, d;
  unsigned int useful_length = useful_length_par;
  unsigned char plog;

  ASSERT (odd == 0 || odd == 1);
  ASSERT (!odd || (amin & 1) == 1);
  ASSERT (!odd || (amax & 1) == 1);
  ASSERT (!odd || (b & 1) == 0);
  ASSERT (amin <= amax);

  if (useful_length > 0)
    useful_length--; /* Make sure we have space for a 0 mark at the end */

  /* Do the sieving */

  while (fb->p > 0)
    {
      ASSERT (fb_entrysize (fb) <= fb->size);
      p = fb->p;
      plog = fb->plog;

      /* Compute amin % p for the a value of the sieve location in 
	 sievearray[0] */
      /* FIXME This modular reduction should be simplified somehow. Do it 
	 once and store it in fb? */

      amin_p = signed_mod_longto32 (amin, p);
      
      for (i = 0; i < fb->nr_roots; i++)
        {
	  d = first_sieve_loc (p, fb->roots[i], amin_p, b, odd);
	  
	  /* Now d is the first index into sievearray where p divides. */

          /* Now update the sieve array. There is some blocking, but no
             bucket sieving atm. */
          for (; d < l; d += p)
	    {
	      unsigned char k;
	      k = sievearray[d] - plog;
	      sievearray[d] = k;
	      if (k <= threshold && useful_length > 0)
		  {
		    *useful_primes++ = p;
		    useful_length--;
		  }
	    }
        }
      
      /* Move on to the next factor base prime */
      fb = fb_next (fb);
    }
  if (useful_length_par > 0)
    *useful_primes++ = 0;
}

#define TYPE_STRING 0
#define TYPE_MPZ 1
#define TYPE_INT 2
#define TYPE_ULONG 3
#define TYPE_DOUBLE 4
#define PARSE_MATCH -1
#define PARSE_ERROR 0
#define PARSE_NOMATCH 1

/* Return value: 1: tag does not match, 0: error occurred, -1 tag did match */

static int
parse_line (void *target, char *line, const char *tag, int *have, 
	    const int type)
{
  char *lineptr = line;

  if (strncmp (lineptr, tag, strlen (tag)) != 0)
    return PARSE_NOMATCH;

  if (have != NULL && *have != 0)
    {
      fprintf (stderr, "parse_line: %s appears twice\n", tag);
      return PARSE_ERROR;
    }
  if (have != NULL)
    *have = 1;
  
  lineptr += strlen (tag);
  if (type == TYPE_STRING) /* character string of length up to 256 */
    {
      strncpy ((char *)target, lineptr, 256);
      ((char *)target)[255] = '0';
    }
  else if (type == TYPE_MPZ)
    {
      mpz_init (*(mpz_t *)target);
      mpz_set_str (*(mpz_t *)target, lineptr, 0);
    }
  else if (type == TYPE_INT)
    *(int *)target = atoi (lineptr);
  else if (type == TYPE_ULONG)
    *(unsigned long *)target = strtoul (lineptr, NULL, 10);
  else if (type == TYPE_DOUBLE)
    *(double *)target = atof (lineptr);
  else return PARSE_ERROR;
  
  return PARSE_MATCH;
}

cado_poly *
read_polynomial (char *filename)
{
  FILE *file;
  const int linelen = 512;
  char line[linelen];
  cado_poly *poly;
  int have_name = 0, have_n = 0, have_Y0 = 0, have_Y1 = 0;
  int i, ok = PARSE_ERROR;

  file = fopen (filename, "r");
  if (file == NULL)
    {
      fprintf (stderr, "read_polynomial: could not open %s\n", filename);
      return NULL;
    }

  poly = (cado_poly *) malloc (sizeof(cado_poly));
  (*poly)->f = (mpz_t *) malloc (MAXDEGREE * sizeof (mpz_t));
  (*poly)->g = (mpz_t *) malloc (2 * sizeof (mpz_t));

  (*poly)->name[0] = '\0';
  (*poly)->degree = -1;
  (*poly)->type[0] = '\0';

  while (!feof (file))
    {
      ok = 1;
      if (fgets (line, linelen, file) == NULL)
	break;
      if (line[0] == '#')
	continue;

      ok *= parse_line (&((*poly)->name), line, "name: ", &have_name, TYPE_STRING);
      ok *= parse_line (&((*poly)->n), line, "n: ", &have_n, TYPE_MPZ);
      ok *= parse_line (&((*poly)->skew), line, "skew: ", NULL, TYPE_DOUBLE);
      ok *= parse_line (&((*poly)->g[0]), line, "Y0: ", &have_Y0, TYPE_MPZ);
      ok *= parse_line (&((*poly)->g[1]), line, "Y1: ", &have_Y1, TYPE_MPZ);
      ok *= parse_line (&((*poly)->type), line, "type: ", NULL, TYPE_STRING);
      ok *= parse_line (&((*poly)->rlim), line, "rlim: ", NULL, TYPE_ULONG);
      ok *= parse_line (&((*poly)->alim), line, "alim: ", NULL, TYPE_ULONG);
      ok *= parse_line (&((*poly)->lpbr), line, "lpbr: ", NULL, TYPE_INT);
      ok *= parse_line (&((*poly)->lpba), line, "lpba: ", NULL, TYPE_INT);
      ok *= parse_line (&((*poly)->mfbr), line, "mfbr: ", NULL, TYPE_INT);
      ok *= parse_line (&((*poly)->mfba), line, "mfba: ", NULL, TYPE_INT);
      ok *= parse_line (&((*poly)->rlambda), line, "rlambda: ", NULL, TYPE_DOUBLE);
      ok *= parse_line (&((*poly)->alambda), line, "alambda: ", NULL, TYPE_DOUBLE);
      ok *= parse_line (&((*poly)->qintsize), line, "qintsize: ", NULL, TYPE_INT);


      if (ok == 1 && line[0] == 'c' && isdigit (line[1]) && line[2] == ':' &&
	       line[3] == ' ')
	{
	  int index = line[1] - '0', i;
	  for (i = (*poly)->degree + 1; i <= index; i++)
	    mpz_init ((*poly)->f[i]);
	  if (index > (*poly)->degree)
	    (*poly)->degree = index;
	  mpz_set_str ((*poly)->f[index], line + 4, 0);
	  ok = -1;
	}

      if (ok == PARSE_NOMATCH)
	{
	  fprintf (stderr, 
		   "read_polynomial: Cannot parse line %s\nIgnoring.\n", 
		   line);
	  continue;
	}
      
      if (ok == PARSE_ERROR)
	break;
    }
  
  if (ok != PARSE_ERROR)
    {
      if (have_n == 0)
	{
	  fprintf (stderr, "n ");
	  ok = PARSE_ERROR;
	}
      if (have_Y0 == 0)
	{
	  fprintf (stderr, "Y0 ");
	  ok = PARSE_ERROR;
	}
      if (have_Y1 == 0)
	{
	  fprintf (stderr, "Y1 ");
	  ok = PARSE_ERROR;
	}
      if (ok == PARSE_ERROR)
	fprintf (stderr, "are missing in polynomial file\n");
    }

  if (ok == PARSE_ERROR)
    {
      for (i = 0; i <= (*poly)->degree; i++)
	mpz_clear ((*poly)->f[i]);
      if (have_n)
	mpz_clear ((*poly)->n);
      if (have_Y0)
	mpz_clear ((*poly)->g[0]);
      if (have_Y1)
	mpz_clear ((*poly)->g[1]);
      free ((*poly)->f);
      free ((*poly)->g);
      free (poly);
      poly = NULL;
    }
  
  return poly;
}

void
print_useful (const fbprime_t *useful_primes, 
	      const unsigned int useful_length)
{
  unsigned int useful_count = 0;

  if (useful_primes == NULL)
    return;

  if (*useful_primes == 0)
    {
      printf ("# There were no useful primes\n");
      return;
    }

  printf ("# Useful primes were: ");
  while (*useful_primes != 0)
    {
      useful_count ++;
      printf (FBPRIME_FORMAT " ", *useful_primes++);
    }
  printf ("\n");
  if (useful_count == useful_length - 1)
    printf ("#Storage for useful primes is full, "
	    "consider increasing useful_length\n");
  fflush (stdout);
}

static unsigned long
find_sieve_reports (const unsigned char *sievearray, long *reports, 
                    const unsigned int reports_len, 
                    const unsigned char reports_threshold, 
                    const long amin, const long amax, const unsigned long b,
                    const int odd, const int verbose)
{
  long long tsc1, tsc2;
  long a;
  unsigned long reports_nr, d;
  const int b3 = (b % 3 == 0); /* We skip over $a$, $3|a$ if $3|b$ to save some
				  space in the reports list */
  
  ASSERT (odd == 0 || odd == 1);
  ASSERT (!odd || (amin & 1) == 1);
  ASSERT (!odd || (amax & 1) == 1);
  ASSERT (amin <= amax);
  ASSERT (b > 0);

  rdtscll (tsc1);
  reports_nr = 0;

  for (a = amin, d = 0; reports_nr < reports_len && a <= amax; 
       a += 1 << odd, d++)
    {
      /* The + 10 is to deal with accumulated rounding error that might have
	 cause the sieve value to drop below 0 and wrap around */
      if ((unsigned char) (sievearray[d] + 10) <= reports_threshold + 10)
        {
	  ASSERT (a == amin + (long) (d << odd));
	  if (!b3 || (a % 3) != 0) /* FIXME: have 2 separate loops */
	    reports[reports_nr++] = a;
        }
    }
  rdtscll (tsc2);
#ifdef HAVE_MSRH
  if (verbose)
    {
      printf ("# There were %lu sieve reports (a,b not necessarily coprime)\n",
	      reports_nr);
      printf ("# Finding sieve reports took %lld clocks\n", tsc2 - tsc1);
    }
#endif

  return reports_nr;
}

unsigned long
sieve_one_side (unsigned char *sievearray, factorbase_t fb,
		long *reports, const unsigned int reports_len, 
		const unsigned char reports_threshold,
		fbprime_t *useful_primes, const unsigned int useful_length,
		const long amin, const long amax, const unsigned long b,
		const unsigned long proj_roots, const double log_scale,
		const double *dpoly, const unsigned int deg, 
		const int verbose)
{
  long long tsc1, tsc2;
  unsigned long reports_nr;
  const int odd = 1 - (b & 1); /* If odd=1, only odd $a$ are sieved */
  const long eff_amin = amin + ((odd && (amin & 1) == 0) ? 1 : 0);
  const long eff_amax = amax - ((odd && (amax & 1) == 0) ? 1 : 0);
  const unsigned long l = ((amax - amin) >> odd) + 1;

  fb_disable_roots (fb->fblarge, b, verbose);
  /* FIXME: need to disable entries in L1fb as well */
  
  compute_norms (sievearray, eff_amin, eff_amax, b, dpoly, deg, proj_roots, 
		 log_scale, odd, verbose);

  /* Sieve L1 primes */
  if (fb->fbL1 != NULL)
    {
      unsigned long blockstart;
      
      fb_initloc_small (fb->fbL1init, fb->fbL1, eff_amin, b, odd);
      
      rdtscll (tsc1);
      for (blockstart = 0; blockstart < l; blockstart += fb->fbL1bound)
	{
	  if (fb->fbL1bound <= l - blockstart)
	    sieve_small_slow (sievearray + blockstart, fb->fbL1init, 
			      fb->fbL1bound);
	  else
	    sieve_small_slow (sievearray + blockstart, fb->fbL1init, 
			      l - blockstart);
	}
      rdtscll (tsc2);
#ifdef HAVE_MSRH
      if (verbose)
	printf ("# Sieving L1 primes took %lld clocks\n", tsc2 - tsc1);
#endif
    }

  rdtscll (tsc1);
  sieve (sievearray, fb->fblarge, eff_amin, eff_amax, b, reports_threshold, 
	 useful_primes, useful_length, odd);
  rdtscll (tsc2);
#ifdef HAVE_MSRH
  if (verbose)
    printf ("# Sieving large fb primes took %lld clocks\n", tsc2 - tsc1);
#endif
  
  fb_restore_roots (fb->fblarge, b, verbose);
  
  if (verbose >= 2)
    print_useful (useful_primes, useful_length);
  
  /* Store the sieve reports */
  reports_nr = find_sieve_reports (sievearray, reports, reports_len, 
				   reports_threshold, eff_amin, eff_amax, b, 
				   odd, verbose);

  return reports_nr;
}

static inline void
add_fbprime_to_list (fbprime_t *list, unsigned int *cursize, 
		     const unsigned int maxsize, const fbprime_t toadd)
{
  if ((*cursize) < maxsize)
    list[(*cursize)++] = toadd;
}

/* Divides the prime q and its powers out of C, appends each q divided out 
   to primes_a, subject to the length restriction *nr_primes < max_nr_primes. 
   Returns the exponent of q that divided (0 if q didn't divide at all). */

static inline int
trialdiv_one_prime (const fbprime_t q, mpz_t C, unsigned int *nr_primes,
                    fbprime_t *primes, const unsigned int max_nr_primes)
{
  int nr_divide = 0;

  while (mpz_divisible_ui_p (C, q))
    {
      nr_divide++;
      add_fbprime_to_list (primes, nr_primes, max_nr_primes, q);
      mpz_tdiv_q_ui(C, C, q);
    }
  return nr_divide;
}

/* Factor an unsigned long into fb primes, adds primes to list. Not fast. */
static inline void 
trialdiv_slow (unsigned long C, fbprime_t *primes, unsigned int *nr_primes, 
	       const unsigned int max_nr_primes)
{
  fbprime_t q;

  while (C > 1)
    {
      q = iscomposite (C);
      if (q == 0)
	{
	  q = (fbprime_t) C;
	  ASSERT ((unsigned long) q == C); /* Check for truncation */
	}
      C /= (unsigned long) q;
      add_fbprime_to_list (primes, nr_primes, max_nr_primes, q);
    }
}

static int
trialdiv_one_side (mpz_t norm, mpz_t *scaled_poly, const long a, 
		   const unsigned long b, int degree, 
		   fbprime_t *primes, unsigned int *nr_primes, 
		   const unsigned int max_nr_primes,
		   const unsigned long proj_divisor, 
		   const unsigned int nr_proj_primes, 
		   const fbprime_t *proj_primes,
		   const fbprime_t *useful_primes,
		   factorbase_degn_t *fullfb,
		   unsigned int *cof_toolarge,
		   unsigned int *lp_toolarge,
		   const int lpb, const int mfb, const double lambda)
{
  unsigned int k;
  fbprime_t largest_useful;
  factorbase_degn_t *fbptr;
  size_t log2size;

  /* 1. Compute norm */
  mp_poly_eval (norm, scaled_poly, degree, a);
  mpz_abs (norm, norm);
  
  *nr_primes = 0;
  
  /* 2. Divide out the primes with projective roots */
  if (proj_divisor > 1)
    {
      unsigned long r;
      ASSERT_ALWAYS (nr_proj_primes > 0);
      r = mpz_tdiv_q_ui (norm, norm, proj_divisor);
      ASSERT_ALWAYS (r == 0);
      
      for (k = 0; k < nr_proj_primes; k++)
	add_fbprime_to_list (primes, nr_primes, max_nr_primes, proj_primes[k]);
    }
  
  /* 3. See if any of the "useful primes" divide this norm */
  largest_useful = 0;
  for (k = 0; useful_primes[k] > 0; k++)
    {
      if (trialdiv_one_prime (useful_primes[k], norm, nr_primes, 
			      primes, max_nr_primes))
	{
	  if (largest_useful == 0 || largest_useful < useful_primes[k])
	    largest_useful = useful_primes[k];
	}
    }
  if (0 && largest_useful == 0)
    {
      fprintf (stderr, "Warning: norm of degree %d polynomial for report "
	       "%lu, %lu had no useful prime\n", degree, a, b);
      gmp_fprintf (stderr, "Norm is %Zd\n", norm);
    }
  
  /* 4. Go through the factor base until the fb prime is equal to the 
     largest "useful" prime we've seen (all larger fb primes that
     divide the norm must also have been useful), or until the 
     end of the factor base. */ 
  for (fbptr = fullfb; 
       fbptr->p != largest_useful; /* This is a trick. If there was a  useful prime, we'll hit that eventually and stop there. If there wasn't, largest_useful equals 0, which is also the end-of-factorbase-marker, so we'll stop at the end. */
       fbptr = fb_next (fbptr))
    trialdiv_one_prime (fbptr->p, norm, nr_primes, primes,
			max_nr_primes);
  
  /* 5. Check if the cofactor is small enough */
  log2size = mpz_sizeinbase (norm, 2) - 1;
  if ((double) log2size > 
      lpb * lambda + SIEVE_PERMISSIBLE_ERROR)
    {
      gmp_fprintf (stderr, 
		   "Sieve report %ld, %lu is not smooth for degree %d poly, "
		   "cofactor is %Zd with %d bits\n", 
		   a, b, degree, norm, mpz_sizeinbase (norm, 2));
      return 0;
    }
  if ((int) log2size > mfb)
    {
      (*cof_toolarge)++;
      return 0;
    }
  
  /* 6. Check if the cofactor is < lpb */
  if ((int) log2size <= lpb)
    {
      fbprime_t q;

      /* We assume that this implies that the cofactor is 1 or prime,
	 i.e. that 2^lbp < rlima^2 */

      if (mpz_cmp_ui (norm, 1UL) > 0)
	{
	  ASSERT (mpz_probab_prime_p (norm, PRP_REPS));
	  q = (fbprime_t) mpz_get_ui (norm);
	  add_fbprime_to_list (primes, nr_primes, max_nr_primes, q);
	  mpz_set_ui (norm, 1UL);
	}
      return 1;
    }
  else
    {
      /* If this cofactor is a prp, since it's > lpb, 
	 skip this report */
      if (mpz_probab_prime_p (norm, PRP_REPS))
	{
	  (*lp_toolarge)++;
	  return 0;
	}
    }
  return 1;
}


static void
trialdiv_and_print (cado_poly *poly, const unsigned long b, 
                    const long *reports_a, const unsigned int reports_a_nr, 
		    fbprime_t *useful_primes_a, const long *reports_r, 
		    const unsigned int reports_r_nr, 
		    fbprime_t *useful_primes_r, 
		    factorbase_t fba, factorbase_t fbr, const int verbose)
{
  unsigned long long tsc1, tsc2;
  unsigned long proj_divisor_a, proj_divisor_r;
  unsigned int i, j, k;
  const unsigned int max_nr_primes = 128, max_nr_proj_primes = 16;
  mpz_t Fab, Gab, scaled_poly_a[MAXDEGREE], scaled_poly_r[2];
  fbprime_t primes_a[max_nr_primes], primes_r[max_nr_primes];
  fbprime_t proj_primes_a[max_nr_primes], proj_primes_r[max_nr_primes];
  unsigned int nr_primes_a, nr_primes_r, nr_proj_primes_a, nr_proj_primes_r;
  unsigned int lp_a_toolarge = 0, lp_r_toolarge = 0, 
    cof_a_toolarge = 0, cof_r_toolarge = 0;
  int ok;

  mpz_init (Fab);
  mpz_init (Gab);
  for (i = 0; i <= (unsigned) (*poly)->degree; i++)
    mpz_init (scaled_poly_a[i]);
  mpz_init (scaled_poly_r[0]);
  mpz_init (scaled_poly_r[1]);

  rdtscll (tsc1);

  /* Multiply f_i (and g_i resp.) by b^(deg-i) and put in scaled_poly */
  mp_poly_scale (scaled_poly_a, (*poly)->f, (*poly)->degree, b, -1); 
  mp_poly_scale (scaled_poly_r, (*poly)->g, 1, b, -1); 

  /* Factor primes with projective roots on algebraic side. Prime factors
     go in the proj_primes_a[] list, the number of prime factors in 
     nr_proj_primes_a, the product is kept in proj_divisor_a. */
  nr_proj_primes_a = 0;
  proj_divisor_a = mpz_gcd_ui (NULL, (*poly)->f[(*poly)->degree], b);
  trialdiv_slow (proj_divisor_a, proj_primes_a, &nr_proj_primes_a, 
		 max_nr_proj_primes);

  /* Same for rational side */
  nr_proj_primes_r = 0;
  proj_divisor_r = mpz_gcd_ui (NULL, (*poly)->g[1], b);
  trialdiv_slow (proj_divisor_r, proj_primes_r, &nr_proj_primes_r, 
		 max_nr_proj_primes);

  for (i = 0, j = 0; i < reports_a_nr && j < reports_r_nr;)
    {
      if (reports_a[i] == reports_r[j] && gcd(labs(reports_a[i]), b) == 1)
	{
          const long a = reports_a[i];
      
	  /* For the algebraic side, then for the rational side, we
	     1. Compute the norm
	     2. Divide out primes with projective roots
	     3. Divide out the useful prime(s), remembering the smallest one
	     4. Trial divide up to the smallest useful prime (or up to the
	        factor base limit if there was no useful prime)
	     5. Check if the cofactor is small enough
	     6. Check if the cofactor is < lpb, if yes report relation
	     7. Check if the cofactor is a prp, if yes skip to next report

	     Factoring composite cofactors into large primes is left to the 
	     next step in the tool chain.
	  */

          /* Do the algebraic side */

	  ok = trialdiv_one_side (Fab, scaled_poly_a, a, b, (*poly)->degree, 
				  primes_a, &nr_primes_a, max_nr_primes,
				  proj_divisor_a, nr_proj_primes_a, 
				  proj_primes_a, useful_primes_a, fba->fullfb, 
				  &cof_a_toolarge, &lp_a_toolarge, 
				  (*poly)->lpba, (*poly)->mfba, 
				  (*poly)->alambda);

	  if (!ok) 
	    goto nextreport;

          /* Now the rational side */
	  
	  ok = trialdiv_one_side (Gab, scaled_poly_r, a, b, 1, 
				  primes_r, &nr_primes_r, max_nr_primes,
				  proj_divisor_r, nr_proj_primes_r, 
				  proj_primes_r, useful_primes_r, fbr->fullfb, 
				  &cof_r_toolarge, &lp_r_toolarge, 
				  (*poly)->lpbr, (*poly)->mfbr, 
				  (*poly)->rlambda);

	  if (!ok) 
	    goto nextreport;

	  /* Now print the relations */
	  printf ("%ld,%lu:", a, b);
	  for (k = 0; k < nr_primes_r; k++)
	    printf ("%x%s", primes_r[k], k+1==nr_primes_r?":":",");
	  for (k = 0; k < nr_primes_a; k++)
	    printf ("%x%s", primes_a[k], k+1==nr_primes_a?"\n":",");
	  fflush (stdout);
	}


    nextreport:
      /* Assumes values in reports_a are sorted and unique, 
	 same for reports_r  */
      if (reports_a[i] < reports_r[j])
	i++;
      else
	j++;
    }
  
  for (i = 0; i <= (unsigned)(*poly)->degree; i++)
    mpz_clear (scaled_poly_a[i]);
  mpz_clear (scaled_poly_r[0]);
  mpz_clear (scaled_poly_r[1]);
      
  mpz_clear (Fab);
  mpz_clear (Gab);
  
  rdtscll (tsc2);
#ifdef HAVE_MSRH
  if (verbose)
  {
    printf ("# Trial factoring/printing took %lld clocks\n", tsc2 - tsc1);
    printf ("# Too large cofactors (discarded in this order): "
	    "alg composite > mfba: %d, alg prp > lpba: %d, "
	    "rat composite > mfbr: %d, rat prp > lpbr: %d\n", 
	    cof_a_toolarge, lp_a_toolarge, cof_r_toolarge, lp_r_toolarge);
  }
#endif
}


int
main (int argc, char **argv)
{
  long amin, amax;
  unsigned long bmin, bmax, b;
  fbprime_t *useful_primes_a, *useful_primes_r;
  factorbase_t fba, fbr;
  long *reports_a, *reports_r;
  char *fbfilename = NULL, *polyfilename = NULL;
  unsigned char *sievearray;
  unsigned int useful_length_a, useful_length_r, reports_a_len, reports_a_nr, 
    reports_r_len, reports_r_nr;
  int verbose = 0;
  unsigned int deg;
  unsigned int i;
  double dpoly_a[MAXDEGREE], dpoly_r[2];
  const double log_scale = 1. / log (2.); /* Lets use log_2() for a start */
  cado_poly *cpoly;
  char report_a_threshold, report_r_threshold;


  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 1 && strcmp (argv[1], "-v") == 0)
	{
	  verbose++;
	  argc--;
	  argv++;
	}
      else if (argc > 2 && strcmp (argv[1], "-fb") == 0)
	{
	  fbfilename = argv[2];
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-poly") == 0)
	{
	  polyfilename = argv[2];
	  argc -= 2;
	  argv += 2;
	}

      else 
	break;
    }

  if (argc < 5)
    {
      fprintf (stderr, "Please specify amin amax bmin bmax\n");
      exit (EXIT_FAILURE);
    }

  amin = atol (argv[1]);
  amax = atol (argv[2]);
  bmin = strtoul (argv[3], NULL, 10);
  bmax = strtoul (argv[4], NULL, 10);

  /* Check command line parameters for sanity */

  if (amin >= amax)
    {
      fprintf (stderr, "amin must be less than amax\n");
      exit (EXIT_FAILURE);
    }

  if (bmin > bmax)
    {
      fprintf (stderr, "bmin must be less than or equal to bmax\n");
      exit (EXIT_FAILURE);
    }

  if (polyfilename == NULL)
    {
      fprintf (stderr, "Please specify a polynomial file with -poly\n");
      exit (EXIT_FAILURE);
    }

  if (fbfilename == NULL)
    {
      fprintf (stderr, 
	       "Please specify a factor base file with the -fb option\n");
      exit (EXIT_FAILURE);
    }

  cpoly = read_polynomial (polyfilename);
  if (cpoly == NULL)
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }
  if (verbose)
  {
      printf ("Read polynomial file %s\n", polyfilename);
      printf ("Polynomials are:\n");
      mp_poly_print ((*cpoly)->f, (*cpoly)->degree, "f(x) =", 0);
      printf ("\n");
      mp_poly_print ((*cpoly)->g, 1, "g(x) =", 0);
      printf ("\n");
  }

#ifdef PARI
  mp_poly_print ((*cpoly)->f, (*cpoly)->degree, "f(x) =");
  printf (" /* PARI */\n");
  printf ("F(a,b) = f(a/b)*b^%d /* PARI */\n", (*cpoly)->degree);
  mp_poly_print ((*cpoly)->g, 1, "g(x) =");
  printf (" /* PARI */\n");
  printf ("G(a,b) = g(a/b)*b /* PARI */\n");
#endif


  /* Read/generate and split the factor bases */

  /* Read algebraic fb */
  fba->fullfb = fb_read (fbfilename, log_scale, verbose);
  if (fba->fullfb == NULL)
    {
      fprintf (stderr, "Could not read factor base\n");
      exit (EXIT_FAILURE);
    }
  fba->fblarge = fba->fullfb;
  /* Extract the primes < L1SIZE into their small prime array */
  if (1)
    fb_extract_small (fba, L1SIZE, verbose);
  else
    fba->fbL1 = NULL;

  /* Generate rational fb */
  fbr->fullfb = fb_make_linear ((*cpoly)->g, (fbprime_t) (*cpoly)->rlim, 
				log_scale, verbose);
  if (fbr == NULL)
    {
      fprintf (stderr, 
	       "Could not generate factor base for linear polynomial\n");
      exit (EXIT_FAILURE);
    }
  fbr->fblarge = fbr->fullfb;
  if (1)
    fb_extract_small (fbr, L1SIZE, verbose);
  else
    fbr->fbL1 = NULL;
  
  deg = (*cpoly)->degree;
  for (i = 0; i <= deg; i++)
      dpoly_a[i] = mpz_get_d ((*cpoly)->f[i]);
  dpoly_r[0] = mpz_get_d ((*cpoly)->g[0]);
  dpoly_r[1] = mpz_get_d ((*cpoly)->g[1]);
  report_a_threshold = (unsigned char) ((double)((*cpoly)->lpba) * log(2.0) * 
					(*cpoly)->alambda * log_scale + 0.5);
  report_r_threshold = (unsigned char) ((double)((*cpoly)->lpbr) * log(2.0) * 
					(*cpoly)->rlambda * log_scale + 0.5);

#if 0
  if (verbose)
    for (s = 0; fb_skip(fb, s)->p != 0; s += fb_entrysize (fb_skip (fb, s)))
      fb_print_entry (fb_skip(fb, s));
#endif

  sievearray = (unsigned char *) malloc ((amax - amin + 1) * sizeof (char));
  useful_length_a = ((amax - amin + 1)) / 100 + 1000;
  useful_length_r = ((amax - amin + 1)) / 10 + 1000;
  useful_primes_a = (fbprime_t *) malloc (useful_length_a* sizeof (fbprime_t));
  useful_primes_r = (fbprime_t *) malloc (useful_length_r* sizeof (fbprime_t));
  ASSERT (useful_primes_a != NULL && useful_primes_r != NULL);

  reports_a_len = useful_length_a;
  reports_a = (long *) malloc (reports_a_len * sizeof (long));
  ASSERT (reports_a != NULL);
  reports_r_len = useful_length_r;
  reports_r = (long *) malloc (reports_r_len * sizeof (long));
  ASSERT (reports_r != NULL);

  for (b = bmin; b <= bmax; b++)
    {
      unsigned long proj_roots; /* = gcd (b, leading coefficient) */
      
      if (verbose)
	printf ("# Sieving line b = %lu\n", b);

      proj_roots = mpz_gcd_ui (NULL, (*cpoly)->f[deg], b);
      if (verbose)
	printf ("# Projective roots for b = %lu on algebtaic side are: %lu\n", 
	        b, proj_roots);

      if (verbose)
	  printf ("# Sieving algebraic side\n");

#ifdef PARI
      mp_poly_print ((*cpoly)->f, (*cpoly)->degree, "P(a,b) = ", 1);
      printf (" /* PARI */\n");
#endif

      reports_a_nr = 
	sieve_one_side (sievearray, fba,  
			reports_a, reports_a_len, report_a_threshold, 
			useful_primes_a, useful_length_a, 
			amin, amax, b, proj_roots, log_scale, 
			dpoly_a, deg, verbose);

      proj_roots = mpz_gcd_ui (NULL, (*cpoly)->g[1], b);
      if (verbose)
	printf ("# Projective roots for b = %lu on rational side are: %lu\n", 
	        b, proj_roots);

      if (verbose)
	  printf ("#Sieving rational side\n");

#ifdef PARI
      mp_poly_print ((*cpoly)->g, 1, "P(a,b) = ", 1);
      printf (" /* PARI */\n");
#endif

      reports_r_nr = 
	sieve_one_side (sievearray, fbr, 
			reports_r, reports_r_len, 
			report_r_threshold, useful_primes_r, 
			useful_length_r, amin, amax, b, proj_roots, log_scale, 
			dpoly_r, 1, verbose);

      if (reports_a_len == reports_a_nr)
	  fprintf (stderr, "Warning: sieve reports list on algebraic side "
		   "full with %u entries for b=%lu\n", reports_a_len, b);

      if (reports_r_len == reports_r_nr)
	  fprintf (stderr, "Warning: sieve reports list on rational side "
		   "full with %u entries for b=%lu\n", reports_r_len, b);

      /* Not trial factor the candidate relations */
      trialdiv_and_print (cpoly, b, 
			  reports_a, reports_a_nr, useful_primes_a, 
			  reports_r, reports_r_nr, useful_primes_r, 
			  fba, fbr, verbose);
    }

#if 0
  /* FIXME implement this */
  fb_clear (fba);
  fb_clear (fbr);
#endif
  free (useful_primes_a);
  free (useful_primes_r);
  free (reports_a);
  free (reports_r);
  free (sievearray);

  return 0;
}
