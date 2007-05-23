/* Apply the sieve updates to sievearray */
/* Output: sievearray. Must have enough memory allocated */
/* Inputs: amin, amax, b the $a$ range and the line $b$ to sieve */
/* factorbase, polynomial */

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
#include <assert.h>
#define ASSERT_ALWAYS(x) assert(x)
#ifdef WANT_ASSERT
  #define ASSERT(x) assert(x)
#else
  #define ASSERT(x)
#endif

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define REFAC_PRP_SIZE_THRES 50
#define REFAC_PRP_THRES 500
#define SIEVE_PERMISSIBLE_ERROR 5

/*****************************************************************
 *                      Functions for calculus                   *
 *****************************************************************/

/* Evaluate poly */

double 
calc_poly_eval (const double *poly, const int deg, const double x)
{
    double r;
    int i;

    ASSERT (deg <= MAXDEGREE);
    r = poly[deg];
    for (i = deg - 1; i >= 0; i--)
	r = r * x + poly[i];

    return r;
}

/* deriv == poly is ok */
void
calc_deriv (double *deriv, const double *poly, const int deg)
{
    int i;
    ASSERT (deg <= MAXDEGREE);
    for (i = 1; i <= deg; i++)
	deriv[i - 1] = (double) i * poly[i];
}

/* Newton root finding */
double
calc_newton_root (const double *poly, const int deg, const double start)
{
    double x, lastx, deriv[MAXDEGREE + 1];

    ASSERT (deg <= MAXDEGREE);
    calc_deriv (deriv, poly, deg);

    x = start;
    lastx = x + 1.;
    while (fabs (x - lastx) > 1.e-5) /* 1.e-5 should not take too long,
					even for repeated roots */
    {
	lastx = x;
	x -= calc_poly_eval(poly, deg, x) / calc_poly_eval(deriv, deg - 1, x);
    }

    /* Two more, with quadratic convergence that should give error < 1.e-20
       which is good enough even for long double. For repeated roots: FIXME.*/
    x -= calc_poly_eval(poly, deg, x) / calc_poly_eval(deriv, deg - 1, x);
    x -= calc_poly_eval(poly, deg, x) / calc_poly_eval(deriv, deg - 1, x);

    return x;
}

/* Divide a polynomial by (x - root). Returns the remainder, which should
   be small if root was indeed a root */
double
calc_divide_root (double *quotient, const double *poly, const int deg, 
		  const double root)
{
    double a, b;
    int i;

    a = poly[deg];
    for (i = deg; i > 0; i--)
    {
	b = a;
	a = poly[i - 1] + a * root;
	quotient[i - 1] = b;
    }
    
    return a;
}

int
calc_find_extrema (double *extrema, const double *poly, const int deg)
{
    double deriv[MAXDEGREE + 1];
    double root;
    double *derivptr = deriv;
    int i, j = 0;
    
    ASSERT(deg <= MAXDEGREE);
    ASSERT(poly[deg] != 0.);
    /* Constant polys have no extrema and would mess up code below */
    if (deg == 0)
	return 0;
    calc_deriv (deriv, poly, deg);
    /* Deal with multiple roots at 0 which may be plenty in SNFS polys */
    if (*derivptr == 0.)
    {
	extrema[j++] = 0.;
	while (*derivptr == 0.)
	    derivptr++; 
    }
    /* Find number of real roots. TBD. */
    /* FIXME for repeated roots */
    for (i = deg - 1; i > 0; i--)
    {
	root = calc_newton_root (deriv, i, 0.);
	extrema[j++] = root;
	printf ("calc_find_extrema: f´(%f) = %f\n", 
		root, calc_poly_eval(deriv, i, root));
	calc_divide_root (deriv, deriv, i, root);
    }
    return j;
}



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
mp_poly_print (mpz_t *poly, int deg, const char *name)
{
  int i;
  printf ("%s", name);
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

/*****************************************************************
 *           Simple functions for (modular) arithmetic           *
 *****************************************************************/


/* Returns 0 if n is prime, otherwise the smallest prime factor of n */
fbprime_t
iscomposite (const fbprime_t n)
{
  fbprime_t i, i2;

  if (n % 2 == 0)
    return (n == 2) ? 0 : 2;

  /* (i + 2)^2 = i^2 + 4*i + 4 */
  for (i = 3, i2 = 9; i2 <= n; i2 += (i+1) * 4, i += 2)
    if (n % i == 0)
	return i;

  return 0;
}

static unsigned long
gcd (unsigned long a, unsigned long b)
{
  unsigned long t;
  uint32_t a32, b32, t32;

  ASSERT (b > 0);
  
  if (a >= b)
    a %= b;

  while (a > UINT32_MAX)
    {
      if (b - a < a)
	t = b - a;
      else
	t = b % a;
      b = a;
      a = t;
    }

  if (a == 0)
    return b;

  a32 = a;
  b32 = b;
  t32 = t;

  while (a32 > 0)
    {
      t32 = b32 % a32; /* The 32 bit DIV is much faster */
      b32 = a32;
      a32 = t32;
    }

  return b32;
}

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
  const long eff_amin = amin + ((odd && amin % 2 == 0) ? 1L : 0L);
  const long eff_amax = amax - ((odd && amax % 2 == 0) ? 1L : 0L);
  long a, a2;
  int i;
  unsigned char n1, n2, nmax;
  const int stride = 128;
  
  ASSERT_ALWAYS (odd == 0 || odd == 1);
  ASSERT_ALWAYS (amin < amax);

  rdtscll (tsc1);
  
  bpow = 1.;
  for (i = deg; i >= 0; i--)
    {
      f[i] = poly[i] * bpow;
      bpow *= (double) b;
    }
  
  n1 = log_norm (f, deg, (double) eff_amin, log_scale, log_proj_roots);
  nmax = n1;

  a = eff_amin;
  a2 = eff_amin;
  /* a == eff_amin + (d << odd), d = (a - eff_amin) >> odd */

  while (a <= eff_amax)
    {
      /* We'll cover [a ... a2[ now */
      ASSERT (a == a2);
      ASSERT (n1 == log_norm (f, deg, (double) a, log_scale, log_proj_roots));
      a2 = MIN(a + stride, eff_amax + (1 << odd));
      ASSERT (!odd || (a2 - a) % 2 == 0);
      
      n2 = log_norm (f, deg, (double) a2, log_scale, log_proj_roots);
      
      if (n1 == n2) 
	{
	  /* Let's assume the log norm is n1 everywhere in this interval */
	  memset (sievearray + ((a - eff_amin) >> odd), n1, (a2 - a) >> odd);
	}
      else
	{
	  /* n1 and n2 are different. Do each a up to (exclusive) a2 
	     individually */
	  /* printf ("log_c(F(%ld, %lu)) == %d != log_c(F(%ld, %lu)) == %d\n",
	     a, b, n1, a2, b, n2); */
	  sievearray[(a - eff_amin) >> odd] = n1;
	  a += 1 << odd;
	  for ( ; a < a2; a += 1 << odd)
	    {
	      unsigned char n = log_norm (f, deg, (double) a, log_scale, 
					  log_proj_roots);
	      sievearray[(a - eff_amin) >> odd] = n;
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
  for (a = eff_amin; a <= eff_amax; a += 1 << odd)
    printf ("if(log_c(%c(%ld, %lu), %f) != %d,"
	    "print (\"log_c(%c(\", %ld,\", \", %lu\"), %f) \", %d)) /* PARI */\n", 
	    deg > 1 ? 'F' : 'G', a, b, log_proj_roots, sievearray[(a - eff_amin) >> odd],
	    deg > 1 ? 'F' : 'G', a, b, log_proj_roots, sievearray[(a - eff_amin) >> odd]); 
#endif

  return nmax;
}

/* sievearray must be long enough to hold amax-amin+1 chars.
   proj_roots is the rounded log of the projective roots */

/* odd: if 1, the sieve array contains only locations for odd a, 
        amin <= a <= amax */

void 
sieve (unsigned char *sievearray, factorbase_t *fb, 
       const long amin, const long amax, const unsigned long b, 
       const unsigned char threshold, fbprime_t *useful_primes_par,
       const fbprime_t useful_threshold, const unsigned int useful_length_par,
       const int odd)
{
  /* The sievearray[0] entry corresponds to (eff_amin, b), and
     sievearray[d] to (eff_amin + d * (1 + odd), b) */
  const long eff_amin = amin + ((odd && amin % 2 == 0) ? 1L : 0L);
  const long eff_amax = amax - ((odd && amax % 2 == 0) ? 1L : 0L);
  fbprime_t *useful_primes = useful_primes_par;
  const uint32_t l = eff_amax - eff_amin + 1;
  uint32_t i, amin_p, p, d;
  unsigned int useful_length = useful_length_par;
  unsigned char plog;

  ASSERT_ALWAYS (amax > amin);

  if (useful_length > 0)
    useful_length--; /* Make sure we have space for a 0 mark at the end */

  /* Do the sieving */

  while (fb->p > 0)
    {
      ASSERT (fb_entrysize (fb) <= fb->size);
      p = fb->p;
      plog = fb->plog;

      /* Compute a % p for the a value of the sieve location in 
	 sievearray[0] */
      /* FIXME This modular reduction should be simplified somehow. Do it 
	 once and store it in fb? */
      if (eff_amin < 0)
	{
	  amin_p = p - ((unsigned long)(-eff_amin) % p);
	  if (amin_p == p)
	    amin_p = 0; /* FIXME ugly */
	}
      else
        amin_p = (unsigned long) eff_amin % p;

      for (i = 0; i < fb->nr_roots; i++)
        {
	  modulus m;
	  residue r1, r2;

          /* Find first index in sievearray where p on this root divides.
             So we want a/b == r (mod p) <=> a == br (mod p). Then we want
             d so that eff_amin + d * (1 + odd) == a (mod p) 
	     <=> d == (a - eff_amin) / (1 + odd) (mod p) */
	  ASSERT (fb->roots[i] < p);
	  ASSERT (b % p != 0);
	  mod_initmod_ul (m, p); /* Most of the mod_*() calls are no-ops */
	  mod_init (r1, m);
	  mod_init (r2, m);
	  mod_set_ul (r1, b, m); /* Modular reduction */
	  mod_set_ul_reduced (r2, fb->roots[i], m);
          mod_mul (r1, r1, r2, m); /* Multiply and mod reduction. If we keep 
	    r_i in Montgomery representation, a single mul/REDC will compute 
	    the product, reduce it mod p and return it as an integer. TBD. */
	  mod_sub_ul (r1, r1, amin_p, m);
	  if (odd)
	    mod_div2 (r1, r1, m);
          d = mod_get_ul (r1, m); 
          ASSERT (d < p);
	  mod_clear (r1, m);
	  mod_clear (r2, m);
	  mod_clearmod (m);
	  
	  /* Now d is the first index into sievearray where p divides. */

          /* Now update the sieve array. There is no partitioning, blocking, 
             bucket sieving or anything atm. */
          for (; d < l; d += p)
	    {
	      unsigned char k;
	      k = sievearray[d] - plog;
	      sievearray[d] = k;
	      if (k <= threshold)
		if (p >= useful_threshold && useful_length > 0)
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
  const long eff_amin = amin + ((odd && amin % 2 == 0) ? 1 : 0);
  const long eff_amax = amax + ((odd && amax % 2 == 0) ? 1 : 0);
  long a;
  unsigned long reports_nr, odd_b = b, d;
  
  ASSERT (odd == 0 || odd == 1);

  rdtscll (tsc1);
  reports_nr = 0;
  if (odd) /* If all $a$ are odd, we can divide 2's out of $b$ for the gcd */
    for (odd_b = b; odd_b % 2 == 0; odd_b >>= 1);

  for (a = eff_amin, d = 0; reports_nr < reports_len && a <= eff_amax; 
       a += 1 << odd, d++)
    {
      /* The + 10 is to deal with accumulated rounding that might have
	 cause the sieve value to drop below 0 and wrap around */
      if ((unsigned char) (sievearray[d] + 10) <= reports_threshold + 10)
        {
	  ASSERT (a = eff_amin + (d << odd));
          if (gcd(labs(a), odd_b) == 1)
	    reports[reports_nr++] = a;
        }
    }
  rdtscll (tsc2);
#ifdef HAVE_MSRH
  if (verbose)
    {
      printf ("# There were %lu sieve reports\n", reports_nr);
      printf ("# Finding sieve reports took %lld clocks\n", tsc2 - tsc1);
    }
#endif

  return reports_nr;
}

unsigned long
sieve_one_side (unsigned char *sievearray, factorbase_t *fb,
		long *reports, const unsigned int reports_len, 
		const unsigned char reports_threshold,
		fbprime_t *useful_primes, const fbprime_t useful_threshold, 
		const unsigned int useful_length,
		const long amin, const long amax, const unsigned long b,
		const unsigned long proj_roots, const double log_scale,
		const double *dpoly, const unsigned int deg, 
		const int verbose)
{
  long long tsc1, tsc2;
  unsigned long reports_nr;
  int odd = (b % 2 == 0) ? 1 : 0; /* If odd = 1, only odd $a$ are sieved */

  fb_disable_roots (fb, b, verbose);
  
  compute_norms (sievearray, amin, amax, b, dpoly, deg, proj_roots, 
		 log_scale, odd, verbose);
  
  rdtscll (tsc1);
  sieve (sievearray, fb, amin, amax, b, reports_threshold, useful_primes, 
	 useful_threshold, useful_length, odd);
  rdtscll (tsc2);
#ifdef HAVE_MSRH
  if (verbose)
    printf ("# Sieving took %lld clocks\n", tsc2 - tsc1);
#endif
  
  fb_restore_roots (fb, b, verbose);
  
  if (verbose >= 2)
    print_useful (useful_primes, useful_length);
  
  /* Store the sieve reports */
  reports_nr = find_sieve_reports (sievearray, reports, reports_len, 
				   reports_threshold, amin, amax, b, odd, 
				   verbose);

  return reports_nr;
}


/* Divides the prime q and its powers out of C, appends each q divided out
   to primes_a. If anything was divided out, returns 1 if cofactor is 1 or
   a probable prime. Otherwise returns 0. */

static int
trialdiv_one_prime (const fbprime_t q, mpz_t C, unsigned int *nr_primes_a,
                    fbprime_t *primes_a, const unsigned int max_nr_primes)
{
  int is_prp = 0;

  if (mpz_divisible_ui_p (C, q))
    {
      do
	{
	  if (*nr_primes_a < max_nr_primes)
	    primes_a[(*nr_primes_a)++] = q;
	  mpz_tdiv_q_ui(C, C, q);
	} while (mpz_divisible_ui_p (C, q));

      if (mpz_cmp_ui (C, 1) == 0 || 
	  (
#if 1
	   q > REFAC_PRP_THRES
#else
	   mpz_sizeinbase (C, 2) <= REFAC_PRP_SIZE_THRES
#endif
	   && mpz_probab_prime_p (C, 1)))
        {
          is_prp = 1;
        }
    }
  return is_prp;
}

static void
trialdiv_and_print (cado_poly *poly, const unsigned long b, 
                    const long *reports_a, const unsigned int reports_a_nr, 
		    fbprime_t *useful_primes_a, const long *reports_r, 
		    const unsigned int reports_r_nr, 
		    fbprime_t *useful_primes_r, 
		    factorbase_t *fba, const int verbose)
{
  unsigned long long tsc1, tsc2;
  unsigned long proj_primes;
  unsigned int i, j, k;
  const unsigned int max_nr_primes = 128;
  mpz_t Fab, Gab, scaled_poly_a[MAXDEGREE], scaled_poly_r[2];
  fbprime_t primes_a[max_nr_primes], primes_r[max_nr_primes];
  fbprime_t q, incrq;
  unsigned int nr_primes_a, nr_primes_r;
  int cof_a_prp, cof_r_prp;
  unsigned int lp_a_toolarge = 0, lp_r_toolarge = 0, 
    cof_a_toolarge = 0, cof_r_toolarge = 0;
  int log2size;

  mpz_init (Fab);
  mpz_init (Gab);
  for (i = 0; i <= (unsigned) (*poly)->degree; i++)
    mpz_init (scaled_poly_a[i]);
  mpz_init (scaled_poly_r[0]);
  mpz_init (scaled_poly_r[1]);

  /* Multiply f_i by b^(deg-i) and put in scaled_poly */
  mp_poly_scale (scaled_poly_a, (*poly)->f, (*poly)->degree, b, -1); 
  mp_poly_scale (scaled_poly_r, (*poly)->g, 1, b, -1); 

  proj_primes = mpz_gcd_ui (NULL, (*poly)->f[(*poly)->degree], b);

  rdtscll (tsc1);
  for (i = 0, j = 0; i < reports_a_nr && j < reports_r_nr;)
    {
      ASSERT (gcd(labs(reports_a[i]), b) == 1);
      ASSERT (gcd(labs(reports_r[j]), b) == 1);

      if (reports_a[i] == reports_r[j])
	{
          const long a = reports_a[i];
          factorbase_t *nextfb;
      
	  mp_poly_eval (Fab, scaled_poly_a, (*poly)->degree, a);
	  mpz_abs (Fab, Fab);
	  mp_poly_eval (Gab, scaled_poly_r, 1, a);
	  mpz_abs (Gab, Gab);
	  
          nr_primes_a = nr_primes_r = 0;
	  cof_a_prp = cof_r_prp = 0;

          /* Do the algebraic side */

	  /* Divide out the primes with projective roots */
	  if (proj_primes > 1)
	    {
	      ASSERT (mpz_tdiv_ui (Fab, proj_primes) == 0);
	      mpz_tdiv_q_ui (Fab, Fab, proj_primes);
	      
	      do 
		{
		  q = iscomposite (proj_primes);
		  if (q == 0)
		    q = proj_primes;
		  proj_primes /= q;
		  if (nr_primes_a < max_nr_primes)
		    primes_a[nr_primes_a++] = q;
		} while (proj_primes > 1);
	    }

	  /* See if any of the "useful primes" divide this norm */
	  for (k = 0; useful_primes_a[k] != 0; k++)
	  {
	      const fbprime_t q = useful_primes_a[k];
	      if (mpz_divisible_ui_p (Fab, q))
	        {
	          if (nr_primes_a < max_nr_primes)
	            primes_a[nr_primes_a++] = q;
		  mpz_tdiv_q_ui (Fab, Fab, q);
	        }
	  }
	  
	  /* Go through the entire factor base */ 
	  for (nextfb = fba; nextfb->p != 0; nextfb = fb_next (nextfb))
            if (trialdiv_one_prime (nextfb->p, Fab, &nr_primes_a, primes_a,
                                    max_nr_primes))
	    {
		if (mpz_cmp_ui (Fab, 1) > 0)
		    cof_a_prp = 1;
		break;
	    }
	  
	  /* Test if the cofactor after trial division is small enough */
	  log2size = mpz_sizeinbase (Fab, 2);
	  if ((double) log2size > 
	      (*poly)->lpba * (*poly)->alambda + SIEVE_PERMISSIBLE_ERROR)
	  {
	      gmp_fprintf (stderr, 
			   "Sieve report %ld, %lu is not smooth on "
			   "algebraic side, cofactor is %Zd with %d bits\n", 
			   a, b, Fab, mpz_sizeinbase (Fab, 2));
	      goto nextreport;
	  }

	  if (cof_a_prp && log2size > (*poly)->lpba)
	  {
	      lp_a_toolarge++;
	      goto nextreport;
	  }

	  if (!cof_a_prp && log2size > (*poly)->mfba)
	  {
	      cof_a_toolarge++;
	      goto nextreport;
	  }

	  /* Do the rational side */

	  while (mpz_even_p (Gab))
	    {
	      if (nr_primes_r < max_nr_primes)
	        primes_r[nr_primes_r++] = 2;
              mpz_div_2exp (Gab, Gab, 1);
	    }
	  
          trialdiv_one_prime (3, Gab, &nr_primes_r, primes_r, max_nr_primes);

	  for (k = 0; useful_primes_r[k] != 0; k++)
	  {
	      const fbprime_t q = useful_primes_r[k];
	      if (mpz_divisible_ui_p (Gab, q))
	        {
	          if (nr_primes_r < max_nr_primes)
	            primes_r[nr_primes_r++] = q;
		  mpz_tdiv_q_ui (Gab, Gab, q);
	        }
	  }
	  
          for (q = 5, incrq = 2; q <= (*poly)->rlim; 
               q += incrq, incrq = 4 - incrq)
            if (trialdiv_one_prime (q, Gab, &nr_primes_r, primes_r,
                                    max_nr_primes))
	    {
		if (mpz_cmp_ui (Gab, 1) > 0)
		    cof_r_prp = 1;
		break;
	    }

	  if ((double) mpz_sizeinbase (Gab, 2) > 
	      (*poly)->lpbr * (*poly)->rlambda + SIEVE_PERMISSIBLE_ERROR)
	  {
	      gmp_fprintf (stderr, 
			   "Sieve report %ld, %lu is not smooth on "
			   "rational side, cofactor is %Zd with %d bits\n", 
			   a, b, Gab, mpz_sizeinbase (Gab, 2));
	  } else if (cof_r_prp && 
		     mpz_sizeinbase (Gab, 2) > (unsigned)(*poly)->lpbr)
	  {
	      lp_r_toolarge++;
	  } else if (!cof_r_prp && 
		     mpz_sizeinbase (Gab, 2) > (unsigned)(*poly)->mfbr)
	  {
	      cof_r_toolarge++;
	  } else {
	      /* Now print the relations */
	      printf ("%ld,%lu:", a, b);
	      for (k = 0; k < nr_primes_r; k++)
		  printf ("%x%s", primes_r[k], k+1==nr_primes_r?":":",");
	      for (k = 0; k < nr_primes_a; k++)
		  printf ("%x%s", primes_a[k], k+1==nr_primes_a?"\n":",");
	  }
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
	    "alg prp %d, alg composite %d, rat prp %d, rat composite %d\n", 
	    lp_a_toolarge, cof_a_toolarge, lp_r_toolarge, cof_r_toolarge);
  }
#endif
}


int
main (int argc, char **argv)
{
  long amin, amax;
  unsigned long bmin, bmax, b;
  fbprime_t *useful_primes_a, *useful_primes_r, useful_threshold;
  factorbase_t *fba, *fbr;
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
      mp_poly_print ((*cpoly)->f, (*cpoly)->degree, "f(x) =");
      printf ("\n");
      mp_poly_print ((*cpoly)->g, 1, "g(x) =");
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

  fba = fb_read (fbfilename, log_scale, verbose);
  if (fba == NULL)
    {
      fprintf (stderr, "Could not read factor base\n");
      exit (EXIT_FAILURE);
    }

  fbr = fb_make_linear ((*cpoly)->g, (fbprime_t) (*cpoly)->rlim, log_scale,
			verbose);
  if (fbr == NULL)
    {
      fprintf (stderr, 
	       "Could not generate factor base for linear polynomial\n");
      exit (EXIT_FAILURE);
    }
  
  if (0)
  {
      factorbase_t *fbptr;
      for (fbptr = fbr; fbptr->p != 0; fbptr = fb_next (fbptr))
	  fb_print_entry(fbptr);
  }

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

  useful_threshold = fb_maxprime (fba) / 10;
  if (verbose)
    printf ("# Threshold for useful primes is " FBPRIME_FORMAT "\n", 
	    useful_threshold);

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
	  printf ("#Sieving algebraic side\n");

      reports_a_nr = 
	sieve_one_side (sievearray, fba, reports_a, reports_a_len, 
			report_a_threshold, useful_primes_a, useful_threshold,
			useful_length_a, amin, amax, b, proj_roots, log_scale, 
			dpoly_a, deg, verbose);

      proj_roots = mpz_gcd_ui (NULL, (*cpoly)->g[1], b);
      if (verbose)
	printf ("# Projective roots for b = %lu on rational side are: %lu\n", 
	        b, proj_roots);

      if (verbose)
	  printf ("#Sieving rational side\n");

      reports_r_nr = 
	sieve_one_side (sievearray, fbr, reports_r, reports_r_len, 
			report_r_threshold, useful_primes_r, useful_threshold,
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
			  fba, verbose);
    }


  free (fba);
  free (useful_primes_a);
  free (useful_primes_r);
  free (reports_a);
  free (reports_r);
  free (sievearray);

  return 0;
}
