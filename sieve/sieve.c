/* Apply the sieve updates to sievearray */
/* Output: sievearray. Must have enough memory allocated */
/* Inputs: amin, amax, b the $a$ range and the line $b$ to sieve */
/* factorbase, polynomial */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <asm/msr.h>
#define MAXDEGREE 10
#ifdef WANT_ASSERT
  #include <assert.h>
  #define ASSERT(x) assert(x)
#else
  #define ASSERT(x)
#endif


typedef uint32_t fbprime_t;
#define FBPRIME_FORMAT "%u"
typedef fbprime_t fbroot_t;
#define FBROOT_FORMAT "%u"
typedef unsigned long largeprime_t;
const char *largeprime_format = "%lu";

typedef struct {
  fbprime_t p;            /* A prime or a prime power */
  unsigned char plog;     /* logarithm (to some suitable base) of this prime */
  unsigned char nr_roots; /* how many roots there are for this prime */
  unsigned char size;     /* The length of the struct in bytes */
  unsigned char dummy[1]; /* For dword aligning the roots */
  fbroot_t roots[0];      /* the actual length of this array is determined
                             by nr_roots */
} factorbase32_t;


/*****************************************************************
 *                      Functions for calculus                   *
 *****************************************************************/

/* Evaluate poly */

double 
calc_poly_eval (const double *poly, const int deg, const double x)
{
    double r;
    int i;

    ASSERT (deg <= MAXDEG);
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
    ASSERT (deg <= MAXDEG);
    for (i = 1; i <= deg; i++)
	deriv[i - 1] = (double) i * poly[i];
}

/* Newton root finding */
double
calc_newton_root (const double *poly, const int deg, const double start)
{
    double x, lastx, deriv[MAXDEGREE + 1];

    ASSERT (deg <= MAXDEG);
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
    
    ASSERT(deg <= MAXDEG);
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


/* Returns 0 if n is prime, otherwise the smallest prime factor of n */
static fbprime_t
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

  if (a == 0)
    return b;

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


static unsigned long
mulmod (const unsigned long a, const unsigned long b, const unsigned long m)
{
  unsigned long _r, _a = a;

  __asm__ ( "mulq %2\n\t"
            "divq %3"
            : "=&d" (_r), "+a" (_a)
            : "rm" (b), "rm" (m)
            : "cc");

  return _r;
}


/*****************************************************************
 *                Functions for the factor base                  *
 *****************************************************************/

/* Hack to get around C's automatic multiplying constants to be added to 
   pointers by the pointer's base data type size */
static inline factorbase32_t *
fb_skip (const factorbase32_t *fb, const size_t s)
{
  return (factorbase32_t *)((char *)fb + s);
}

static inline factorbase32_t *
fb_next (const factorbase32_t *fb)
{
  return (factorbase32_t *)((char *)fb + fb->size);
}


static size_t
fb_entrysize (const factorbase32_t *fb)
{
  return (sizeof (factorbase32_t) + fb->nr_roots * sizeof (fbroot_t));
}

void 
print_fb_entry (factorbase32_t *fb)
{
  int i;
  printf (FBPRIME_FORMAT ": ", fb->p);
  for (i = 0; i + 1 < fb->nr_roots; i++)
    printf (FBROOT_FORMAT ",", fb->roots[i]);
  printf (FBROOT_FORMAT "\n", fb->roots[i]);
}


/* Add fb_add to (void *)fb + fbsize. If a realloc failed, returns NULL.
   fb_add->size need not be set by caller, this function does it */

factorbase32_t *
add_to_fb (factorbase32_t *fb, size_t *fbsize, size_t *fballoc,
	   const size_t allocblocksize, factorbase32_t *fb_add)
{
  const size_t fb_addsize  = fb_entrysize (fb_add); 
  factorbase32_t *newfb = fb;

  ASSERT(fb_addsize <= allocblocksize); /* Otherwise we still might not have
					   enough mem after the realloc */

  /* Do we need more memory for fb? */
  if (*fballoc < *fbsize + fb_addsize)
    {
      *fballoc += allocblocksize;
      newfb = (factorbase32_t *) realloc (fb, *fballoc);
      if (newfb == NULL)
	{
	  fprintf (stderr, 
		   "Could not reallocate factor base to %lu bytes\n",
		   (unsigned long) *fballoc);
	  return NULL;
	}
    }
  memcpy (fb_skip(newfb, *fbsize), fb_add, fb_addsize);
  fb_skip(newfb, *fbsize)->size = fb_addsize;
  
  *fbsize += fb_addsize;

  return newfb;
}

static factorbase32_t *
fb_find_p (factorbase32_t *fb, const fbprime_t p)
{
  while (fb->p > 0 && fb->p < p) /* Assumes fb is sorted in asc. order */
    fb = fb_next (fb);

  if (fb->p == p)
    return fb;

  return NULL;
}

factorbase32_t *
read_fb (const char *filename, const double log_scale)
{
  factorbase32_t *fb = NULL, *fb_cur, *fb_new;
  FILE *fbfile;
  size_t fbsize = 0, fballoc = 0;
  const size_t linesize = 255;
  size_t linelen;
  char line[linesize];
  char *lineptr;
  unsigned int i;
  int ok;
  unsigned long linenr = 0;
  const size_t allocblocksize = 1<<20; /* Allocate in MB chunks */
  fbprime_t p;

  fbfile = fopen (filename, "r");
  if (fbfile == NULL)
    {
      fprintf (stderr, "Could not open file %s for reading\n", filename);
      return NULL;
    }

  fb_cur = (factorbase32_t *) malloc (sizeof (factorbase32_t) + 
				      MAXDEGREE * sizeof(fbroot_t));
  if (fb_cur == NULL)
    {
      fclose (fbfile);
      return NULL;
    }

  while (!feof(fbfile))
    {
      if (fgets (line, linesize, fbfile) == NULL)
	break;
      linenr++;
      linelen = strlen (line);
      if (linelen > 0 && line[linelen - 1] == '\n')
	linelen--; /* Remove newline */
      for (i = 0; i < linelen; i++) /* Skip comments */
	if (line[i] == '#')
	  break;
      linelen = i;
      if (linelen == 0) /* Skip empty/comment lines */
	continue;
      line[linelen] = '\0';

      /* Parse the line */
      lineptr = line;
      fb_cur->p = strtoul (lineptr, &lineptr, 10);
      ok = 0;
      if (fb_cur->p == 0)
	fprintf (stderr, "read_fb: prime is not an integer on line %lu\n", 
		 linenr);
      else if (*lineptr != ':')
	fprintf (stderr, "read_fb: prime is not followed by colon on line %lu",
		 linenr);
      else
	ok = 1;
      if (!ok)
	continue;

      p = iscomposite (fb_cur->p);
      if (p == 0) /* It's a prime, do a normal log */
	fb_cur->plog = (unsigned char) floor (log ((double) fb_cur->p) * 
					      log_scale + 0.5);
      else
	{
	  /* Take into account the logs of the smaller powers of p that
	     have been added already. We should not just use 
	     fb_cur->plog = log(p) as that would allow rounding errors to 
	     accumulate. */
	  double oldlog;
	  oldlog = floor (log ((double) (fb_cur->p / p)) * log_scale + 0.5);
	  fb_cur->plog = (unsigned char) floor (log ((double) fb_cur->p) * 
						log_scale - oldlog + 0.5);
	}

      lineptr++; /* Skip colon */
      ok = 1;
      fb_cur->nr_roots = 0;
      /* Read roots */
      while (ok && *lineptr != '\0' && fb_cur->nr_roots <= MAXDEGREE)
	{
	  fb_cur->roots[fb_cur->nr_roots++] = strtoul (lineptr, &lineptr, 10);
	  if (*lineptr != '\0' && *lineptr != ',')
	    ok = 0;
	  if (*lineptr == ',')
	    lineptr++;
	}
      
      fb_new = add_to_fb (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
      fb = fb_new;
    }      

  fb_cur->p = 0;
  fb->nr_roots = 0;
  
  fb_new = add_to_fb (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
  if (fb_new == NULL)
    free (fb);
  fb = fb_new;
  
  fclose (fbfile);
  free (fb_cur);
  return fb;
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

  n = (unsigned char) floor (log (fabs (r)) * log_scale 
			     - log_proj_roots + 0.5);
  /* printf ("Norm at x = %.0f is %.0f, rounded log is %d\n", x, r, (int) n); */
  return n;
}


/* Very slow but thorough way of computing norms */

void
compute_norms (unsigned char *sievearray, const long amin, const long amax, 
	       const unsigned long b, const double *poly, const int deg, 
	       const double proj_roots, const double log_scale)
{
    double f[MAXDEGREE + 1]; /* Poly in a for a given fixed b */
    double bpow;
    const double log_proj_roots = log(proj_roots) * log_scale;
    long a, a2;
    int i;
    unsigned char n1, n2;
    const int stride = 256;

    bpow = 1.;
    for (i = deg; i >= 0; i--)
    {
	f[i] = poly[i] * bpow;
	bpow *= (double) b;
    }

    n1 = log_norm (f, deg, (double) amin, log_scale, log_proj_roots);
    
    a2 = amin;
    for (a = amin; a <= amax;)
    {
      ASSERT (a = a2);
      ASSERT (n1 = log_norm (f, deg, (double) a, log_scale, log_proj_roots));
      /* We'll cover [a ... a2 - 1] now */
      a2 = (a + stride <= amax + 1) ? a + stride : amax + 1;

      n2 = log_norm (f, deg, (double) a2, log_scale, log_proj_roots);

      if (n1 == n2) /* Let's assume the log norm is n1 everywhere in this
		       interval */
	{
	  /* printf ("log_c(F(%ld, %lu)) == log_c(F(%ld, %lu)) == %d\n",
	     a, b, a2, b, n1); */
	  memset (sievearray + a - amin, n1, a2 - a);
	}
      else
	{
	  /* n1 and n2 are different. Do each a up to (exclusive) a2 
	     individually */
	  /* printf ("log_c(F(%ld, %lu)) == %d != log_c(F(%ld, %lu)) == %d\n",
	     a, b, n1, a2, b, n2); */
	  sievearray[a - amin] = n1;
	  for ( ; a < a2; a++)
	    sievearray[a - amin] = log_norm (f, deg, (double) a, log_scale, 
					     log_proj_roots);
	}
      a = a2;
      n1 = n2;
    }
}

/* sievearray must be long enough to hold amax-amin+1 chars.
   proj_roots is the rounded log of the projective roots */

void 
sieve (unsigned char *sievearray, factorbase32_t *fb, 
       const long amin, const long amax, const unsigned long b, 
       const unsigned char threshold, fbprime_t *useful_primes,
       const fbprime_t useful_threshold)
{
  uint32_t i, amin_p, a, p, d;
  unsigned char plog;
  const unsigned long l = amax - amin + 1;

  ASSERT (amax > amin);

  /* Do the sieving */

  while (fb->p > 0)
    {
      ASSERT (fb->size <= fb_entrysize(*fb))
      p = fb->p;
      plog = fb->plog;

      /* This modular reduction should be simplified somehow. Do it once
         and store it in fb? */
      if (amin < 0)
        amin_p = p - ((unsigned long)(-amin) % p);
      else
        amin_p = (unsigned long) amin % p;

      for (i = 0; i < fb->nr_roots; i++)
        {
          /* Find first index in sievearray where p on this root divides.
             So we want a/b == r (mod p) <=> a == br (mod p). Then we want
             d so that amin + d == a (mod p) <=> d == a - amin (mod p)
          */

          a = mulmod (b, fb->roots[i], p); /* We have r_i < p, hence b*r_i/p 
	  its in register and there is no division overflow. If we keep r_i
          in Montgomery representation, a single mul/REDC will compute the
          product, reduce it mod p and return it as an integer. TBD. */

          d = (a >= amin_p) ? a - amin_p : a + p - amin_p; 
	  /* Now d is the first index into sievearray where p divides. */

          /* Now update the sieve array. There is no partitioning, blocking, 
             bucket sieving or anything atm. */
          for (; d < l; d += p)
	    {
	      unsigned char k;
	      k = sievearray[d] - plog;
	      sievearray[d] = k;
	      if (p >= useful_threshold && k <= threshold)
		*useful_primes++ = p;
	    }
        }
      
      /* Move on to the next factor base prime */
      fb = fb_next (fb);
    }
  *useful_primes++ = 0;
}



int
main (int argc, char **argv)
{
  long amin, amax, a;
  unsigned long proj_roots; /* = gcd (b, leading coefficient) */
  unsigned long bmin, bmax, b;
  long long tsc1, tsc2;
  fbprime_t *useful_primes;
  factorbase32_t *fb;
  char *fbfilename = NULL;
  unsigned char *sievearray;
  int threshold = 10;
  int verbose = 0;
  int deg;
  int i;
  double dpoly[MAXDEGREE];
  const double log_scale = 1. / log (2.); /* Lets use log_2() for a start */
  mpz_t poly[MAXDEGREE], scaled_poly[MAXDEGREE], Fab;


  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 1 && strcmp (argv[1], "-v") == 0)
	{
	  verbose++;
	  argc--;
	  argv++;
	}
      else if (argc > 2 && strcmp (argv[1], "-thres") == 0)
	{
	  threshold = atoi (argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-fb") == 0)
	{
	  fbfilename = argv[2];
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

  if (fbfilename == NULL)
    {
      fprintf (stderr, 
	       "Please specify a factor base file with the -fb option\n");
      exit (EXIT_FAILURE);
    }

  fb = read_fb (fbfilename, log_scale);

  if (fb == NULL)
    {
      fprintf (stderr, "Could not read factor base\n");
      exit (EXIT_FAILURE);
    }

  deg = 5;
  /* FIXME: Read poly from command line or file */
  for (i = 0; i <= deg; i++)
    {
      mpz_init (poly[i]);
      mpz_init (scaled_poly[i]);
    }
  mpz_init (Fab);
  mpz_set_si (poly[4], 5017309194362523L);
  mpz_set_si (poly[3], -1406293661386525L);
  mpz_set_si (poly[2], -1131155401311965L);
  mpz_set_si (poly[1], 4737694118287353L);
  mpz_set_si (poly[0], -3415040824020545L);
  for (i = 0; i <= deg; i++)
    dpoly[i] = mpz_get_d (poly[i]);


#if 0
  if (verbose)
    for (s = 0; fb_skip(fb, s)->p != 0; s += fb_entrysize (fb_skip (fb, s)))
      print_fb_entry (fb_skip(fb, s));
#endif

  sievearray = (unsigned char *) malloc ((amax - amin + 1) * sizeof (char));
  useful_primes = (fbprime_t *) malloc ((amax - amin + 1) * sizeof (fbprime_t));

  for (b = bmin; b <= bmax; b++)
    {
      unsigned long t;

      /* See what projective roots we have */
      proj_roots = mpz_gcd_ui (NULL, poly[deg], b);
      if (verbose)
	printf ("# Projective roots for b = %lu are %lu\n", b, proj_roots);

      /* Remove the roots of primes that divide b, and of the powers of those 
	 primes, from factor base */
      rdtscll (tsc1);
      t = b;
      while (t > 1)
	{
	  factorbase32_t *fb_del;
	  unsigned long p, ppow;
	  p = iscomposite (t);
	  if (p == 0)
	    p = t;
	  do t /= p; while (t % p == 0);
	  ppow = p;
	  while ((fb_del = fb_find_p (fb, ppow)) != NULL)
	    {
	      if (verbose)
		{
		  printf ("# Temporarily removing from factor base ");
		  print_fb_entry (fb_del);
		}
	      fb_del->nr_roots = 0;
	      ppow *= p;
	    }
	}
      rdtscll (tsc2);
      if (verbose)
	printf ("# Removing primes from fb took %lld clocks\n", tsc2 - tsc1);

      rdtscll (tsc1);
      compute_norms (sievearray, amin, amax, b, dpoly, deg, proj_roots, 
		     log_scale);
      rdtscll (tsc2);
      if (verbose)
	printf ("# Computing norms took %lld clocks\n", tsc2 - tsc1);

      rdtscll (tsc1);
      sieve (sievearray, fb, amin, amax, b, threshold, useful_primes, 500);
      rdtscll (tsc2);
      if (verbose)
	printf ("# Sieving took %lld clocks\n", tsc2 - tsc1);
      
      
      /* Put the roots back in */
      rdtscll (tsc1);
      t = b;
      while (t > 1)
	{
	  factorbase32_t *fb_restore;
	  unsigned long p, ppow;
	  p = iscomposite (t);
	  if (p == 0)
	    p = t;
	  t /= p;
	  ppow = p;
	  while ((fb_restore = fb_find_p (fb, ppow)) != NULL)
	    {
	      fb_restore->nr_roots = 
		(fb_restore->size - sizeof (factorbase32_t)) 
		/ sizeof (fbroot_t);
	      if (verbose)
		{
		  printf ("# Restored to factor base ");
		  print_fb_entry (fb_restore);
		}
	      ppow *= p;
	    }
	}
      rdtscll (tsc2);
      if (verbose)
	printf ("# Restoring primes to fb took %lld clocks\n", tsc2 - tsc1);
      
      if (verbose && *useful_primes != 0)
	{
	  printf ("# Useful primes were: ");
	  while (*useful_primes != 0)
	    printf (FBPRIME_FORMAT " ", *useful_primes++);
	  printf ("\n");
	}


      /* Not trial factor the candidate relations */

      /* Multiply f_i by b^(deg-i) and put in scaled_poly */
      mp_poly_scale (scaled_poly, poly, deg, b, -1); 
      /* The &poly[0] is ugly, but otherwise gcc complains about not being 
	 able to convert pointer type. FIXME: find out what exactly irks gcc */

      rdtscll (tsc1);
      for (a = amin; a <= amax; a++)
	{
	  /* The + 10 is to deal with accumulated rounding that might have
	     cause the sieve value to drop below 0 and wrap around */
	  if ((unsigned char) (sievearray[a - amin] + 10) <= threshold + 10 && 
	      gcd(labs(a), b) == 1)
	    {
	      factorbase32_t *nextfb;
	      fbprime_t *nextuseful;
	      int first_print = 1;

	      mp_poly_eval (Fab, scaled_poly, deg, a);
	      mpz_abs (Fab, Fab);

	      printf ("%ld %lu: ", a, b);

	      /* See if any of the "useful primes" divide this norm */
	      for (nextuseful = useful_primes; *nextuseful != 0; nextuseful++)
		{
		  if (mpz_tdiv_ui(Fab, *nextuseful) == 0)
		    {
		      printf ("%s" FBPRIME_FORMAT, first_print ? "" : ", ",
			      *nextuseful);
		      mpz_tdiv_q_ui (Fab, Fab, *nextuseful);
		      first_print = 0;
		    }
		}

	      /* Go through the entire factor base */
	      for (nextfb = fb; nextfb->p != 0; 
		   nextfb = fb_next (nextfb))
		{
		  if (mpz_divisible_ui_p (Fab, nextfb->p))
		    {
		      while (mpz_divisible_ui_p (Fab, nextfb->p))
			{
			  printf ("%s" FBPRIME_FORMAT, 
				  first_print ? "" : ", ", nextfb->p);
			  mpz_tdiv_q_ui(Fab, Fab, nextfb->p);
			  first_print = 0;
			}
		      if (mpz_cmp_ui (Fab, 1) == 0 || 
			  mpz_probab_prime_p (Fab, 1))
			break;
		    }
		}
	      if (mpz_cmp_ui (Fab, 1) == 0)
		gmp_printf ("; #log res = %d\n", (int) sievearray[a - amin]);
	      else if (mpz_probab_prime_p (Fab, 1))
		gmp_printf ("%s%Zd; #log res = %d\n", 
			    first_print ? "" : ", ", Fab, 
			    (int) sievearray[a - amin]);
	      else
		{
		  /* FIXME: Factor into large primes */
		  gmp_printf ("; # cof. = %Zd, log res = %d\n", 
			      Fab, (int) sievearray[a - amin]);
		}
	    }
	}
      rdtscll (tsc2);
      if (verbose)
	printf ("# Trial factoring/printing took %lld clocks\n", tsc2 - tsc1);
    }


  free (fb);
  free (sievearray);

  return 0;
}
