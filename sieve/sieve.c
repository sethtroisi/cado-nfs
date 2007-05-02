/* Apply the sieve updates to sievearray */
/* Output: sievearray. Must have enough memory allocated */
/* Inputs: amin, amax, b the $a$ range and the line $b$ to sieve */
/* factorbase, polynomial */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#define MAXDEGREE 10
#ifdef WANT_ASSERT
  #include <assert.h>
  #define ASSERT(x) assert(x)
#else
  #define ASSERT(x)
#endif


typedef uint32_t fbprime_t;
typedef fbprime_t fbroot_t;
typedef unsigned long largeprime_t;

typedef struct {
  fbprime_t p;            /* A prime or a prime power */
  unsigned char plog;     /* logarithm (to some suitable base) of this prime */
  unsigned char nr_roots; /* how many roots there are for this prime */
  unsigned char size;     /* The length of the struct in bytes */
  unsigned char dummy[1]; /* For dword aligning the roots */
  fbroot_t roots[0];      /* the actual length of this array is determined
                             by nr_roots */
} factorbase32_t;


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
    while (abs (x - lastx) > 1.e-5) /* 1.e-5 should not take too long,
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


/* Hack to get around C's automatic multiplying constants to be added to 
   pointers by the pointer's base data type's size */
static inline factorbase32_t *
fb_skip (const factorbase32_t *fb, const size_t s)
{
  return (factorbase32_t *)((void *)fb + s);
}

static size_t
fb_entrysize (const factorbase32_t *fb)
{
  return (sizeof (factorbase32_t) + fb->nr_roots * sizeof (fbroot_t));
}

static factorbase32_t *
fb_find_p (factorbase32_t *fb, const fbprime_t p)
{
  while (fb->p > 0 && fb->p < p) /* Assumes fb is sorted in asc. order */
    fb = fb_skip (fb, fb->size);

  if (fb->p == p)
    return fb;

  return NULL;
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

/* Very slow but thorough way of computing norms */

void
compute_norms (unsigned char *sievearray, const long amin, const long amax, 
	       const unsigned long b, const double *poly, const int deg, 
	       const double proj_roots, const double log_scale)
{
    double f[MAXDEGREE + 1]; /* Poly in a for a given fixed b */
    double bpow, r;
    const double log_proj_roots = log(proj_roots) * log_scale;
    long a;
    int i;

    bpow = 1.;
    for (i = deg; i >= 0; i--)
    {
	f[i] = poly[i] * bpow;
	bpow *= (double) b;
    }

    for (a = amin; a <= amax; a++)
    {
	unsigned char n;
	r = f[deg];
	for (i = deg - 1; i >= 0; i--)
	    r = r * (double) a + f[i];
	n = round (log(fabs(r)) * log_scale - log_proj_roots);
	sievearray[a - amin] = n;
	if (0)
	    printf ("Norm at (%ld, %lu) = %f, rounded log = %d\n", 
		    a, b, r, (int) n);
    }
}

/* sievearray must be long enough to hold amax-amin+1 chars.
   proj_roots is the rounded log of the projective roots */

void 
sieve (unsigned char *sievearray, factorbase32_t *fb, 
       const long amin, const long amax, const unsigned long b)
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
            sievearray[d] -= plog;
        }
      
      /* Move on to the next factor base prime */
      fb = fb_skip (fb, fb->size);
    }
}


void 
print_fb_entry (factorbase32_t *fb)
{
  int i;
  printf ("%lu: ", (unsigned long) fb->p);
  for (i = 0; i + 1 < fb->nr_roots; i++)
    printf ("%lu,", (unsigned long) fb->roots[i]);
  printf ("%lu\n", (unsigned long) fb->roots[i]);
}


/* Add fb_add to (void *)fb + fb_size. If a realloc failed, returns NULL.
   fb_add->size need not be set by caller, this function does it */

void *
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
      newfb = realloc (fb, *fballoc);
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

  fb_cur = malloc (sizeof (factorbase32_t) + MAXDEGREE * sizeof(fbroot_t));
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
	fb_cur->plog = round (log ((double) fb_cur->p) * log_scale);
      else
	{
	  /* Take into account the logs of the smaller powers of p that
	     have been added already. We should not just use 
	     fb_cur->plog = log(p) as that would allow rounding errors to 
	     accumulate. */
	  double oldlog;
	  oldlog = round (log ((double) (fb_cur->p / p)) * log_scale);
	  fb_cur->plog = round (log ((double) fb_cur->p) * log_scale - oldlog);
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


int
main (int argc, char **argv)
{
  factorbase32_t *fb;
  char *fbfilename = NULL;
  long amin, amax, i;
  unsigned long bmin, bmax, b;
  unsigned char *sievearray;
  const double log_scale = 1. / log (2.); /* Lets use log_2() for a start */
  int threshold = 10;
  int verbose = 0;
  unsigned long leading_coeff = 1;
  unsigned long proj_roots; /* = gcd (b, leading coefficient) */
  double poly[MAXDEGREE] = {-3., 0., 0., 0., 0., 0., 1., 0., 0., 0.};
  int deg = 6;

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
      else if (argc > 2 && strcmp (argv[1], "-leading") == 0)
	{
	  leading_coeff = strtoul(argv[2], NULL, 10);;
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

#if 0
  if (verbose)
    for (s = 0; fb_skip(fb, s)->p != 0; s += fb_entrysize (fb_skip (fb, s)))
      print_fb_entry (fb_skip(fb, s));
#endif

  sievearray = malloc ((amax - amin + 1) * sizeof (char));

  for (b = bmin; b <= bmax; b++)
    {
      unsigned long t;

      /* See what projective roots we have */
      proj_roots = gcd (b, leading_coeff);
      if (verbose)
	printf ("Projective roots for b = %lu are %lu\n", b, proj_roots);

      /* Remove the roots of primes that divide b, and of the powers of those 
	 primes, from factor base */
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
		  printf ("Temporarily removing from factor base ");
		  print_fb_entry (fb_del);
		}
	      fb_del->nr_roots = 0;
	      ppow *= p;
	    }
	}

      compute_norms (sievearray, amin, amax, b, poly, deg, proj_roots, 
		     log_scale);
      sieve (sievearray, fb, amin, amax, b);
      
      /* Put the roots back in */
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
		  printf ("Restored to factor base ");
		  print_fb_entry (fb_restore);
		}
	      ppow *= p;
	    }
	}
      
      for (i = amin; i <= amax; i++)
	  if ( ((int) (sievearray[i - amin] + (unsigned char) 10) - 10)
	       <= threshold && gcd(labs(i), b) == 1)
	      printf ("%ld, %lu: %d\n", 
		      i, b, (int) sievearray[i - amin]);
    }


  free (fb);
  free (sievearray);

  return 0;
}
