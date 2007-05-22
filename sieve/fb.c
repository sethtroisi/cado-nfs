/*****************************************************************
 *                Functions for the factor base                  *
 *****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#define ASSERT_ALWAYS(x) assert(x)
#ifdef WANT_ASSERT
  #define ASSERT(x) assert(x)
#else
  #define ASSERT(x)
#endif
#ifdef HAVE_MSRH
  #include <asm/msr.h>
#else
  #define rdtscll(x)
#endif
#include "fb.h"
#include "mod_ul.c"

/* Some prototypes for functions scattered in various files with no 
   corresponding header files. FIXME: clean this up!*/
fbprime_t iscomposite (const fbprime_t);
void mp_poly_print (mpz_t *, int, const char *);

void 
fb_print_entry (factorbase_t *fb)
{
  int i;
  printf ("Prime " FBPRIME_FORMAT " with rounded log %d and roots ", 
	  fb->p, (int) fb->plog);
  for (i = 0; i + 1 < fb->nr_roots; i++)
    printf (FBROOT_FORMAT ", ", fb->roots[i]);
  printf (FBROOT_FORMAT "\n", fb->roots[i]);
}


/* Add fb_add to (void *)fb + fbsize. If a realloc failed, returns NULL.
   fb_add->size need not be set by caller, this function does it */

factorbase_t *
fb_add_to (factorbase_t *fb, size_t *fbsize, size_t *fballoc,
	   const size_t allocblocksize, factorbase_t *fb_add)
{
  const size_t fb_addsize  = fb_entrysize (fb_add); 
  factorbase_t *newfb = fb;

  ASSERT(fb_addsize <= allocblocksize); /* Otherwise we still might not have
					   enough mem after the realloc */

  /* Do we need more memory for fb? */
  if (*fballoc < *fbsize + fb_addsize)
    {
      *fballoc += allocblocksize;
      newfb = (factorbase_t *) realloc (fb, *fballoc);
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

factorbase_t *
fb_find_p (factorbase_t *fb, const fbprime_t p)
{
  while (fb->p > 0 && fb->p < p) /* Assumes fb is sorted in asc. order */
    fb = fb_next (fb);

  if (fb->p == p)
    return fb;

  return NULL;
}

unsigned char
fb_log (double n, double log_scale, double offset)
{
  return (unsigned char) floor (log (n) * log_scale + offset + 0.5);
}

/* Generate a factor base with primes <= bound for a linear polynomial */

factorbase_t *
fb_make_linear (mpz_t *poly, const fbprime_t bound, const double log_scale, 
		const int verbose)
{
  long long tsc1, tsc2;
  fbprime_t p;
  factorbase_t *fb = NULL, *fb_cur, *fb_new;
  size_t fbsize = 0, fballoc = 0;
  const size_t allocblocksize = 1<<20;

  rdtscll (tsc1);
  fb_cur = (factorbase_t *) malloc (fb_entrysize_uc (1));
  ASSERT (fb_cur != NULL);

  fb_cur->nr_roots = 1;
  fb_cur->size = fb_entrysize_uc (1);

  if (verbose)
    {
      printf ("# Making factor base for polynomial g(x) = ");
      mp_poly_print (poly, 1, "");
    }

  p = 2;
  while (p <= bound)
    {
      modulus m;
      residue r1, r2;

      fb_cur->p = p;
      fb_cur->plog = fb_log (fb_cur->p, log_scale, 0.);

      mod_initmod_ul (m, p);
      mod_init (r1, m);
      mod_init (r2, m);
      
      mod_set_ul_reduced (r2, mpz_fdiv_ui (poly[1], p), m);
      /* If p | g1 and !(p | g0), p|G(a,b) <=> p|b and that will be handeled
	 by the projective roots */
      if (mod_get_ul (r2, m) == 0)
	continue;
      mod_set_ul_reduced (r1, mpz_fdiv_ui (poly[0], p), m);

      /* We want g_1 * a + g_0 * b == 0 <=> a/b == - g0 / g1 */
      mod_inv (r2, r2, m); /* r2 = 1 / g1 */

      mod_mul (r2, r1, r2, m); /* r2 = g0 / g1 */
      mod_neg (r2, r2, m); /* r2 = - g0 / g1 */

      fb_cur->roots[0] = mod_get_ul (r2, m);

#ifdef WANT_ASSERT
      {
	residue r3;
	mod_init (r3, m);
	mod_set_ul_reduced (r3, mpz_fdiv_ui (poly[1], p), m);
	mod_mul (r3, r3, r2, m); /* g1 * (- g0 / g1) */
	mod_add (r3, r3, r1, m); /* g1 * (- g0 / g1) + g0 */
	ASSERT (mod_get_ul (r3, m) == 0);
	mod_clear (r3, m);
      }
#endif

      mod_clear (r1, m);
      mod_clear (r2, m);
      mod_clearmod (m);

      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
      /* fb_print_entry (fb_cur); */
      fb = fb_new;
      
      /* FIXME: handle prime powers */

      /* Skip to next prime. FIXME: use a sieve */
      do {p++;} while (iscomposite (p));
    }

  if (fb != NULL) /* If nothing went wrong so far, put the end-of-fb mark */
    {
      fb_cur->p = 0;
      fb_cur->nr_roots = 0;
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	free (fb);
      fb = fb_new;
    }

  rdtscll (tsc2);
#ifdef HAVE_MSRH
  if (verbose)
    printf ("# Generating rational factor base took %lld clocks\n", 
	    tsc2 - tsc1);
#endif

  return fb;
}


factorbase_t *
fb_read (const char *filename, const double log_scale, const int verbose)
{
  factorbase_t *fb = NULL, *fb_cur, *fb_new;
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
  fbprime_t maxprime = 0;
  unsigned long nr_primes = 0;

  fbfile = fopen (filename, "r");
  if (fbfile == NULL)
    {
      fprintf (stderr, "Could not open file %s for reading\n", filename);
      return NULL;
    }

  fb_cur = (factorbase_t *) malloc (sizeof (factorbase_t) + 
				      MAXDEGREE * sizeof(fbroot_t));
  if (fb_cur == NULL)
    {
      fprintf (stderr, "Could not allocate memory for factor base\n");
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
      while (linelen > 0 && isspace (line[i - 1]))
	linelen--; /* Skip whitespace at end of line */
      if (linelen == 0) /* Skip empty/comment lines */
	continue;
      line[linelen] = '\0';

      /* Parse the line */
      lineptr = line;
      fb_cur->p = strtoul (lineptr, &lineptr, 10);
      ok = 0;
      if (fb_cur->p == 0)
	fprintf (stderr, "fb_read: prime is not an integer on line %lu\n", 
		 linenr);
      else if (*lineptr != ':')
	fprintf (stderr, "fb_read: prime is not followed by colon on line %lu",
		 linenr);
      else
	ok = 1;
      if (!ok)
	continue;

      if (fb_cur->p <= maxprime)
	{
	  fprintf (stderr, "Error, primes in factor base file are "
		   "not in increasing order\n");
	  free (fb);
	  fb = NULL;
	  break;
	}

      p = iscomposite (fb_cur->p);
      if (p == 0) /* It's a prime, do a normal log */
	fb_cur->plog = fb_log (fb_cur->p, log_scale, 0.);
      else
	{
	  /* Take into account the logs of the smaller powers of p that
	     have been added already. We should not just use 
	     fb_cur->plog = log(p) as that would allow rounding errors to 
	     accumulate. */
	  double oldlog;
	  oldlog = fb_log (fb_cur->p / p, log_scale, 0.);
	  fb_cur->plog = fb_log (fb_cur->p, log_scale, - oldlog);
	}

      lineptr++; /* Skip colon */
      ok = 1;
      fb_cur->nr_roots = 0;
      /* Read roots */
      while (ok && *lineptr != '\0' && fb_cur->nr_roots <= MAXDEGREE)
	{
	  fbroot_t r = strtoul (lineptr, &lineptr, 10);
	  fb_cur->roots[fb_cur->nr_roots++] = r;
	  if (r >= fb_cur->p) /* Check if root is properly reduced */
	    ok = 0;
	  if (*lineptr != '\0' && *lineptr != ',')
	    ok = 0;
	  if (*lineptr == ',')
	    lineptr++;
	}

      if (!ok)
	{
	  fprintf (stderr, "Incorrect format in factor base file line %lu\n",
		   linenr);
	  continue;
	}

      if (fb_cur->nr_roots > MAXDEGREE)
	{
	  printf ("Error, too many roots for prime " FBPRIME_FORMAT 
		  "in factor base line %lu\n", p, linenr);
	  continue;
	}
      
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
      /* fb_print_entry (fb_cur); */
      fb = fb_new;
      maxprime = fb_cur->p;
      nr_primes++;
    }      

  if (fb != NULL) /* If nothing went wrong so far, put the end-of-fb mark */
    {
      fb_cur->p = 0;
      fb_cur->nr_roots = 0;
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	free (fb);
      fb = fb_new;
    }

  if (fb != NULL && verbose)
    {
      printf ("Factor base sucessfully read, %lu primes, largest was "
	      FBPRIME_FORMAT "\n", nr_primes, maxprime);
    }
  
  fclose (fbfile);
  free (fb_cur);
  return fb;
}

/* For all primes p that divide b, disable p and powers of p in fb */

void
fb_disable_roots (factorbase_t *fb, const unsigned long b, const int verbose)
{
  long long tsc1, tsc2;
  unsigned long t;
  
  /* Remove the roots of primes that divide b, and of the powers of those 
     primes, from factor base */
  rdtscll (tsc1);
  t = b;
  while (t > 1)
    {
      factorbase_t *fb_del;
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
	      fb_print_entry (fb_del);
	    }
	  fb_del->nr_roots = 0;
	  ppow *= p;
	}
    }
  rdtscll (tsc2);
#ifdef HAVE_MSRH
  if (verbose)
    printf ("# Removing primes from fb took %lld clocks\n", tsc2 - tsc1);
#endif
}

/* For all primes p that divide b, re-enable p and powers of p in fb */

void
fb_restore_roots (factorbase_t *fb, const unsigned long b, const int verbose)
{
  long long tsc1, tsc2;
  unsigned long t;

  /* Put the roots back in */
  rdtscll (tsc1);
  t = b;
  while (t > 1)
    {
      factorbase_t *fb_restore;
      unsigned long p, ppow;
      p = iscomposite (t);
      if (p == 0)
	p = t;
      t /= p;
      ppow = p;
      while ((fb_restore = fb_find_p (fb, ppow)) != NULL)
	{
	  fb_restore->nr_roots = 
	    (fb_restore->size - sizeof (factorbase_t)) 
	    / sizeof (fbroot_t);
	  if (verbose)
	    {
	      printf ("# Restored to factor base ");
	      fb_print_entry (fb_restore);
	    }
	  ppow *= p;
	}
    }
  rdtscll (tsc2);
#ifdef HAVE_MSRH
  if (verbose)
    printf ("# Restoring primes to fb took %lld clocks\n", tsc2 - tsc1);
#endif
}

fbprime_t
fb_maxprime (factorbase_t *fb)
{
    fbprime_t maxp = 0;
    factorbase_t *fbptr = fb;
    while (fbptr->p != 0)
    {
	ASSERT (maxp <= fbptr->p);
	maxp = fbptr->p;
	fbptr = fb_next (fbptr);
    }
    
    return maxp;
}

