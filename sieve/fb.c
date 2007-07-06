/*****************************************************************
 *                Functions for the factor base                  *
 *****************************************************************/

#include "config.h"
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
#include "basicnt.h"
#include "fb.h"
#include "mod_ul.c"
#include "sieve_aux.h"

/* Some prototypes for functions scattered in various files with no 
   corresponding header files. FIXME: clean this up!*/
void mp_poly_print (mpz_t *, int, const char *);

void 
fb_print_entry (factorbase_degn_t *fb)
{
  int i;
  printf ("Prime " FBPRIME_FORMAT " with rounded log %d and roots ", 
	  fb->p, (int) fb->plog);
  for (i = 0; i < fb->nr_roots; i++)
    {
      printf (FBROOT_FORMAT, fb->roots[i]);
      if (i + 1 < fb->nr_roots)
	printf (", ");
    }
  printf ("\n");
}


/* Add fb_add to (void *)fb + fbsize. If a realloc failed, returns NULL.
   fb_add->size need not be set by caller, this function does it */

factorbase_degn_t *
fb_add_to (factorbase_degn_t *fb, size_t *fbsize, size_t *fballoc,
	   const size_t allocblocksize, factorbase_degn_t *fb_add)
{
  const size_t fb_addsize  = fb_entrysize (fb_add); 
  factorbase_degn_t *newfb = fb;

  ASSERT(fb_addsize <= allocblocksize); /* Otherwise we still might not have
					   enough mem after the realloc */

  /* Do we need more memory for fb? */
  if (*fballoc < *fbsize + fb_addsize)
    {
      *fballoc += allocblocksize;
      newfb = (factorbase_degn_t *) realloc (fb, *fballoc);
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

factorbase_degn_t *
fb_find_p (factorbase_degn_t *fb, const fbprime_t p)
{
  while (fb->p > 0 && fb->p < p) /* Assumes fb is sorted in asc. order */
    {
      ASSERT (fb->p < fb_next (fb)->p || fb_next (fb)->p == 0);
      fb = fb_next (fb);
    }

  if (fb->p == p)
    return fb;

  return NULL;
}


/* Sort n primes in array *primes into ascending order */

void
fb_sortprimes (fbprime_t *primes, const unsigned int n)
{
  unsigned int k, l, m;
  fbprime_t t;

  for (l = n; l > 1; l = m)
    for (k = m = 0; k < l - 1; k++)
      if (primes[k] > primes[k + 1])
	{
	  t = primes[k];
	  primes[k] = primes[k + 1];
	  primes[k + 1] = t;
	  m = k + 1;
	}
#ifdef WANT_ASSERT
  for (k = 1; k < n; k++)
    ASSERT(primes[k - 1] <= primes[k]);
#endif
}


unsigned char
fb_log (double n, double log_scale, double offset)
{
  return (unsigned char) floor (log (n) * log_scale + offset + 0.5);
}

/* Generate a factor base with primes <= bound for a linear polynomial */

factorbase_degn_t *
fb_make_linear (mpz_t *poly, const fbprime_t bound, const double log_scale, 
		const int verbose)
{
  long long tsc1, tsc2;
  fbprime_t p;
  factorbase_degn_t *fb = NULL, *fb_cur, *fb_new;
  size_t fbsize = 0, fballoc = 0;
  const size_t allocblocksize = 1<<20;

  rdtscll (tsc1);
  fb_cur = (factorbase_degn_t *) malloc (fb_entrysize_uc (1));
  ASSERT (fb_cur != NULL);

  fb_cur->nr_roots = 1;
  fb_cur->size = fb_entrysize_uc (1);

  if (verbose)
    {
      printf ("# Making factor base for polynomial g(x) = ");
      mp_poly_print (poly, 1, "");
      printf ("\n");
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

#ifdef REDC_ROOTS
      if (p % 2 != 0)
	fb_cur->invp = - mod_invmodlong (m);
#endif
      fb_cur->roots[0] = mod_get_ul (r2, m);

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


factorbase_degn_t *
fb_read (const char *filename, const double log_scale, const int verbose)
{
  factorbase_degn_t *fb = NULL, *fb_cur, *fb_new;
  FILE *fbfile;
  size_t fbsize = 0, fballoc = 0;
  const size_t linesize = 255;
  size_t linelen;
  char line[linesize];
  char *lineptr;
  int ok;
  unsigned int i;
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

  fb_cur = (factorbase_degn_t *) malloc (sizeof (factorbase_degn_t) + 
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
          linelen = i;
      while (linelen > 0 && isspace (line[linelen - 1]))
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
      if (p != 0)
        {
          fbprime_t q = fb_cur->p;
          while (q % p == 0)
            q /= p;
          if (q != 1)
            {
              fprintf (stderr, "Error, " FBPRIME_FORMAT " on line %lu in "
                       "factor base is not a prime or prime power\n",
                       fb_cur->p, linenr);
              break;                       
            }
        }
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
	  fbroot_t root = strtoul (lineptr, &lineptr, 10);
	  fb_cur->roots[fb_cur->nr_roots++] = root;
	  if (root >= fb_cur->p) /* Check if root is properly reduced */
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

      /* Sort the roots into ascending order. We assume that fbprime_t and 
         fbroot_t are typecaste compatible. This is ugly. */
      ASSERT_ALWAYS (sizeof (fbprime_t) == sizeof (fbroot_t));
      fb_sortprimes ((fbprime_t *) fb_cur->roots, fb_cur->nr_roots);
      /* Eliminate duplicate roots. Pari's polrootsmod() can produce them */
      {
        unsigned int i = 0, j;
        for (j = 1; j < (unsigned int) fb_cur->nr_roots; j++)
          {
            if (fb_cur->roots[i] == fb_cur->roots[j])
              {
                fprintf (stderr, "Warning, factor base has repeated root " 
                         FBROOT_FORMAT " for prime " FBPRIME_FORMAT ". "
                         "Keeping only one.\n", fb_cur->roots[i], fb_cur->p);
              }
            else
              fb_cur->roots[++i] = fb_cur->roots[j];
          }
          fb_cur->nr_roots = i + 1;
      }
      
#ifdef REDC_ROOTS
      /* Compute invp */
      if (fb_cur->p % 2 != 0)
	{
	  modulus m;
	  
	  mod_initmod_ul (m, fb_cur->p);
	  fb_cur->invp = - mod_invmodlong (m);
	  mod_clearmod (m);
	}
#endif

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
      printf ("# Factor base sucessfully read, %lu primes, largest was "
	      FBPRIME_FORMAT "\n", nr_primes, maxprime);
    }
  
  fclose (fbfile);
  free (fb_cur);
  return fb;
}

/* For all primes p that divide b, disable p and powers of p in fb */

void
fb_disable_roots (factorbase_degn_t *fb, const unsigned long b, const int verbose)
{
  long long tsc1, tsc2;
  unsigned long t;
  
  /* Remove the roots of primes that divide b, and of the powers of those 
     primes, from factor base */
  rdtscll (tsc1);
  t = b;
  while (t > 1)
    {
      factorbase_degn_t *fb_del;
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
fb_restore_roots (factorbase_degn_t *fb, const unsigned long b, const int verbose)
{
  long long tsc1, tsc2;
  unsigned long t;

  /* Put the roots back in */
  rdtscll (tsc1);
  t = b;
  while (t > 1)
    {
      factorbase_degn_t *fb_restore;
      unsigned long p, ppow;
      p = iscomposite (t);
      if (p == 0)
	p = t;
      t /= p;
      ppow = p;
      while ((fb_restore = fb_find_p (fb, ppow)) != NULL)
	{
	  fb_restore->nr_roots = 
	    (fb_restore->size - sizeof (factorbase_degn_t)) 
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
fb_maxprime (factorbase_degn_t *fb)
{
    fbprime_t maxp = 0;
    factorbase_degn_t *fbptr = fb;
    while (fbptr->p != 0)
    {
	ASSERT (maxp <= fbptr->p);
	maxp = fbptr->p;
	fbptr = fb_next (fbptr);
    }
    
    return maxp;
}

/* Extracts primes < bound from fb->fblarge and store them in
   factorbase_small_t format. fb->fblarge gets updated to point
   at smallest prime not extracted. */

void
fb_extract_small (factorbase_t fb, const unsigned int bound, 
		  const int lvl, const int verbose)
{
  factorbase_degn_t *fbptr;
  unsigned int i, j, size;

  ASSERT_ALWAYS (lvl >= 0 && lvl <= SIEVE_BLOCKING);

  /* Count how many roots there are among the primes <= bound, i.e. how
     many entries the small fb will need */
  size = 0;
  for (fbptr = fb->fblarge; fbptr->p <= bound && fbptr->p > 0; 
       fbptr = fb_next (fbptr))
    size += fbptr->nr_roots;

  /* Add one entry for the stop marker */
  size++;

  /* Allocate memory for the small prime factor base and the small prime
     sieve-ready factor base (which may have fewer primes, if some are not 
     coprime to b) */
  fb->fbsmall[lvl] = (factorbase_small_t *) 
                       malloc (size * sizeof (factorbase_small_t));
  fb->fbinit[lvl] = (factorbase_small_inited_t *) 
                        malloc (size * sizeof (factorbase_small_t));
  ASSERT_ALWAYS (fb->fbsmall[lvl] != NULL);
  ASSERT_ALWAYS (fb->fbinit[lvl] != NULL);
  fb->fbsmallsize[lvl] = size;
  fb->fbsmallbound[lvl] = bound;
  
  /* Copy info into the small prime factor base */
  j = 0;
  for (fbptr = fb->fblarge; fbptr->p < bound && fbptr->p > 0; 
       fbptr = fb_next (fbptr))
    {
      for (i = 0; i < fbptr->nr_roots; i++, j++)
	{
	  fb->fbsmall[lvl][j].p = fbptr->p;
	  fb->fbsmall[lvl][j].root_and_log = fbptr->roots[i] + 
	    (fbptr->plog << 24);
	}
    }

  /* Put end marker */
  fb->fbsmall[lvl][j].p = 0;
  fb->fbsmall[lvl][j].root_and_log = 0;
  j++;

  /* Update the pointer to the large prime factor base */
  
  fb->fblarge = fbptr;

  ASSERT (j == size);

  if (verbose)
    printf ("# There are %u entries in the level %d small factor base\n", 
            j, lvl + 1);
}

void
fb_initloc_small (factorbase_small_inited_t *initfb,
		  factorbase_small_t *smallfb, 
		  const long amin, const unsigned long b, const int odd)
{
  uint32_t amin_p;
  fbprime_t p; 
  fbroot_t root;
  uint32_t d;
  unsigned char plog;
  
  for ( ; smallfb->p > 0; smallfb++)
    {
      p = smallfb->p;
      ASSERT (p < 0xffffff);
      root = smallfb->root_and_log & 0xffffff;
      ASSERT (root < p);
      plog = smallfb->root_and_log >> 24;

      if (gcd(p, b) == 1) /* FIXME: Speed up this gcd */
	{
	  amin_p = signed_mod_longto32 (amin, p);
	  /* Replace low 24 bits by first sieve location */
          d = first_sieve_loc (p, root, amin_p, b, odd);
	  initfb->p = p;
	  initfb->loc_and_log = (plog << 24) | d;
#if 0
	  printf ("F(%ld, %lu) %% %d == 0 /* PARI fb_initloc_small */\n",
		  amin + ((fb_inited->loc_and_log & 0xffffff) << odd), b, 
		  fb_inited->p);
#endif
	  initfb ++;
	}
    }

  /* Place stop marker */
  initfb->p = 0;
  initfb->loc_and_log = 0;
}

/* Test if the data in factorbase is correct. If side == 0, checks the roots
   for the algebraic polynomial; otherwise for the rational one */

int
fb_check (factorbase_t fb, cado_poly poly, int side)
{
  factorbase_degn_t *fbptr;
  fbprime_t p, lastp;
  unsigned long q;
  unsigned char i;

  lastp = 1;
  for (fbptr = fb->fullfb; fbptr->p != 0; fbptr = fb_next (fbptr))
    {
      p = fbptr->p;

      if (p <= lastp)
	{
	  fprintf (stderr, "Prime " FBPRIME_FORMAT " in factorbase is not "
		   "greater than previous prime " FBPRIME_FORMAT " \n", 
		   p, lastp);
	  return 0;
	}

      if ((q = iscomposite (p)) != 0)
	{
	  fbprime_t cofac = p;
	  while (cofac > 1 && cofac % (fbprime_t) q == 0)
	    cofac /= (fbprime_t) q;
	  if (cofac != 1)
	    {
	      fprintf (stderr, "Prime " FBPRIME_FORMAT " in factorbase is not "
		       "a prime or prime power\n", p);
	      return 1;
	    }
	}

#ifdef REDC_ROOTS
      if (p % 2 != 0 && fbptr->invp * (unsigned long)p != ~(0UL))
	{
	  fprintf (stderr, "For prime " FBPRIME_FORMAT " in factorbase, "
		   "invp %lu is not -1/p (mod 2^wordsize)\n", p, fbptr->invp);
	  return 1;
	}
#endif

      if (fb_entrysize (fbptr) != fbptr->size)
	{
	  fprintf (stderr, "For prime " FBPRIME_FORMAT " in factorbase, size "
		   "entry is %d instead of %d\n", 
		   p, (int) fbptr->size, (int) fb_entrysize (fbptr));
	  return 1;
	}

      for (i = 0; i < fbptr->nr_roots; i++)
	{
	  modulus m;
	  residue val, r, c;
	  unsigned long res;

	  mod_initmod_ul (m, p);
	  mod_init (r, m);
	  mod_init (c, m);
	  mod_init_set0 (val, m);
	  mod_set_ul (r, fbptr->roots[i], m);

	  if (side == 0)
	    {
	      int j;
	      mod_set_ul_reduced (val, mpz_fdiv_ui (poly->f[poly->degree], p), 
				  m);
	      for (j = 1; j <= poly->degree; j++)
		{
		  mod_mul (val, val, r, m);
		  mod_set_ul_reduced 
		    (c, mpz_fdiv_ui (poly->f[poly->degree - j], p), m);
		  mod_add (val, val, c, m);
		}
	    }
	  else
	    {
	      mod_set_ul (val, mpz_fdiv_ui (poly->g[1], p), m);
	      mod_mul (val, val, r, m);
	      mod_set_ul_reduced (c, mpz_fdiv_ui (poly->g[0], p), m);
	      mod_add (val, val, c, m);
	    }

	  res = mod_get_ul (val, m);

	  mod_clear (c, m);
	  mod_clear (r, m);
	  mod_clear (r, m);
	  mod_clearmod (m);

	  if (res != 0)
	    {
	      fprintf (stderr, "False root in factor base: %c(" FBROOT_FORMAT 
		       ") = %lu (mod " FBPRIME_FORMAT ")\n", 
		       (side == 0 ? 'f' : 'g'), fbptr->roots[i], res, p);
	      return 1;
	    }
	}

      /* FIXME: add the tests for the small factor bases */
    }
  return 0;
}
