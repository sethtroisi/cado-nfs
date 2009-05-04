/*****************************************************************
 *                Functions for the factor base                  *
 *****************************************************************/

#include "sieve_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#define rdtscll(x)
#include "basicnt.h"
#include "fb.h"
#include "utils.h"


void 
fb_fprint_entry (FILE *fd, const factorbase_degn_t *fb)
{
  int i;
  fprintf (fd, "Prime " FBPRIME_FORMAT " with rounded log %d and roots ", 
	   fb->p, (int) fb->plog);
  for (i = 0; i < fb->nr_roots; i++)
    {
      fprintf (fd, FBROOT_FORMAT, fb->roots[i]);
      if (i + 1 < fb->nr_roots)
	fprintf (fd, ", ");
    }
  fprintf (fd, "\n");
}

void 
fb_fprint (FILE *fd, const factorbase_degn_t *fb)
{
  while (fb->p != FB_END)
    {
      fb_fprint_entry (fd, fb);
      fb = fb_next (fb);
    }
}

/* Add fb_add to (void *)fb + fbsize. If a realloc failed, returns NULL.
   fb_add->size need not be set by caller, this function does it */

static factorbase_degn_t *
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

static factorbase_degn_t *
fb_find_p (factorbase_degn_t *fb, const fbprime_t p)
{
  while (fb->p != FB_END && fb->p < p) /* Assumes fb is sorted in asc. order */
    {
      ASSERT (fb->p < fb_next (fb)->p || fb_next (fb)->p == FB_END);
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
#ifndef NDEBUG
  for (k = 1; k < n; k++)
    ASSERT(primes[k - 1] <= primes[k]);
#endif
}

#if 0
static double
fb_log_double (double n, double log_scale, double offset)
{
  return floor (log (n) * log_scale + offset + 0.5);
}
#endif

unsigned char
fb_log (double n, double log_scale, double offset)
{
  return (unsigned char) floor (log (n) * log_scale + offset + 0.5);
}

void
fb_init_firstlog (factorbase_t fb)
{
  int i;
  factorbase_degn_t *fbp;

  for (i = 0; i < 256; i++)
    fb->firstlog[i] = NULL;

  for (fbp = fb->fullfb; fbp->p != FB_END; fbp = fb_next (fbp))
    {
      if (fb->firstlog[fbp->plog] == NULL)
	fb->firstlog[fbp->plog] = fbp;
    }
}


/* Generate a factor base with primes <= bound for a linear polynomial */
/* If projective != 0, adds projective roots (for primes that divide 
   leading coefficient */

factorbase_degn_t *
fb_make_linear (mpz_t *poly, const fbprime_t bound, const double log_scale, 
		const int verbose, const int projective, FILE *output)
{
  fbprime_t p;
  factorbase_degn_t *fb = NULL, *fb_cur, *fb_new;
  size_t fbsize = 0, fballoc = 0;
  const size_t allocblocksize = 1 << 20;
  int had_proj_root = 0;
  
  rdtscll (tsc1);
  fb_cur = (factorbase_degn_t *) malloc (fb_entrysize_uc (1));
  ASSERT (fb_cur != NULL);

  fb_cur->nr_roots = 1;
  fb_cur->size = fb_entrysize_uc (1);

  if (verbose)
    gmp_fprintf (output, 
		"# Making factor base for polynomial g(x) = %Zd * x + %Zd\n",
		poly[1], poly[0]);

  for (p = 2 ; p <= bound; p = getprime (p))
    {
      modulusul_t m;
      residueul_t r1, r2;

      fb_cur->p = p;
      fb_cur->plog = fb_log (fb_cur->p, log_scale, 0.);

      modul_initmod_ul (m, p);
      modul_init_noset0 (r1, m);
      modul_init_noset0 (r2, m);
      
      /* We want f_1 * a + f_0 * b == 0 (mod p^k) 
	 with a, b coprime and f_1, f_0 coprime

	 p^m = gcd (f_1, p^k)

	 p^k | F(a,b) <=>
	 ( f_0 * b == 0 (mod p^m) <=> b == 0 (mod p^m)  && 
	   (f_1/p^m) * a + f_0 * (b/p^m) == 0 (mod p^(k-m)) )
	 
	 p^m = (p^k, f_1) => 
	   p^k | F(a,b) <=> p^m | b && f(a/b) == 0 (mod p^(k-m))
	 
      */

      modul_set_ul_reduced (r2, mpz_fdiv_ui (poly[1], p), m);
      if (modul_is0 (r2, m))
	{
	  /* If p | g1 and !(p | g0), p|G(a,b) <=> p|b so that's a 
	     projective root */
	  if (!projective)
	    continue; /* If we don't do projective roots, just move on to 
			 the next prime */
	  if (verbose)
	    {
	      if (!had_proj_root)
		{
		  fprintf (output, "# Primes with projective roots:");
		  had_proj_root = 1;
		}
	      fprintf (output, " " FBPRIME_FORMAT , p);
	    }
	  fb_cur->roots[0] = p; /* The root is 1/0, which we store as 0/1 + p */
	  if (p % 2 != 0)
	    fb_cur->invp = - modul_invmodlong (modul_getmod_ul (m));
	}
      else
	{
	  /* Affine root */
	  modul_set_ul_reduced (r1, mpz_fdiv_ui (poly[0], p), m);

	  /* We want g_1 * a + g_0 * b == 0 <=> a/b == - g0 / g1 */
	  modul_inv (r2, r2, m); /* r2 = 1 / g1 */

	  modul_mul (r2, r1, r2, m); /* r2 = g0 / g1 */
	  modul_neg (r2, r2, m); /* r2 = - g0 / g1 */

#ifndef NDEBUG
	  {
	    residueul_t r3;
	    modul_init_noset0 (r3, m);
	    modul_set_ul_reduced (r3, mpz_fdiv_ui (poly[1], p), m);
	    modul_mul (r3, r3, r2, m); /* g1 * (- g0 / g1) */
	    modul_add (r3, r3, r1, m); /* g1 * (- g0 / g1) + g0 */
	    ASSERT (modul_get_ul (r3, m) == 0);
	    modul_clear (r3, m);
	  }
#endif

	  if (p % 2 != 0)
	    fb_cur->invp = - modul_invmodlong (modul_getmod_ul (m));
	  fb_cur->roots[0] = modul_get_ul (r2, m);
	}

      modul_clear (r1, m);
      modul_clear (r2, m);
      modul_clearmod (m);

      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
      /* fb_fprint_entry (stdout, fb_cur); */
      fb = fb_new;
      
      /* FIXME: handle prime powers */
    }
  getprime (0); /* free prime iterator */

  if (fb != NULL) /* If nothing went wrong so far, put the end-of-fb mark */
    {
      fb_cur->p = FB_END;
      fb_cur->invp = -1L;
      fb_cur->nr_roots = 0;
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	free (fb);
      fb = fb_new;
    }

  free (fb_cur);

  if (had_proj_root)
    fprintf (output, "\n");

  rdtscll (tsc2);

  return fb;
}

#if 0
static int fb_comp(const void *a, const void *b) {
  factorbase_degn_t *pa, *pb;
  pa = (factorbase_degn_t *)a;
  pb = (factorbase_degn_t *)b;
  if (pa->p < pb->p) 
    return -1;
  if (pa->p > pb->p) 
    return 1;
  return 0;
}
#endif

#if 0
/* Generate a factor base with prime powers <= bound for a linear polynomial */
static factorbase_degn_t *
fb_make_linear_powers (mpz_t *poly, const fbprime_t bound, const double log_scale, 
		const int verbose)
{
  fbprime_t p;
  factorbase_degn_t *fb = NULL, *fb_cur, *fb_new;
  size_t fbsize = 0, fballoc = 0;
  const size_t allocblocksize = 1 << 20;
  
  rdtscll (tsc1);
  fb_cur = (factorbase_degn_t *) malloc (fb_entrysize_uc (1));
  ASSERT (fb_cur != NULL);

  fb_cur->nr_roots = 1;
  fb_cur->size = fb_entrysize_uc (1);

  if (verbose)
    {
      printf ("# Making factor base for polynomial g(x) = ");
      printf ("\n");
    }

  for (p = 2 ; p <= bound; p = getprime (p))
    {
      modulusul_t m;
      residueul_t r1, r2;
      uint64_t llq;
      fbprime_t q;
      unsigned char old_logpk;

      llq = (uint64_t)p;
      old_logpk = 0;
      while (llq <= (uint64_t)bound) {
	q = (fbprime_t)llq;
	if (q > p)
	  break;

	fb_cur->p = q;
	// TODO: adjust that!
	if (old_logpk == 0) {
	  old_logpk = fb_log (fb_cur->p, log_scale, 0.);
	  fb_cur->plog = old_logpk;
	} else {
	  double newlog;
	  newlog = fb_log_double (fb_cur->p, log_scale, 0.) - (double)old_logpk;
	  old_logpk = (unsigned char) newlog;
	  fb_cur->plog = old_logpk;
	}

	modul_initmod_ul (m, q);
	modul_init_noset0 (r1, m);
	modul_init_noset0 (r2, m);

	/* We want f_1 * a + f_0 * b == 0 (mod p^k) 
	   with a, b coprime and f_1, f_0 coprime

	   p^m = gcd (f_1, p^k)

	   p^k | F(a,b) <=>
	   ( f_0 * b == 0 (mod p^m) <=> b == 0 (mod p^m)  && 
	   (f_1/p^m) * a + f_0 * (b/p^m) == 0 (mod p^(k-m)) )

	   p^m = (p^k, f_1) => 
	   p^k | F(a,b) <=> p^m | b && f(a/b) == 0 (mod p^(k-m))

*/

	modul_set_ul_reduced (r2, mpz_fdiv_ui (poly[1], q), m);
	/* If p | g1 and !(p | g0), p|G(a,b) <=> p|b and that will be handeled
	   by the projective roots */
	if (modul_get_ul (r2, m) == 0)
	  continue;
	modul_set_ul_reduced (r1, mpz_fdiv_ui (poly[0], q), m);

	/* We want g_1 * a + g_0 * b == 0 <=> a/b == - g0 / g1 */
	if (!modul_inv (r2, r2, m)) /* r2 = 1 / g1 */
	  continue;

	modul_mul (r2, r1, r2, m); /* r2 = g0 / g1 */
	modul_neg (r2, r2, m); /* r2 = - g0 / g1 */

#ifndef NDEBUG
	{
	  residueul_t r3;
	  modul_init_noset0 (r3, m);
	  modul_set_ul_reduced (r3, mpz_fdiv_ui (poly[1], q), m);
	  modul_mul (r3, r3, r2, m); /* g1 * (- g0 / g1) */
	  modul_add (r3, r3, r1, m); /* g1 * (- g0 / g1) + g0 */
	  ASSERT (modul_get_ul (r3, m) == 0);
	  modul_clear (r3, m);
	}
#endif

	if (p % 2 != 0)
	  fb_cur->invp = - modul_invmodlong (m);
	fb_cur->roots[0] = modul_get_ul (r2, m);

	modul_clear (r1, m);
	modul_clear (r2, m);
	modul_clearmod (m);

	fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
	if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
	/* fb_fprint_entry (stdout, fb_cur); */
	fb = fb_new;

	/* get next power of p */
	llq *= (uint64_t)p;
      }
    }
  getprime (0); /* free prime iterator */

  /* Sort the factor base */
  qsort((void *)fb, fbsize / fb_entrysize_uc(1), fb_entrysize_uc(1), fb_comp);


  if (fb != NULL) /* If nothing went wrong so far, put the end-of-fb mark */
    {
      fb_cur->p = FB_END;
      fb_cur->invp = -1L;
      fb_cur->nr_roots = 0;
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	free (fb);
      fb = fb_new;
    }

  free (fb_cur);

  rdtscll (tsc2);

  return fb;
}
#endif

/* Assume q is a prime power p^k with k>=1, return p if k > 1, 0 otherwise. */
static uint32_t
is_prime_power (uint32_t q)
{
  unsigned int maxk, k;
  uint32_t p;

  for (maxk = 0, p = q; p > 1; p /= 2, maxk ++);
  for (k = maxk; k >= 2; k--)
    {
      p = (uint32_t) (pow ((double) q, 1.0 / (double) k) + 0.5);
      if (q % p == 0)
        return p;
    }
  return 0;
}

factorbase_degn_t *
fb_read_addproj (const char *filename, const double log_scale, 
                 const int verbose, const fbprime_t *proj_primes)
{
  factorbase_degn_t *fb = NULL, *fb_cur, *fb_new;
  FILE *fbfile;
  size_t fbsize = 0, fballoc = 0;
  // too small linesize led to a problem with rsa768;
  // it would probably be a good idea to get rid of fgets
  const size_t linesize = 1000;
  size_t linelen;
  char line[linesize];
  char *lineptr;
  int ok;
  unsigned int i;
  unsigned long linenr = 0;
  const size_t allocblocksize = 1<<20; /* Allocate in MB chunks */
  fbprime_t p, q; /* q is factor base entry q = p^k */
  fbprime_t maxprime = 0;
  unsigned long nr_primes = 0;

  fbfile = fopen (filename, "r");
  if (fbfile == NULL)
    {
      fprintf (stderr, "# Could not open file %s for reading\n", filename);
      return NULL;
    }

  fb_cur = (factorbase_degn_t *) malloc (sizeof (factorbase_degn_t) + 
				      MAXDEGREE * sizeof(fbroot_t));
  if (fb_cur == NULL)
    {
      fprintf (stderr, "# Could not allocate memory for factor base\n");
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
      q = strtoul (lineptr, &lineptr, 10);
      ok = 0;
      if (q == 0)
	fprintf (stderr, "# fb_read: prime is not an integer on line %lu\n", 
		 linenr);
      else if (*lineptr != ':')
	fprintf (stderr, "# fb_read: prime is not followed by colon on line %lu",
		 linenr);
      else
	ok = 1;
      if (!ok)
	continue;

      if (q <= maxprime)
	{
	  fprintf (stderr, "# Error, primes in factor base file are "
		   "not in increasing order\n");
	  free (fb);
	  fb = NULL;
	  break;
	}

      /* we assume q is a prime or a prime power */
      p = is_prime_power (q);
      if (p != 0)
        {
          fbprime_t cof = q;
          while (cof % p == 0)
            cof /= p;
          if (cof != 1)
            {
              fprintf (stderr, "# Error, " FBPRIME_FORMAT " on line %lu in "
                       "factor base is not a prime or prime power\n",
                       q, linenr);
              break;                       
            }
        }
      else
	p = q; /* If q is prime, p = q */

      /* We have a valid prime or prime power for q. 
         Maybe there's one (or several) prime with only a projective root 
	 smaller than q that we need to add to the factor base first. */
      while (proj_primes != NULL && *proj_primes != 0 && *proj_primes < q)
        {
          fb_cur->p = *proj_primes++;
          fb_cur->plog = fb_log (fb_cur->p, log_scale, 0.);
          fb_cur->nr_roots = 1;
          fb_cur->roots[0] = fb_cur->p;
	  /* Compute invp */
	  if (fb_cur->p % 2 != 0)
	    {
	      modulusul_t m;
	      
	      modul_initmod_ul (m, fb_cur->p);
	      fb_cur->invp = - modul_invmodlong (modul_getmod_ul (m));
	      modul_clearmod (m);
	    }
	  
          fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
          if (fb_new == NULL)
            {
              free (fb);
              fb = NULL;
              break;
            }
          /* fb_fprint_entry (stdout, fb_cur); */
          fb = fb_new;
          maxprime = fb_cur->p;
          nr_primes++;
        }

      fb_cur->p = q;
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
      while (ok && *lineptr != '\0')
	{
	  if (fb_cur->nr_roots == MAXDEGREE)
            {
              fprintf (stderr, "# Error, too many roots for " FBPRIME_FORMAT 
                      "in factor base line %lu\n", fb_cur->p, linenr);
              ok = 0;
              break;
            }
	  fbroot_t root = strtoul (lineptr, &lineptr, 10);
	  fb_cur->roots[fb_cur->nr_roots++] = root;
	  /* Check if root is properly reduced. Projective roots s/t (mod p^k)
	     with p|t are stored as p + (t/s mod p^k). */
	  if (root >= 2 * fb_cur->p || (root >= fb_cur->p && root % p != 0))
	    ok = 0;
	  if (*lineptr != '\0' && *lineptr != ',')
	    ok = 0;
	  if (*lineptr == ',')
	    lineptr++;
	}

      if (ok && proj_primes != NULL && fb_cur->p == proj_primes[0])
        {
          int i;
          /* Was the projective root already given in the factor base file? */
          for (i = 0; i < fb_cur->nr_roots; i++)
            if (fb_cur->roots[i] == fb_cur->p)
              break;
          if (i == fb_cur->nr_roots) /* No, need to add it */
            {
              if (fb_cur->nr_roots == MAXDEGREE)
                fprintf (stderr, "# Error, can't add projective root for " 
                         FBPRIME_FORMAT ", already has %d roots\n", 
                        fb_cur->p, fb_cur->nr_roots);
              else
                fb_cur->roots[fb_cur->nr_roots++] = fb_cur->p;
            }
          proj_primes++;
        }

      if (!ok)
	{
	  fprintf (stderr, "# Incorrect format in factor base file line %lu\n",
		   linenr);
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
      
      /* Compute invp */
      if (fb_cur->p % 2 != 0)
	{
	  modulusul_t m;
	  
	  modul_initmod_ul (m, fb_cur->p);
	  fb_cur->invp = - modul_invmodlong (modul_getmod_ul (m));
	  modul_clearmod (m);
	}

      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
      /* fb_fprint_entry (stdout, fb_cur); */
      fb = fb_new;
      maxprime = fb_cur->p;
      nr_primes++;
    }      

  if (fb != NULL) /* If nothing went wrong so far, put the end-of-fb mark */
    {
      fb_cur->p = FB_END;
      fb_cur->invp = -1L;
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

/* return the total size (in bytes) used by fb */
size_t
fb_size (const factorbase_degn_t *fb)
{
  size_t mem = 0;
  while (fb->p != FB_END)
    {
      mem += sizeof (fbprime_t) + sizeof (unsigned long)
        + 4 * sizeof (unsigned char)
        + sizeof (fbroot_t) * fb->nr_roots;
      fb = fb_next (fb);
    }
  return mem;
}

factorbase_degn_t *
fb_read (const char *filename, const double log_scale, const int verbose)
{
  return fb_read_addproj (filename, log_scale, verbose, NULL);
}

/* For all primes p that divide b, disable p and powers of p in fb */

void
fb_disable_roots (factorbase_degn_t *fb, const unsigned long b, const int verbose)
{
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
	      fb_fprint_entry (stdout, fb_del);
	    }
	  fb_del->nr_roots = 0;
	  ppow *= p;
	}
    }
  rdtscll (tsc2);
}

/* For all primes p that divide b, re-enable p and powers of p in fb */

void
fb_restore_roots (factorbase_degn_t *fb, const unsigned long b, const int verbose)
{
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
	      fb_fprint_entry (stdout, fb_restore);
	    }
	  ppow *= p;
	}
    }
  rdtscll (tsc2);
}

#if 0
static fbprime_t
fb_maxprime (factorbase_degn_t *fb)
{
    fbprime_t maxp = 0;
    factorbase_degn_t *fbptr = fb;
    while (fbptr->p != FB_END)
    {
	ASSERT (maxp <= fbptr->p);
	maxp = fbptr->p;
	fbptr = fb_next (fbptr);
    }
    
    return maxp;
}
#endif

/* Extracts primes < bound from fb->fblarge and store them in
   factorbase_small_t format. fb->fblarge gets updated to point
   at smallest prime not extracted. */

void
fb_extract_small (factorbase_t fb, const fbprime_t bound, 
		  const int lvl, const int verbose)
{
  factorbase_degn_t *fbptr;
  unsigned int i, j, size;

  ASSERT_ALWAYS (lvl >= 0 && lvl <= SIEVE_BLOCKING);

  /* Count how many roots there are among the primes <= bound, i.e. how
     many entries the small fb will need */
  size = 0;
  for (fbptr = fb->fblarge; fbptr->p <= bound && fbptr->p != FB_END; 
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
  for (fbptr = fb->fblarge; fbptr->p < bound && fbptr->p != FB_END; 
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
  fb->fbsmall[lvl][j].p = FB_END;
  fb->fbsmall[lvl][j].root_and_log = 0;
  j++;

  /* Update the pointer to the large prime factor base */
  
  fb->fblarge = fbptr;

  ASSERT (j == size);

  if (verbose)
    printf ("# There are %u entries in the level %d small factor base\n", 
            j, lvl + 1);
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
  for (fbptr = fb->fullfb; fbptr->p != FB_END; fbptr = fb_next (fbptr))
    {
      p = fbptr->p;

      if (p <= lastp)
	{
	  fprintf (stderr, "Prime " FBPRIME_FORMAT " in factorbase is not "
		   "greater than previous prime " FBPRIME_FORMAT " \n", 
		   p, lastp);
	  return 1;
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

      if (p % 2 != 0 && fbptr->invp * (unsigned long)p != ~(0UL))
	{
	  fprintf (stderr, "For prime " FBPRIME_FORMAT " in factorbase, "
		   "invp %lu is not -1/p (mod 2^wordsize)\n", p, fbptr->invp);
	  return 1;
	}

      if (fb_entrysize (fbptr) != fbptr->size)
	{
	  fprintf (stderr, "For prime " FBPRIME_FORMAT " in factorbase, size "
		   "entry is %d instead of %d\n", 
		   p, (int) fbptr->size, (int) fb_entrysize (fbptr));
	  return 1;
	}

      for (i = 0; i < fbptr->nr_roots; i++)
	{
	  modulusul_t m;
	  residueul_t val, r, c;
	  unsigned long res;

	  modul_initmod_ul (m, p);
	  modul_init_noset0 (r, m);
	  modul_init_noset0 (c, m);
	  modul_init (val, m);
	  modul_set_ul (r, fbptr->roots[i], m);

	  if (side == 0)
	    {
	      int j;
	      modul_set_ul_reduced (val, mpz_fdiv_ui (poly->f[poly->degree], p), 
				  m);
	      for (j = 1; j <= poly->degree; j++)
		{
		  modul_mul (val, val, r, m);
		  modul_set_ul_reduced 
		    (c, mpz_fdiv_ui (poly->f[poly->degree - j], p), m);
		  modul_add (val, val, c, m);
		}
	    }
	  else
	    {
	      modul_set_ul (val, mpz_fdiv_ui (poly->g[1], p), m);
	      modul_mul (val, val, r, m);
	      modul_set_ul_reduced (c, mpz_fdiv_ui (poly->g[0], p), m);
	      modul_add (val, val, c, m);
	    }

	  res = modul_get_ul (val, m);

	  modul_clear (c, m);
	  modul_clear (r, m);
	  modul_clear (r, m);
	  modul_clearmod (m);

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

void
fb_clear (factorbase_t fb)
{
  int i;
  free (fb->fullfb);
  fb->fullfb = NULL;
  fb->fblarge = NULL;
  for (i = 0; i < SIEVE_BLOCKING; i++)
    {
      free (fb->fbsmall[i]);
      fb->fbsmall[i] = NULL;
      free (fb->fbinit[i]);
      fb->fbinit[i] = NULL;
    }
}

/* Extracts primes p <= plim with p/nr_roots <= costlim. 
   List ends with FB_END. Allocates memory */
fbprime_t *
fb_extract_bycost (const factorbase_degn_t *fb, const fbprime_t plim, 
                   const fbprime_t costlim)
{
  const factorbase_degn_t *fb_ptr;
  fbprime_t *primes;
  int i;
  
  fb_ptr = fb;
  i = 0;
  while (fb_ptr->p != FB_END && fb_ptr->p <= plim)  
    {
      if (fb_ptr->p <= costlim * fb_ptr->nr_roots)
        i++;
      fb_ptr = fb_next (fb_ptr);
    }

  primes = (fbprime_t *) malloc ((i + 1) * sizeof (fbprime_t));
  ASSERT (primes != NULL);

  fb_ptr = fb;
  i = 0;
  while (fb_ptr->p != FB_END && fb_ptr->p <= plim)  
    {
      if (fb_ptr->p <= costlim * fb_ptr->nr_roots)
        primes[i++] = fb_ptr->p;
      fb_ptr = fb_next (fb_ptr);
    }

  primes[i] = FB_END;
  
  return primes;
}
