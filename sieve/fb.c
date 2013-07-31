/*****************************************************************
 *                Functions for the factor base                  *
 *****************************************************************/

#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#define rdtscll(x)
#include "fb.h"
#include "portability.h"
#include "utils.h"
#include "ularith.h"

typedef struct {
  factorbase_degn_t *fb;
  size_t size, alloc, blocksize;
} fb_buffer_t;

/* strtoul(), but with const char ** for second argument.
   Otherwise it's not possible to do, e.g., strtoul(p, &p, 10) when p is
   of type const char *
*/
static unsigned long int
strtoul_const(const char *nptr, const char **endptr, const int base)
{
  char *end;
  unsigned long r;
  r = strtoul(nptr, &end, base);
  *endptr = end;
  return r;
}

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

/* Initialise a factor base buffer to empty */
static void
fb_buffer_init(fb_buffer_t *fb_buf, const size_t blocksize)
{
    fb_buf->fb = NULL; /* Set to NULL so realloc() allocates */
    fb_buf->size = 0;
    fb_buf->alloc = 0;
    fb_buf->blocksize = blocksize;
}

/* Extend a factorbase buffer, if necessary, to have room for at least
   addsize additional bytes at the end */
static int
fb_buffer_extend(fb_buffer_t *fb_buf, const size_t addsize)
{
  /* Do we need more memory for fb? */
  if (fb_buf->alloc < fb_buf->size + addsize)
    {
      factorbase_degn_t *newfb;
      size_t newalloc = fb_buf->alloc + fb_buf->blocksize;
      ASSERT(addsize <= fb_buf->blocksize); /* Otherwise we still might not have
					       enough mem after the realloc */
      newfb = (factorbase_degn_t *) realloc (fb_buf->fb, newalloc);
      if (newfb == NULL)
	{
	  fprintf (stderr,
		   "Could not reallocate factor base to %zu bytes\n",
		   fb_buf->alloc);
	  return 0;
	}
      fb_buf->alloc = newalloc;
      fb_buf->fb = newfb;
    }
  return 1;
}

/* Add the factor base enty "fb_add" to the factor base buffer "fb_buf". 
   Return 1 on success, 0 if realloc failed.
   fb_add->size need not be set by caller, this function does it */

static int 
fb_buffer_add (fb_buffer_t *fb_buf, const factorbase_degn_t *fb_add)
{
  const size_t fb_addsize = fb_entrysize (fb_add);

  if (!fb_buffer_extend(fb_buf, fb_addsize))
    return 0;

  /* Append the new entry at the end of the factor base */
  factorbase_degn_t *fb_end_ptr = fb_skip(fb_buf->fb, fb_buf->size);
  memcpy (fb_end_ptr, fb_add, fb_addsize);
  fb_end_ptr->size = fb_addsize;

  fb_buf->size += fb_addsize;

  return 1;
}

static int
fb_buffer_finish(fb_buffer_t *fb_buf)
{
  const size_t fb_addsize = fb_entrysize_uc (0);

  if (!fb_buffer_extend(fb_buf, fb_addsize))
    return 0;

  /* Add the end-of-factor-base marker */
  factorbase_degn_t *fb_end_ptr = fb_skip(fb_buf->fb, fb_buf->size);
  fb_end_ptr->p = FB_END;
  fb_end_ptr->invp = -(redc_invp_t)1;
  fb_end_ptr->nr_roots = 0;
  fb_buf->size += fb_addsize;
  return 1;
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

unsigned char
fb_log (double n, double log_scale, double offset)
{
  return (unsigned char) floor (log (n) * log_scale + offset + 0.5);
}

/* Returns floor(log_2(n)) for n > 0, and 0 for n == 0 */
static unsigned int
fb_log_2 (fbprime_t n)
{
  unsigned int k;
  for (k = 0; n > 1; n /= 2, k++);
  return k;
}

/* Return p^e. Trivial exponentiation for small e, no check for overflow */
static fbprime_t
fb_pow (const fbprime_t p, const unsigned int e)
{
    fbprime_t r = 1;

    for (unsigned int i = 0; i < e; i++)
      r *= p;
    return r;
}

/* Let k be the largest integer with q = p^k, return p if k > 1,
   and 0 otherwise */
fbprime_t
fb_is_power (fbprime_t q)
{
  unsigned int maxk, k;
  uint32_t p;

  maxk = fb_log_2(q);
  for (k = maxk; k >= 2; k--)
    {
      double dp = pow ((double) q, 1.0 / (double) k);
      double rdp = trunc(dp + 0.5);
      if (fabs(dp - rdp) < 0.001) {
        p = (uint32_t) rdp ;
        if (q % p == 0) {
          ASSERT (fb_pow (p, k) == q);
          return p;
        }
      }
    }
  return 0;
}

/* Make one factor base entry for a linear polynomial poly[1] * x + poly[0]
   and the prime (power) q. We assume that poly[0] and poly[1] are coprime.
   Non-projective roots a/b such that poly[1] * a + poly[0] * b == 0 (mod q)
   with gcd(poly[1], q) = 1 are stored as a/b mod q.
   If do_projective != 0, also stores projective roots with gcd(q, f_1) > 1,
   but stores the reciprocal root plus p, p + b/a mod p^k.
   If no root was stored (projective root with do_projective = 0) returns 0,
   returns 1 otherwise. */

static int
fb_make_linear1 (factorbase_degn_t *fb_entry, const mpz_t *poly,
		 const fbprime_t p, const unsigned char logp,
		 const int do_projective)
{
  modulusul_t m;
  residueul_t r0, r1;
  modintul_t gcd;
  int is_projective, rc;

  fb_entry->p = p;
  fb_entry->plog = logp;

  modul_initmod_ul (m, p);
  modul_init_noset0 (r0, m);
  modul_init_noset0 (r1, m);

  modul_set_ul_reduced (r0, mpz_fdiv_ui (poly[0], p), m);
  modul_set_ul_reduced (r1, mpz_fdiv_ui (poly[1], p), m);
  modul_gcd (gcd, r1, m);
  is_projective = (modul_intcmp_ul (gcd, 1UL) > 0);

  if (is_projective)
    {
      if (!do_projective)
	{
	  modul_clear (r0, m);
	  modul_clear (r1, m);
	  modul_clearmod (m);
	  return 0; /* If we don't do projective roots, just return */
	}
      modul_swap (r0, r1, m);
    }

  /* We want poly[1] * a + poly[0] * b == 0 <=>
     a/b == - poly[0] / poly[1] */
  rc = modul_inv (r1, r1, m); /* r1 = 1 / poly[1] */
  ASSERT_ALWAYS (rc != 0);
  modul_mul (r1, r0, r1, m); /* r1 = poly[0] / poly[1] */
  modul_neg (r1, r1, m); /* r1 = - poly[0] / poly[1] */

  fb_entry->roots[0] = modul_get_ul (r1, m) + (is_projective ? p : 0);
  if (p % 2 != 0) {
    ASSERT(sizeof(unsigned long) >= sizeof(redc_invp_t));
    fb_entry->invp = - ularith_invmod (modul_getmod_ul (m));
  }

  modul_clear (r0, m);
  modul_clear (r1, m);
  modul_clearmod (m);

  return 1;
}


/* Generate a factor base with primes <= bound and prime powers <= powbound
   for a linear polynomial. If projective != 0, adds projective roots
   (for primes that divide leading coefficient) */

factorbase_degn_t *
fb_make_linear (const mpz_t *poly, const fbprime_t bound,
		const fbprime_t powbound, const double log_scale,
		const int verbose, const int projective, FILE *output)
{
  fb_buffer_t fb_buf;
  fbprime_t p;
  factorbase_degn_t *fb_cur;
  size_t pow_len = 0;
  const size_t allocblocksize = 1 << 20;
  unsigned long logp;
  int had_proj_root = 0;
  fbprime_t *powers = NULL, min_pow = 0; /* List of prime powers that yet
					    need to be included, and the
					    minimum among them */

  rdtscll (tsc1);
  fb_cur = (factorbase_degn_t *) malloc (fb_entrysize_uc (1));
  ASSERT (fb_cur != NULL);

  fb_cur->nr_roots = 1;
  fb_cur->size = fb_entrysize_uc (1);
  fb_buffer_init (&fb_buf, allocblocksize);

  if (verbose)
    gmp_fprintf (output,
		 "# Making factor base for polynomial g(x) = %Zd*x%s%Zd,\n"
		 "# including primes up to " FBPRIME_FORMAT
		 " and prime powers up to " FBPRIME_FORMAT ".\n",
		 poly[1], (mpz_cmp_ui (poly[0], 0) >= 0) ? "+" : "",
                 poly[0], bound, powbound);

  p = 2;
  while (p <= bound)
    {
      fbprime_t q;
      /* Handle any prime powers that are smaller than p */
      if (pow_len > 0 && min_pow < p)
	{
	  size_t i;
	  /* Find this p^k in the list of prime powers and replace it
	     by p^(k+1) if p^(k+1) does not exceed powbound, otherwise
	     remove p^k from the list */
	  for (i = 0; i < pow_len && powers[i] != min_pow; i++);
	  ASSERT_ALWAYS (i < pow_len);
	  q = fb_is_power (min_pow);
	  ASSERT (min_pow / q >= q && min_pow % (q*q) == 0);
	  logp = fb_log (min_pow, log_scale, 0.) -
	         fb_log (min_pow / q, log_scale, 0.);
	  if (powers[i] <= powbound / q)
	    powers[i] *= q; /* Increase exponent */
	  else
	    {
	      for ( ; i < pow_len - 1; i++) /* Remove from list */
		powers[i] = powers[i + 1];
	      pow_len--;
	    }
	  q = min_pow;
	  /* Find new minimum among the prime powers */
	  if (pow_len > 0)
	    min_pow = powers[0];
	  for (i = 1; i < pow_len; i++)
	    if (powers[i] < min_pow)
	      min_pow =  powers[i];
	}
      else
	{
	  q = p;
	  logp = fb_log (q, log_scale, 0.);
	  /* Do we need to add this prime to the prime powers list? */
	  if (q <= powbound / q)
	    {
	      size_t i;
	      powers = realloc (powers, (++pow_len) * sizeof (fbprime_t));
	      ASSERT_ALWAYS(powers != NULL);
	      powers[pow_len - 1] = q*q;
	      /* Find new minimum among the prime powers */
	      min_pow = powers[0];
	      for (i = 1; i < pow_len; i++)
		if (powers[i] < min_pow)
		  min_pow =  powers[i];
	    }
	  p = getprime (p);
	}

      if (!fb_make_linear1 (fb_cur, poly, q, logp, projective))
	continue; /* If root is projective and we don't want those,
		     skip to next prime */

      if (verbose && fb_cur->p > q)
	{
	  if (!had_proj_root)
	    {
	      fprintf (output, "# Primes with projective roots:");
	      had_proj_root = 1;
	    }
	  fprintf (output, " " FBPRIME_FORMAT , q);
	}

      if (!fb_buffer_add (&fb_buf, fb_cur))
	{
	  free (fb_buf.fb);
	  fb_buf.fb = NULL;
	  break;
	}
      /* fb_fprint_entry (stdout, fb_cur); */

      /* FIXME: handle prime powers */
    }

  getprime (0); /* free prime iterator */

  if (fb_buf.fb != NULL) /* If nothing went wrong so far, put the end-of-fb 
			    mark */
    {
      if (!fb_buffer_finish(&fb_buf)) {
	free (fb_buf.fb);
	fb_buf.fb = NULL;
      }
    }

  free (fb_cur);
  free (powers);

  if (had_proj_root)
    fprintf (output, "\n");

  rdtscll (tsc2);

  return fb_buf.fb;
}

/* return the total size (in bytes) used by fb */
size_t
fb_size (const factorbase_degn_t * const fb)
{
  const factorbase_degn_t *fb_cur = fb;
  size_t mem = 0;
  while (1)
    {
      mem += fb_entrysize (fb_cur);
      if (fb_cur->p == FB_END)
        break;
      fb_cur = fb_next (fb_cur);
    }
  return mem;
}

/* return the total size (in number of roots) for this fb */
size_t
fb_nroots_total (const factorbase_degn_t *fb)
{
  size_t n = 0;
  for( ; fb->p != FB_END ; fb = fb_next (fb))
      n += fb->nr_roots;
  return n;
}


/* Extracts primes (not prime powers) p with p/nr_roots <= costlim up to and 
   excluding the first prime p > plim. Prime powers in the factor base are
   ignored. List ends with FB_END. Allocates memory */
fbprime_t *
fb_extract_bycost (const factorbase_degn_t *fb, const fbprime_t plim,
                   const fbprime_t costlim)
{
  fbprime_t *primes = NULL;

  /* First pass counts primes and allocates memory, second pass writes 
     primes to the allocated memory */
  for (int pass = 0; pass < 2; pass++) {
    const factorbase_degn_t *fb_ptr;
    size_t i = 0;
    fbprime_t old_prime = 1;
    for (fb_ptr = fb; fb_ptr->p != FB_END; fb_ptr = fb_next (fb_ptr)) {
      /* Prime powers p^k are neither added to the array of extracted primes, 
         nor do they stop the loop if p^k > plim */
      if (fb_is_power(fb_ptr->p))
        continue;
      if (fb_ptr->p > plim)
        break;
      if (fb_ptr->p <= costlim * fb_ptr->nr_roots && old_prime!=fb_ptr->p) {
        if (pass == 1)
          primes[i] = fb_ptr->p;
        i++;
        old_prime = fb_ptr->p;
      }
    }
    if (pass == 0) {
      primes = (fbprime_t *) malloc ((i + 1) * sizeof (fbprime_t));
      ASSERT_ALWAYS (primes != NULL);
    } else {
      primes[i] = FB_END;
    }
  }

  return primes;
}

/* Remove newline, comment, and trailing space from a line. Write a
   '\0' character to the line at the position where removed part began (i.e.,
   line gets truncated).
   Return length in characters or remaining line, without trailing '\0'
   character.
*/
static size_t
fb_read_strip_comment (char *const line)
{
    size_t linelen, i;

    linelen = strlen (line);
    if (linelen > 0 && line[linelen - 1] == '\n')
        linelen--; /* Remove newline */
    for (i = 0; i < linelen; i++) /* Skip comments */
        if (line[i] == '#') {
            linelen = i;
            break;
        }
    while (linelen > 0 && isspace((int)(unsigned char)line[linelen - 1]))
        linelen--; /* Skip whitespace at end of line */
    line[linelen] = '\0';

    return linelen;
}

/* Read roots from a factor base file line and store them in fb_entry.
   line must point at the first character of the first root on the line.
   linenr is used only for printing error messages in case of parsing error.
*/
static void
fb_read_roots (factorbase_degn_t * const fb_entry, const char *lineptr,
               const unsigned long linenr)
{
    fb_entry->nr_roots = 0;
    while (*lineptr != '\0')
    {
        if (fb_entry->nr_roots == MAXDEGREE) {
            fprintf (stderr,
                    "# Error, too many roots for prime " FBPRIME_FORMAT
                    " in factor base line %lu\n", fb_entry->p, linenr);
            exit(EXIT_FAILURE);
        }
        fbroot_t root = strtoul_const (lineptr, &lineptr, 10);

        if (fb_entry->nr_roots > 0
                && root <= fb_entry->roots[fb_entry->nr_roots-1]) {
            fprintf (stderr,
                "# Error, roots must be sorted in the fb file, line %lu\n",
                linenr);
            exit(EXIT_FAILURE);
        }

        fb_entry->roots[fb_entry->nr_roots++] = root;
        if (*lineptr != '\0' && *lineptr != ',') {
            fprintf(stderr,
                    "# Incorrect format in factor base file line %lu\n",
                    linenr);
            exit(EXIT_FAILURE);
        }
        if (*lineptr == ',')
            lineptr++;
    }

    if (fb_entry->nr_roots == 0) {
        fprintf (stderr, "# Error, no root for prime " FBPRIME_FORMAT
                " in factor base line %lu\n", fb_entry->p, linenr - 1);
        exit(EXIT_FAILURE);
    }
}

/* Parse a factorbase line. Return 1 if the line could be parsed and,
   if the FB entry is for a prime power, the power does not exceed powlim.
   Otherwise return 0. */
static int
fb_parse_line (factorbase_degn_t *const fb_cur, const char * lineptr,
               const double log_scale, const unsigned long linenr,
               const fbprime_t powlim)
{
    long nlogp, oldlogp = 0;
    fbprime_t p, q; /* q is factor base entry q = p^k */

    q = strtoul_const (lineptr, &lineptr, 10);
    if (q == 0) {
        fprintf(stderr, "# fb_read: prime is not an integer on line %lu\n",
                linenr);
        exit (EXIT_FAILURE);
    } else if (*lineptr != ':') {
        fprintf(stderr,
                "# fb_read: prime is not followed by colon on line %lu",
                linenr);
        exit (EXIT_FAILURE);
    }

    lineptr++; /* Skip colon after q */
    int longversion;
    if (strchr(lineptr, ':') == NULL)
        longversion = 0;
    else
        longversion = 1;

    /* we assume q is a prime or a prime power */
    /* NB: a short version is not possible for a prime power, so we
     * do the test only for !longversion */
    p = (longversion) ? fb_is_power (q) : 0;
    ASSERT(isprime(p != 0 ? p : q));
    fb_cur->p = q;
    if (p != 0) {
        if (powlim && fb_cur->p >= powlim)
            return 0;
    } else {
        p = q; /* If q is prime, p = q */
    }

    /* read the multiple of logp, if any */
    /* this must be of the form  q:nlogp,oldlogp: ... */
    /* if the information is not present, it means q:1,0: ... */
    nlogp = 1;
    oldlogp = 0;
    if (longversion) {
        nlogp = strtoul_const (lineptr, &lineptr, 10);
        /*
        if (nlogp == 0) {
            fprintf(stderr, "# Error in fb_read: could not parse the integer after the colon of prime " FBPRIME_FORMAT "\n", q);
            exit (EXIT_FAILURE);
        }*/
        if (*lineptr != ',') {
            fprintf(stderr, "# fb_read: nlogp is not followed by comma on line %lu", linenr);
            exit (EXIT_FAILURE);
        }
        lineptr++; /* skip comma */
        oldlogp = strtoul_const (lineptr, &lineptr, 10);
        /*
        if (oldlogp == 0) {
            fprintf(stderr, "# Error in fb_read: could not parse the integer after the comma of prime " FBPRIME_FORMAT "\n", q);
            exit (EXIT_FAILURE);
        }*/
        if (*lineptr != ':') {
            fprintf(stderr, "# fb_read: oldlogp is not followed by colon on line %lu", linenr);
            exit (EXIT_FAILURE);
        }
        lineptr++; /* skip colon */
    }
    ASSERT (nlogp > oldlogp);

    if (nlogp == 1) /* typical case: this is the first occurrence of p */
        fb_cur->plog = fb_log (p, log_scale, 0.);
    else
      /* p already occurred before, and was taken into account to the power
         'oldlogp', with bias 'ol' since it was rounded to an integer.
         We now want to take into account the extra contribution from
         p^oldlogp to p^nlogp. */
      {
        double ol = fb_log (fb_pow (p, oldlogp), log_scale, 0.);
        fb_cur->plog = fb_log (fb_pow (p, nlogp), log_scale, - ol);
      }

    /* Read roots */
    fb_read_roots(fb_cur, lineptr, linenr);

    return 1;
}

/* Read a factor base file, splitting it into pieces.
   
   Primes and prime powers up to smalllim go into fb_small. If smalllim is 0,
   all primes go into fb_small, and nothing is written to fb_pieces.
   
   If smalllim is not 0, then nr_pieces separate factor bases are made for
   primes/powers > smalllim; factor base entries from the file are written to 
   these pieces in round-robin manner.

   Pointers to the allocated memory of the factorbases are written to fb_small 
   and, if smalllim > 0, to fb_pieces[0, ..., nr_pieces-1].

   Returns 1 if everything worked, and 0 if not (i.e., if the file could not be 
   opened, or memory allocation failed)
*/

int 
fb_read_split (factorbase_degn_t **fb_small, factorbase_degn_t **fb_pieces, 
               const char * const filename, const double log_scale, 
               const fbprime_t smalllim, const int nr_pieces, 
               const int verbose, const fbprime_t lim, const fbprime_t powlim)
{
    fb_buffer_t *fb_bufs; /* fb_bufs[0] is for primes <= smalllim, rest is for 
			     larger primes. If smalllim == 0, then fb_bufs[0] 
			     is used for all primes. */
    factorbase_degn_t *fb_cur;
    size_t nr_buffers;
    FILE *fbfile;
    // too small linesize led to a problem with rsa768;
    // it would probably be a good idea to get rid of fgets
    const size_t linesize = 1000;
    char line[linesize];
    unsigned long linenr = 0;
    const size_t allocblocksize = 1<<20; /* Allocate in MB chunks */
    fbprime_t maxprime = 0;
    unsigned long nr_primes = 0;
    int nextpiece = 0, error = 0;

    fbfile = fopen (filename, "r");
    if (fbfile == NULL) {
        fprintf (stderr, "# Could not open file %s for reading\n", filename);
        return 0;
    }

    fb_cur = (factorbase_degn_t *) malloc (fb_entrysize_uc(MAXDEGREE));
    if (fb_cur == NULL) {
        fprintf (stderr, "# Could not allocate memory for factor base\n");
        fclose (fbfile);
        return 0;
    }

    nr_buffers = (smalllim == 0) ? 1 : nr_pieces + 1;
    fb_bufs = (fb_buffer_t *) malloc (nr_buffers * sizeof(fb_buffer_t));
    if (fb_bufs == NULL) {
        fprintf (stderr, "# %s(): could not allocate memory for fb_bufs\n", 
		 __func__);
	free (fb_cur);
	fclose (fbfile);
	return 0;
    }
    for (size_t i = 0; i < nr_buffers; i++)
      fb_buffer_init(&fb_bufs[i], allocblocksize);

    while (!feof(fbfile)) {
        if (fgets (line, linesize, fbfile) == NULL)
            break;
        linenr++;
        if (fb_read_strip_comment(line) == (size_t) 0) {
            /* Skip empty/comment lines */
            continue;
        }

        if (fb_parse_line (fb_cur, line, log_scale, linenr, powlim) == 0)
            continue;
        if (lim && fb_cur->p >= lim)
            break;

        /* Compute invp */
        if (fb_cur->p % 2 != 0) {
            ASSERT(sizeof(unsigned long) >= sizeof(redc_invp_t));
            fb_cur->invp = - ularith_invmod ((unsigned long) fb_cur->p);
        }

        /* To which factor base do we add? */
	size_t add_to = 0;
	if (smalllim != 0 && fb_cur->p > smalllim) {
	    add_to = nextpiece + 1;
            if (++nextpiece == nr_pieces)
                nextpiece = 0;
	}
	if (!fb_buffer_add (&fb_bufs[add_to], fb_cur)) {
	    error = 1;
	    break;
	}

        /* fb_fprint_entry (stdout, fb_cur); */
	if (fb_cur->p > maxprime)
	    maxprime = fb_cur->p;
        nr_primes++;
    }

    /* If nothing went wrong so far, put the end-of-fb markers */
    for (size_t i = 0; !error && i < nr_buffers; i++) {
        if (!fb_buffer_finish (&fb_bufs[i]))
	    error = 1;
    }
    
    /* If there was any error, free all the allocated memory */
    if (error) {
        for (size_t i = 0; i < nr_buffers; i++) {
            free(fb_bufs[i].fb);
            fb_bufs[i].fb = NULL; /* Thus caller's values will be set to NULL */
        }
    }
    if (!error && verbose)
    {
        printf ("# Factor base sucessfully read, %lu primes, largest was "
                FBPRIME_FORMAT "\n", nr_primes, maxprime);
    }

    *fb_small = fb_bufs[0].fb;
    for (size_t i = 1; i < nr_buffers; i++) {
      fb_pieces[i - 1] = fb_bufs[i].fb;
    }

    free (fb_bufs);
    fclose (fbfile);
    free (fb_cur);

    return error ? 0 : 1;
}


factorbase_degn_t *
fb_read (const char * const filename, const double log_scale, 
         const int verbose, const fbprime_t lim, const fbprime_t powlim)
{
    factorbase_degn_t *fb;
    if (!fb_read_split (&fb, NULL, filename, log_scale, 0, 0, verbose, lim, 
                        powlim))
        return NULL;
    return fb;
}


void 
fb_dump_degn (const factorbase_degn_t *fb, const char *filename)
{
    FILE *f = fopen(filename, "wb");
    if (f == NULL) {
        fprintf (stderr, "Could not open %s for writing: %s\n",
                 filename, strerror(errno));
        exit(EXIT_FAILURE);
    }
    size_t size = fb_size (fb);
    size_t written = fwrite (fb, 1, size, f);
    if (written != size) {
        fprintf (stderr, "Could not write %zu bytes to %s: %s\n",
                 size, filename, strerror(errno));
        exit(EXIT_FAILURE);
    }
    fclose(f);
}

/* Returns NULL iff the file could not be opened, and MAP_FAILED in case of
   an error (incl. when mmap() is not available) */
factorbase_degn_t *
fb_mmap(const char *fbcache MAYBE_UNUSED)
{
#ifdef HAVE_MMAP
    FILE *f = fopen(fbcache, "rb");
    if (f == NULL)
        return NULL;
    if (fseek (f, 0, SEEK_END) != 0) {
        fprintf (stderr, "fseek on %s failed: %s", fbcache, strerror(errno));
        return MAP_FAILED;
    }
    int64_t size = ftell (f);
    ASSERT_ALWAYS(size >= 0);
    if (fseek (f, 0, SEEK_SET) != 0) {
        fprintf (stderr, "fseek on %s failed: %s", fbcache, strerror(errno));
        return MAP_FAILED;
    }
    int fd = fileno (f);
    factorbase_degn_t *fb = mmap (NULL, (size_t) size, PROT_READ,
                                  MAP_PRIVATE, fd, 0);
    return fb;
#else
    return MAP_FAILED;
#endif
}

