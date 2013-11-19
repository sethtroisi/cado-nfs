/*****************************************************************
*                Functions for the factor base                  *
*****************************************************************/

#ifndef FB_H
#define FB_H

#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#ifdef HAVE_SYS_MMAN_H
#include <sys/mman.h>
#else
#define MAP_FAILED ((void *) -1)
#endif
#include <gmp.h>
#include "las-config.h"

/* Data types */

typedef unsigned int fbprime_t; /* 32 bits should be enough for everyone */
#define FBPRIME_FORMAT "%u"
#define FBPRIME_MAX UINT_MAX
#define FBPRIME_BITS 32
typedef fbprime_t fbroot_t;
#define FBROOT_FORMAT "%u"
typedef unsigned long largeprime_t; /* On IA32 they'll only get 32 bit 
                                       large primes */
#define LARGEPRIME_FORMAT "%lu"

/* If SUPPORT_LARGE_Q is defined, 64-bit redc is used in the function that
   converts roots to the p-lattice, and the redc code needs a 64-bit 
   precomputed inverse. If SUPPORT_LARGE_Q is not defined, we store only a
   32-bit inverse to conserve memory. */
#if defined(SUPPORT_LARGE_Q)
typedef uint64_t redc_invp_t;
#else
typedef uint32_t redc_invp_t;
#endif

/* The following format takes 16+4k bytes per prime with k roots. Since
 * the expected number of roots for primes with at least one root is
 * 1.58 (for generic Galois group), we are slightly above 14 bytes per
 * root.
 */

/* Factor base entry with (possibly) several roots */
typedef struct {
  fbprime_t p;            /* A prime or a prime power */
  unsigned char nr_roots; /* how many roots there are for this prime */
  unsigned char plog;     /* logarithm (to some suitable base) of this prime */
  unsigned char size;     /* The length of the struct in bytes */
  unsigned char dummy[1]; /* For dword aligning the roots. In fact uneeded, C99
                             guarantees proper alignment of roots[]. It's only a
                             precaution against -fpack-struct or other
                             ABI-breaking behaviours */
  redc_invp_t invp;       /* -1/p (mod 2^wordsize) for REDC */
  /* Note that invp may have a stronger alignment constraint than p, thus must
   * not appear before the tiny fields plog and nr_roots which can easily
   * fit inbetween the two */
  fbroot_t roots[0];      /* the actual length of this array is determined
                             by nr_roots */
} factorbase_degn_t;

#define FB_END ((fbprime_t) 1)

void		fb_fprint_entry (FILE *, const factorbase_degn_t *);
void            fb_fprint (FILE *, const factorbase_degn_t *);
void            fb_sortprimes (fbprime_t *, const unsigned int);
unsigned char	fb_log (double, double, double);
fbprime_t       fb_is_power (fbprime_t);
int             fb_make_linear (factorbase_degn_t **, factorbase_degn_t ***, 
                                const mpz_t *, fbprime_t, fbprime_t, int, 
                                fbprime_t, double, int, int, FILE *);
int             fb_read_split (factorbase_degn_t **, factorbase_degn_t ***, 
                               const char *, double, fbprime_t, int, int, 
                               fbprime_t, fbprime_t);
fbprime_t	*fb_extract_bycost (const factorbase_degn_t *, 
                                    const fbprime_t, const fbprime_t costlim);
size_t          fb_size (const factorbase_degn_t *);                   
size_t          fb_nroots_total (const factorbase_degn_t *fb);
void            fb_dump_degn (const factorbase_degn_t *, const char *);
factorbase_degn_t *	fb_mmap(const char *);

/* Some inlined functions which need to be fast */
  
/* Hack to get around C's automatic multiplying constants to be added to 
   pointers by the pointer's base data type size */

__attribute__ ((unused))
static inline factorbase_degn_t *
fb_skip (const factorbase_degn_t *fb, const size_t s)
{
  return (factorbase_degn_t *)((char *)fb + s);
}
  
__attribute__ ((unused))
static inline factorbase_degn_t *
fb_next (const factorbase_degn_t *fb)
{
  return (factorbase_degn_t *)((char *)fb + fb->size);
}

__attribute__ ((unused))
static size_t
fb_entrysize (const factorbase_degn_t *fb)
{
  return (sizeof (factorbase_degn_t) + fb->nr_roots * sizeof (fbroot_t));
}

__attribute__ ((unused))
static size_t
fb_entrysize_uc (const unsigned char n)
{
  return (sizeof (factorbase_degn_t) + n * sizeof (fbroot_t));
}

/* Write the end-of-factor-base marker at *fb */
__attribute__ ((unused))
static void
fb_write_end_marker (factorbase_degn_t *fb)
{
  fb->p = FB_END;
  fb->invp = -(redc_invp_t)1;
  fb->nr_roots = 0;
}

/* Most often, we're in fact considering const iterators, but we refuse
 * to bother with having two interfaces */

struct fb_iterator_s {
    factorbase_degn_t * fb;
    int i;
};

typedef struct fb_iterator_s fb_iterator[1];
typedef struct fb_iterator_s *fb_iterator_ptr;
typedef const struct fb_iterator_s *fb_iterator_srcptr;

static inline void fb_iterator_init_set_fb(fb_iterator_ptr t, factorbase_degn_t * fb)
{
    memset(t, 0, sizeof(*t));
    t->fb = fb;
}

static inline void fb_iterator_init_set(fb_iterator_ptr t, fb_iterator_srcptr u)
{
    memcpy(t, u, sizeof(*u));
}

static inline void fb_iterator_clear(fb_iterator_ptr t)
{
    memset(t, 0, sizeof(*t));
}

static inline void fb_iterator_next(fb_iterator_ptr t)
{
    if (++(t->i) < t->fb->nr_roots)
        return;
    t->fb = fb_next(t->fb);
    t->i = 0;
}

static inline fbprime_t fb_iterator_get_r(fb_iterator_srcptr t)
{
    return t->fb->roots[t->i];
}

static inline int fb_iterator_lessthan_fb(fb_iterator_srcptr t, const factorbase_degn_t * fb)
{
    return (char*)(t->fb) < (char*)fb;
}

static inline int fb_iterator_lessthan(fb_iterator_srcptr t, fb_iterator_srcptr u)
{
    int r = (char*)(u->fb) - (char*)(t->fb);
    if (r > 0) return 1;
    if (r < 0) return 0;
    return t->i < u->i;
}

/* Computes t - u. Assumes for the primary branch that u <= t. If t < u,
 * symmetrises by swapping arguments */
static inline int fb_iterator_diff(fb_iterator_srcptr t, fb_iterator_srcptr u)
{
    fb_iterator q;
    if (fb_iterator_lessthan(t, u)) {
        return -fb_iterator_diff(u, t);
    }
    fb_iterator_init_set(q, u);
    int n = -q->i;
    q->i = 0;
    for( ; fb_iterator_lessthan_fb(q, t->fb) ; ) {
        n += q->fb->nr_roots;
        q->fb = fb_next(q->fb);
    }
    n += t->i;
    fb_iterator_clear(q);
    return n;
}
static inline int fb_iterator_diff_fb(fb_iterator_srcptr t, factorbase_degn_t * u)
{
    fb_iterator qu;
    fb_iterator_init_set_fb(qu, u);
    int n = fb_iterator_diff(t, qu);
    fb_iterator_clear(qu);
    return n;
}

static inline int fb_diff(factorbase_degn_t * t, factorbase_degn_t * u)
{
    fb_iterator qu;
    fb_iterator_init_set_fb(qu, u);
    fb_iterator qt;
    fb_iterator_init_set_fb(qt, t);
    int n = fb_iterator_diff(qt, qu);
    fb_iterator_clear(qt);
    fb_iterator_clear(qu);
    return n;
}

static inline int fb_iterator_over(fb_iterator_srcptr t)
{
    return t->fb->p == FB_END;
}

#endif
