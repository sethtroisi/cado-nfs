/*****************************************************************
*                Functions for the factor base                  *
*****************************************************************/

#ifndef __FB_H
#define __FB_H

#include "config.h"
#include "cado.h"

#define FB_END ((fbprime_t) 1)

void		fb_print_entry (factorbase_degn_t *);
void            fb_sortprimes (fbprime_t *, const unsigned int);
unsigned char	fb_log (double, double, double);
void            fb_init_firstlog (factorbase_t);
factorbase_degn_t * 	fb_make_linear (mpz_t *, const fbprime_t, const double, 
                                const int);
factorbase_degn_t *	fb_read (const char *, const double, const int);
void		fb_disable_roots (factorbase_degn_t *, const unsigned long, 
                                  const int);
void		fb_restore_roots (factorbase_degn_t *, const unsigned long, 
                                  const int);
void 		fb_extract_small (factorbase_t, const unsigned int, const int,
                                  const int);
int             fb_check (factorbase_t, cado_poly, int);
void            fb_clear (factorbase_t);
fbprime_t	*fb_extract_bycost (const factorbase_degn_t *, 
                                    const fbprime_t, const fbprime_t costlim);
                   

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

#endif
