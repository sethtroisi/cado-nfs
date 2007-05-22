/*****************************************************************
*                Functions for the factor base                  *
*****************************************************************/

#include "cado.h"

void		fb_print_entry (factorbase_t *);
factorbase_t *	fb_add_to (factorbase_t *, size_t *, size_t *, const size_t, 
                           factorbase_t *);
factorbase_t *	fb_find_p (factorbase_t *, const fbprime_t);
unsigned char	fb_log (double, double, double);
factorbase_t * 	fb_make_linear (mpz_t *, const fbprime_t, const double, 
                                const int);
factorbase_t *	fb_read (const char *, const double, const int);
void		fb_disable_roots (factorbase_t *, const unsigned long, 
                                  const int);
void		fb_restore_roots (factorbase_t *, const unsigned long, 
                                  const int);
fbprime_t       fb_maxprime (factorbase_t *);

/* Some inlined functions which need to be fast */
  
/* Hack to get around C's automatic multiplying constants to be added to 
   pointers by the pointer's base data type size */

__attribute__ ((unused))
static inline factorbase_t *
fb_skip (const factorbase_t *fb, const size_t s)
{
  return (factorbase_t *)((char *)fb + s);
}
  
__attribute__ ((unused))
static inline factorbase_t *
fb_next (const factorbase_t *fb)
{
  return (factorbase_t *)((char *)fb + fb->size);
}

__attribute__ ((unused))
static size_t
fb_entrysize (const factorbase_t *fb)
{
  return (sizeof (factorbase_t) + fb->nr_roots * sizeof (fbroot_t));
}

__attribute__ ((unused))
static size_t
fb_entrysize_uc (const unsigned char n)
{
  return (sizeof (factorbase_t) + n * sizeof (fbroot_t));
}
