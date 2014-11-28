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
#define FBPRIME_FORMAT "u"
#define FBPRIME_MAX UINT_MAX
#define FBPRIME_BITS 32
typedef fbprime_t fbroot_t;
#define FBROOT_FORMAT "u"
typedef unsigned long largeprime_t; /* On IA32 they'll only get 32 bit 
                                       large primes */
#define LARGEPRIME_FORMAT "lu"

/* If SUPPORT_LARGE_Q is defined, 64-bit redc is used in the function that
   converts roots to the p-lattice, and the redc code needs a 64-bit 
   precomputed inverse. If SUPPORT_LARGE_Q is not defined, we store only a
   32-bit inverse to conserve memory. */
#if defined(SUPPORT_LARGE_Q)
typedef uint64_t redc_invp_t;
#else
typedef uint32_t redc_invp_t;
#endif

#define LOG_SCALE 1.4426950408889634 /* 1/log(2) to 17 digits, rounded to
                                        nearest. This is enough to uniquely
                                        identify the corresponding IEEE 754
                                        double precision number */

/* The following format takes 16+4k bytes per prime with k roots. Since
 * the expected number of roots for primes with at least one root is
 * 1.58 (for generic Galois group), we are slightly above 14 bytes per
 * root.
 */

unsigned char	fb_log (double, double, double);
fbprime_t       fb_pow (fbprime_t, unsigned long);
fbprime_t       fb_is_power (fbprime_t, unsigned long *);
unsigned char   fb_make_steps(fbprime_t *, fbprime_t, double);

/* Some inlined functions which need to be fast */
  
#endif
