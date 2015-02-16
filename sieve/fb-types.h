#ifndef FB_TYPES_H
#define FB_TYPES_H

/* Elementary data types for the factor base */

#include <stdint.h>
#include "las-config.h"

typedef unsigned int fbprime_t; /* 32 bits should be enough for everyone */
#define FBPRIME_FORMAT "u"
#define FBPRIME_MAX UINT_MAX
#define FBPRIME_BITS 32
typedef fbprime_t fbroot_t;
#define FBROOT_FORMAT "u"
typedef unsigned long largeprime_t; /* On IA32 they'll only get 32 bit
                                       large primes */
#define LARGEPRIME_FORMAT "lu"

#define FB_MAX_PARTS 4

/* Each slice in a factor base has a unique index */
typedef size_t slice_index_t;
/* Each factor base entry withing a slice has a unique offset */
typedef uint16_t slice_offset_t;

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

#endif
