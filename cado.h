/* Common header file for the CADO project
 
Author: Paul Zimmermann

Copyright 2007 INRIA

This file is part of the CADO project.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/

#ifndef _CADO_H
#define _CADO_H

#include "gmp.h"

/* Common asserting/debugging defines */
#include <assert.h>
#define ASSERT_ALWAYS(x) assert(x)
#ifdef WANT_ASSERT
  #define ASSERT(x) assert(x)
#else
  #define ASSERT(x)
#endif
#ifdef WANT_ASSERT_EXPENSIVE
  #define ASSERT_EXPENSIVE(x) assert(x)
#else
  #define ASSERT_EXPENSIVE(x)
#endif

/* The maximum degree of polynomials supported. Used for statically 
   allocating storage (i.e. "mpz_t poly[MAXDEGREE]") */
#define MAXDEGREE 10
/* default degrees for polynomial selection, entries are in digits */
#define DEFAULT_DEGREES {{70, 3}, {90, 4}, {ULONG_MAX, 5}}
#define DEFAULT_DEGREES_LENGTH 3

/* default rational/algebraic factor base bounds */
#define DEFAULT_RLIM {{70, 200000}, {80, 350000}, {100, 1800000}, \
                      {130, 3497867}, {ULONG_MAX, 4294967295UL}}
#define DEFAULT_RLIM_LENGTH 5
#define DEFAULT_ALIM {{70, 200000}, {80, 500000}, {100, 1800000}, \
                      {130, 11380951}, {ULONG_MAX, 4294967295UL}}
#define DEFAULT_ALIM_LENGTH DEFAULT_RLIM_LENGTH

/* default large prime bounds */
#define DEFAULT_LPBR {{70, 23}, {90, 24}, {100, 26}, {130, 28}, \
		      {ULONG_MAX, 25}}
#define DEFAULT_LPBR_LENGTH 5
#define DEFAULT_LPBA DEFAULT_LPBR
#define DEFAULT_LPBA_LENGTH DEFAULT_LPBR_LENGTH

/* default factor-residual bounds */
#define DEFAULT_MFBR {{70, 35}, {90, 37}, {100, 48}, {130, 56}, \
		      {ULONG_MAX, 39}}
#define DEFAULT_MFBR_LENGTH 5
#define DEFAULT_MFBA DEFAULT_MFBR
#define DEFAULT_MFBA_LENGTH DEFAULT_MFBR_LENGTH

/* default lambda values */
#define DEFAULT_RLAMBDA {{70, 1.5}, {90, 1.7}, {100, 2.5}, {130, 2.7}, \
                         {ULONG_MAX, 2.7}}
#define DEFAULT_RLAMBDA_LENGTH 5
#define DEFAULT_ALAMBDA DEFAULT_RLAMBDA
#define DEFAULT_ALAMBDA_LENGTH DEFAULT_RLAMBDA_LENGTH

/* default sieving block lengths */
#define DEFAULT_QINT {{70, 5000}, {90, 10000}, {100, 100000}, {130, 100000}, \
                      {ULONG_MAX, 1000000}}
#define DEFAULT_QINT_LENGTH 5

typedef struct
{
  mpz_t n;    /* number to factor */
  int degree; /* (optional) wanted degree */
} __cado_input_struct;

typedef __cado_input_struct cado_input[1];

typedef struct
{
  char name[256]; /* name */
  mpz_t n;        /* number to factor */
  double skew;    /* skewness */
  int degree;     /* (algebraic) degree */
  mpz_t *f;       /* algebraic coefficients */
  mpz_t *g;       /* rational coefficients */
  mpz_t m;        /* common root of f and g mod n */
  char type[256]; /* type (gnfs or snfs) */
  unsigned long rlim; /* rational  factor base bound */
  unsigned long alim; /* algebraic factor base bound */
  int lpbr;           /* rational  large prime bound is 2^lpbr */
  int lpba;           /* algebraic large prime bound is 2^lpba */
  int mfbr;           /* bound for rational  residuals is 2^mfbr */
  int mfba;           /* bound for algebraic residuals is 2^mfba */
  double rlambda;     /* lambda sieve parameter on the rational  side */
  double alambda;     /* lambda sieve parameter on the algebraic side */
  int qintsize;       /* sieve block range */
} __cado_poly_struct;

typedef __cado_poly_struct cado_poly[1];

/* Data types */

typedef unsigned int fbprime_t; /* 32 bits should be enough for everyone */
#define FBPRIME_FORMAT "%u"
typedef fbprime_t fbroot_t;
#define FBROOT_FORMAT "%u"
typedef unsigned long largeprime_t; /* On IA32 they'll only get 32 bit 
                                       large primes */
#define LARGEPRIME_FORMAT "%lu"

/* Factor base entry with (possibly) several roots */
typedef struct {
  fbprime_t p;            /* A prime or a prime power */
  unsigned char plog;     /* logarithm (to some suitable base) of this prime */
  unsigned char nr_roots; /* how many roots there are for this prime */
  unsigned char size;     /* The length of the struct in bytes */
  unsigned char dummy[1]; /* For dword aligning the roots */
  fbroot_t roots[0];      /* the actual length of this array is determined
                             by nr_roots */
} factorbase_degn_t;

/* Factorbase entry with prime < 2^24 and exactly 1 root */
typedef struct {
  fbprime_t p;            /* A prime or a prime power < 2^24 */
  fbroot_t root_and_log;  /* The root and approxomate log (upper 8 bits) */
} factorbase_small_t;

/* Factorbase entry with prime < 2^24 and index to next sieve array
   location where this prime divides  */
typedef struct {
  fbprime_t p;            /* A prime or a prime power < 2^24 */
  fbprime_t loc_and_log;  /* Next location in sieve array where this prime
			     will divide (low 24 bits) and approxomate log
			     (upper 8 bits) */
} factorbase_small_inited_t;


typedef struct {
  factorbase_degn_t *fullfb; /* The complete factor base */
  factorbase_small_t *fbL1, *fbL2;
  factorbase_small_inited_t *fbL1init, *fbL2init;
  factorbase_degn_t *fblarge; /* Pointer to the entries in fullfb with primes
  				 > fbL2bound */
  unsigned int fbL1size;  /* Number of entries in fbL1, incl. stop marker */
  unsigned int fbL1bound; /* Upper bound on primes in fbL1 */
  unsigned int fbL2size;  /* Number of entries in fbL2 */
  unsigned int fbL2bound; /* Upper bound on primes in fbL2 */
} factorbase_t[1];

typedef struct {
  long a;		/* only a is allowed to be negative */
  unsigned long b;
  int nb_rp;		/* number of rational primes */
  int nb_ap;		/* number of algebraic primes */
  unsigned long * rp;	/* array of rational primes */
  unsigned long * ap;	/* array of algebraic primes */
  unsigned long * ar;	/* array of corresponding root (optional, this can
			   be garbage) */
} relation_t;
#endif
