Create a config.h file with contents like


/*************************************************/
#ifndef __CONFIG_H__
#define __CONFIG_H__

#define HAVE_MSRH 1
#define WANT_ASSERT 1
#define WANT_ASSERT_EXPENSIVE 0

/* Number of blocking levels for small factor base primes, should
   correspond to cache levels. Sieving will be done in SIEVE_BLOCKING + 1
   passes: SIEVE_BLOCKING passes updating directly, and one pass with
   bucket sorting. Bucket sorting not implementest atm. */
#define SIEVE_BLOCKING 1
static const unsigned long CACHESIZES[SIEVE_BLOCKING] = {16384};

/* Controls whether the code is allowed to skip forward in the factor
   base during trial division. Mostly meant for speed tests/debugging */
#define TRIALDIV_SKIPFORWARD 1

/* Define this to the a value of the relation you want to trace */
/* #define TRACE_RELATION_A -1234L */

#endif
/*************************************************/


Compile with "make".

Make a factor base with, for example,

echo "rootfind(5017309194362523*x^4 -1406293661386525*x^3 -1131155401311965*x^2 +4737694118287353*x -3415040824020545, 1000000)" | gp -q -p 1000000 rootfind.gp >c80.roots.txt

Run the siever with

./sieve -poly c80.poly -fb c80.roots.txt -1000000 1000000 2001 2101

where c80.poly is the polynomial file.

Parameters for the siever are
-v : verbose output, including timing info
-reports_a_len : how many entries to allocate for reports list on 
                 algebraic side
-reports_r_len : how many entries to allocate for reports list on 
                 rational side
