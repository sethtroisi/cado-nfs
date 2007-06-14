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
   bucket sorting. */
#define SIEVE_BLOCKING 1
static const unsigned long CACHESIZES[SIEVE_BLOCKING] = {16384};

#endif
/*************************************************/


Compile with

gcc -I ../ -O0 -g -o sieve -Wall -Wextra sieve.c fb.c -lm -lgmp

(Note: linking GMP statically with /usr/lib/libgmp.a for example may yield a
speedup up to 7%.)

or, if you're like me and have compiled a static GMP with full frame pointer 
support for debugging on 32 bit machines, use i.e.

gcc -I ../ -O0 -g -o sieve -Wall -Wextra sieve.c -lm /usr/local/lib/libgmp_fp.a

Make a factor base with, for example,

echo "rootfind(5017309194362523*x^4 -1406293661386525*x^3 -1131155401311965*x^2 +4737694118287353*x -3415040824020545, 1000000)" | gp -q -p 1000000 rootfind.gp >c80.roots.txt

Run the siever with

./sieve -poly c80.poly -fb c80.roots.txt -1000000 1000000 2001 2101

where c80.poly is the polynomial file.
