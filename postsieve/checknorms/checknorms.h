/**
 * \file    checknorms.h
 * \author  Jerome Milan (heavily based on Paul Zimmermann's checknorms program)
 * \date    Fri Jan 25 2008
 * \brief   Check relations produced by siever and factor residues.
 * \version 0.1
 *
 * This is a slight modification of Paul Zimmermann's checknorms program
 * to use TIFA's factoring primitives instead of ECM. Usage is unchanged:
 *
 * <tt>checknorms -poly c80.poly c80.rels1 [c80.rels2 ...]</tt>
 *
 * Some options have been added. For more information, type
 * <tt>checknorms -h</tt>.
 */

#if !defined(_CADO_POSTSIEVE_CHECKNORMS_H_)
#define _CADO_POSTSIEVE_CHECKNORMS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "cado.h"   // For cado_poly type
#include "tifa.h"   // For TIFA's *_array_t types.

//------------------------------------------------------------------------------
#define VERBOSE             1    // to print warnings
#define PERFORM_CHECK_PRIME 1    // to check size of prime
#define DEBUG               0    // to run step by step (with lots of messages!)
#define PRINT_SUMMARY       1    // to print reporting info at end of program
//------------------------------------------------------------------------------
#define NORM_IS_ONE           0  //
#define NORM_IS_ACCEPTED      1  // Error codes used by the
#define NORM_IS_REJECTED      2  // factor_completely function
#define ERROR_FACTORING_NORM  3  //
//-----------------------------------------------------------------------------
#define NFACTORS        8        // (initial) size of array for norm's factors
#define MAXNMSG         10       // max number of warning messages
#define REL_MAXLENGTH   1024     // max length in characters of a relation
#define NMILLER_RABIN   20       // number of iterations in composition test
#define DFLT_MAX_NLP    4        // default number of large primes accepted
#define DFLT_TDMAX      100      // default bound for trial division
//------------------------------------------------------------------------------
#define IS_PRIME(X)     (0 != mpz_probab_prime_p((X), NMILLER_RABIN))
#define IS_COMPOSITE(X) (0 == mpz_probab_prime_p((X), NMILLER_RABIN))
#define IS_SQUARE(X)    (0 != mpz_perfect_square_p((X)))
#define BITSIZE(X)      (mpz_sizeinbase((X), 2))
//------------------------------------------------------------------------------
#if DEBUG
    //
    // _WARNING_: Setting DEBUG to one will prompt for user input before
    //            reading a relation and print a truckload of messages.
    //
    #define MSGDEBUG(...) gmp_printf(__VA_ARGS__);fflush(stdout);
#else
    #define MSGDEBUG(...) /* voluntarily left empty */
#endif
#if VERBOSE
    #define WARNING(...) gmp_fprintf(stderr, "WARNING: " __VA_ARGS__)
#else
    #define WARNING(...) /* voluntarily left empty */
#endif
#if PERFORM_CHECK_PRIME
    #define CHECK_PRIME(P, SB, A, B) check_prime((P), (SB), (A), (B))
#else
    #define CHECK_PRIME(P, SB, A, B) /* voluntarily left empty */
#endif
#define ERROR(...)  gmp_fprintf(stderr, "ERROR: " __VA_ARGS__)
#define MSG(...)    gmp_fprintf(stderr, __VA_ARGS__)
//------------------------------------------------------------------------------
//
// Global variables (bouh!) used for reporting... 
//
static unsigned long count_too_many_factors     = 0;
static unsigned long count_too_large_prime_norm = 0;
static unsigned long count_too_large_factor     = 0;
static unsigned long count_no_factor_found      = 0;
static unsigned long count_partial_fact_found   = 0;
static unsigned long count_perfect_power        = 0;
static unsigned long count_check_prime          = 0;
static unsigned long count_large_algprime       = 0;
static unsigned long count_large_ratprime       = 0;
static unsigned long count_not_coprime          = 0;
static unsigned long count_not_dividing         = 0;
static unsigned long count_larger_mfba          = 0;
static unsigned long count_larger_mfbr          = 0;
static unsigned long total_nrelations_read      = 0;
//-----------------------------------------------------------------------------
#define NPIX      46        // number of entries in pi(x) table
#define MAX_X_PIX 1000000   // largest x for pi(x) table
//-----------------------------------------------------------------------------
//
// Approximate values of pi(x) for small x.
// Values taken from: http://www.trnicely.net/pi/tabpi.html.
//
static const uint64_t pi_table_x[NPIX] = {
        10,     20,     30,     40,     50,     60,     70,     80,     90,
       100,    200,    300,    400,    500,    600,    700,    800,    900,
      1000,   2000,   3000,   4000,   5000,   6000,   7000,   8000,   9000,
     10000,  20000,  30000,  40000,  50000,  60000,  70000,  80000,  90000,
    100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000,
   1000000
};
//-----------------------------------------------------------------------------
static const uint32_t pi_table_pi[NPIX] = {
        4,     8,    10,    12,    15,    17,    19,    22,    24,
       25,    46,    62,    78,    95,   109,   125,   139,   154, 
      168,   303,   430,   550,   669,   783,   900,  1007,  1117,
     1229,  2262,  3245,  4203,  5133,  6057,  6935,  7837,  8713,
     9592, 17984, 25997, 33860, 41538, 49098, 56543, 63951, 71274,
    78498
};
//-----------------------------------------------------------------------------
//
// Returns an approximation of pi(x) for x <= MAX_X_PIX.
// Returned approximation is always greater than or equals to pi(x).
//
uint32_t pi_approx(uint64_t x);
//-----------------------------------------------------------------------------
//
// Set 'norm' to the norm of g(a,b), i.e:
//   norm <- |g(a,b)| = abs(g[1]*a + g[0]*b)
//
void get_rational_norm(mpz_t norm, mpz_t *g, long a, unsigned long b);
//-----------------------------------------------------------------------------
//
// Set 'norm' to the norm of f(a,b), i.e:
//   norm <- |f(a,b)| = abs(f[d]*a^d + ... + f[0]*b^d)
//
void get_algebraic_norm(mpz_t norm, mpz_t *f, int d, long a, unsigned long b);
//-----------------------------------------------------------------------------
//
// Prints a warning message if p <= sb
//
inline void check_prime(mpz_t p, unsigned long sb, long a, unsigned long b);
//-----------------------------------------------------------------------------
//
// Attempts to factor norm and stores the found factors in the array
// factors and the multiplicities in the array multis. It returns:
//
//          NORM_IS_ONE: if norm == 1
//     NORM_IS_ACCEPTED: if norm has exactly or less than maxnlp prime
//                       factors and if each of these prime are less
//                       than 2^lp
//     NORM_IS_REJECTED: if norm has more than maxnlp prime or if one of
//                       its prime factor is greater than (or equal to) 2^lp
// ERROR_FACTORING_NORM: if norm could not be factored completely or if
//                       other apocalyptic errors occur.
//
unsigned long factor_completely(
    mpz_array_t*    const factors,
    uint32_array_t* const multis,
    const mpz_t norm,
    const size_t lp,
    const unsigned int maxnlp,
    const unsigned int cmult
);
//-----------------------------------------------------------------------------
//
// Returns true iff a and b are coprime.
//
inline bool coprime(long a, unsigned long b);
//-----------------------------------------------------------------------------
//
// Checks, filters, completes relations read from a relation file and prints
// surviving relations on stdout.
//
//       f : name of relation file
//   cpoly : polynomial used
// verbose : verbosity switch
//  maxnlp : maximum number of large prime accepted in norms' residuals.
//   cmult : if non-null, considers p^m as m large primes instead of just one. 
//    mfbr : bound for rational  residues is 2^mfbr
//    mfba : bound for algebraic residues is 2^mfba
//  primes : small primes used for trial division (should contains primes
//           omitted in relations)
//     npa : number of small primes to trial divide by on rational side
//     npa : number of small primes to trial divide by on algebraic side
//
unsigned long checkrels(
    char *f,
    cado_poly cpoly,
    int verbose,
    unsigned int maxnlp,
    unsigned int cmult,
    size_t mfbr,
    size_t mfba,
    unsigned long *primes,
    unsigned long npr,
    unsigned long npa
);
//-----------------------------------------------------------------------------
//
// Prints usage message on command line.
//
void print_usage(char* progname);
//-----------------------------------------------------------------------------
//
// Prints help message if -h option passed on command line.
//
void print_help(char* progname);
//-----------------------------------------------------------------------------
// ckn_mu_secs: ChecKNorms's MicroSECondS
//
// Returns current time expressed as microseconds. Used to get elapsed time
// intervals taking into account CPU time _and_ system time spent on behalf of
// the calling process. microseconds() from src/utils/timing.c discards system
// time.
//
uint64_t ckn_mu_secs();
//-----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

#endif

