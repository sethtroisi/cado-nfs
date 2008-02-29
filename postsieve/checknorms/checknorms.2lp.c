/**
 * \file    checknorms.c
 * \author  Jerome Milan (heavily based on Paul Zimmermann's checknorms program)
 * \date    Tue Jan 8 2008
 * \brief   Check relations produced by siever and factor residues.
 * \version 0.2.1
 *
 * This is a slight modification of Paul Zimmermann's checknorms program
 * to use TIFA's factoring primitives instead of ECM. Usage is unchanged:
 *
 * <tt>checknorms -poly c80.poly c80.rels1 [c80.rels2 ...]</tt>
 */

/*
 * History:
 * --------
 *   0.2.1: Tue Jan 22 2008 by JM
 *          - Jump to next relation if a prime does not divide norm (was
 *            previously aborting the whole program)
 *          - Fixed extra comma appearing when appending primes to empty list.
 *          - Fixed compilation warning (hopefully).
 *          - Fixed bug in control flow.
 *          - Raised tdmax to 20 when using option '-factorall'
 *     0.2: Tue Jan  8 2008 by JM
 *          - Added '-factorall' option.
 *     0.1: Thu Dec 20 2007 by JM
 *          - Initial version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/time.h>     /* for getrusage */
#include <sys/resource.h> /* for getrusage */

#include "cado.h"
#include "utils.h"
#include "tifa.h"
#include "funcs.h"        // "Private" header file from the TIFA library

//------------------------------------------------------------------------------
#define VERBOSE             1  // to print warnings
#define PERFORM_CHECK_PRIME 1  // to check size of prime
#define DEBUG               0  // to run step by step (with lots of messages!)
#define PRINT_SUMMARY       1  // to print reporting info at end of program
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
#define ERROR(...)   gmp_fprintf(stderr, "ERROR:   " __VA_ARGS__)
#define MSG(...)     gmp_fprintf(stderr, __VA_ARGS__)
//------------------------------------------------------------------------------
#define NFACTORS        8    /* (initial) size of array for norm's factors */
#define MAXNMSG         10   /* max number of warning messages             */
#define REL_MAXLENGTH   1024 /* max length in characters of a relation     */
#define NMILLER_RABIN   20   /* number of iterations in composition test   */
//------------------------------------------------------------------------------
#define IS_PRIME(X)     (0 != mpz_probab_prime_p((X), NMILLER_RABIN))
#define IS_COMPOSITE(X) (0 == mpz_probab_prime_p((X), NMILLER_RABIN))
#define IS_SQUARE(X)    (0 != mpz_perfect_square_p((X)))
#define BITSIZE(X)      (mpz_sizeinbase((X), 2))
//------------------------------------------------------------------------------
//
// Global variables (bouh!) used for error reporting...
//
static unsigned long count_more_two_factors     = 0;
static unsigned long count_too_large_norm       = 0;
static unsigned long count_bad_factor           = 0;
static unsigned long count_too_large_prime_norm = 0;
static unsigned long count_no_factor_found      = 0;
static unsigned long count_partial_fact_found   = 0;
static unsigned long count_factorization_pb     = 0;
static unsigned long count_perfect_square       = 0;
static unsigned long count_perfect_power        = 0;
static unsigned long check_prime_reports        = 0;
//-----------------------------------------------------------------------------
uint64_t ckn_mu_secs() {
    // ckn_mu_secs: ChecKNorms's MicroSECondS
    //
    // Returns current time expressed as microseconds. Used to get elapsed time
    // intervals taking into account CPU time _and_ system time spent on behalf
    // of the calling process. microseconds() from src/utils/timing.c discards 
    // system time.
    //
    struct rusage res[1];
    getrusage(RUSAGE_SELF, res);
    
    uint64_t r;
    r  = (uint64_t) res->ru_utime.tv_sec + res->ru_stime.tv_sec;
    r *= (uint64_t) 1000000UL;
    r += (uint64_t) res->ru_utime.tv_usec + res->ru_stime.tv_usec;
    return r;
}
//-----------------------------------------------------------------------------
void get_rational_norm(mpz_t norm, mpz_t *g, long a, unsigned long b) {
    //
    // Set 'norm' to the norm of g(a,b), i.e:
    //   norm <- |g(a,b)| = abs(g[1]*a + g[0]*b)
    //
    mpz_mul_si(norm, g[1], a);
    mpz_addmul_ui(norm, g[0], b);
    mpz_abs(norm, norm);
}
//-----------------------------------------------------------------------------
void get_algebraic_norm(mpz_t norm, mpz_t *f, int d, long a, unsigned long b) {
    //
    // Set 'norm' to the norm of f(a,b), i.e:
    //   norm <- |f(a,b)| = abs(f[d]*a^d + ... + f[0]*b^d)
    //
    mpz_t B; /* powers of b */
    mpz_init_set_ui(B, 1);
    mpz_set(norm, f[d]);

    while (d-- > 0) {
        mpz_mul_si(norm, norm, a);
        mpz_mul_ui(B, B, b);
        mpz_addmul(norm, f[d], B);
    }
    mpz_clear(B);
    mpz_abs(norm, norm);
}
//-----------------------------------------------------------------------------
inline void check_prime(mpz_t p, unsigned long sb, long a, unsigned long b) {
    //
    // Prints a warning message if p <= sb
    //
    // _WARNING_: Uses and modifies the following global variable:
    //            - check_prime_reports
    //
    if (mpz_cmp_ui(p, sb) <= 0) {
        check_prime_reports++;
        if (check_prime_reports <= MAXNMSG) {
            WARNING("(a=%ld, b=%lu): prime %Zd smaller "
                    "than factor base bound %lu\n", a, b, p, sb);
        }
    }
}
//-----------------------------------------------------------------------------
unsigned long factor(mpz_t p1, mpz_t p2,
                     unsigned int* const m1 __attribute__ ((unused)),
                     unsigned int* const m2 __attribute__ ((unused)),
                     const mpz_t norm,
                     const size_t lp) {

    MSGDEBUG("factor: entered\n");
    //
    // Factor norm into at most 2 primes smaller than 2^lp:
    //
    //     0) if norm = 1, return 0
    //     1) if norm is prime, norm < 2^lp, put it into p1 and return 1
    //     2) if norm = p1*p2 with p1,p2 primes < 2^lp, return 2
    //     3) otherwise if norm has 3 prime factors of more,
    //        or if one prime factor is >= 2^lp, return 3.
    //
    // _TO_DO_: Check that the rare case p1 == p2 doesn't destroy anything.
    //
    // _WARNING_: Uses and modifies the following global variables:
    //            - count_more_two_factors
    //            - count_bad_factor
    //            - count_factorization_pb
    //            - count_too_large_prime_norm
    //            - count_too_large_norm
    //            - count_no_factor_found
    //            - count_perfect_square
    //
    unsigned long sn = BITSIZE(norm);

    if (mpz_cmp_ui(norm, 1) == 0) {
        MSGDEBUG("norm is one!\n", norm);
        return 0;
    }
    if (IS_PRIME(norm)) {
        if (sn > lp) {
            //
            // This happens very frequently
            //
            count_too_large_prime_norm++;
            MSGDEBUG("Prime norm %Zd is too large\n", norm);
            return 3;
        }
        mpz_set(p1, norm);
        return 1;
    }
    //
    // Now, we know that norm has at least 2 prime factors
    //
    if (sn > 2 * lp) {
        //
        // norm has either exactly 2 prime factors with (at least) one larger
        // than 2^lp or norm has more than 2 factors. Either way, we discard the
        // relation...
        //
        count_too_large_norm++;
        if (count_too_large_norm < MAXNMSG) {
            WARNING("Composite norm %Zd is too large\n", norm);
        }
        return 3;
    }
    if (IS_SQUARE(norm)) {
        //
        // _TO_DO_: Check that the rare case p1 == p2 doesn't destroy anything.
        //
        count_perfect_square++;
        if (count_perfect_square < MAXNMSG) {
            WARNING("norm is a perfect square!\n");
        }
        mpz_sqrt(p1, norm);
        if (IS_PRIME(p1)) {
            mpz_set(p2, p1);
            return 2;
        }
        count_more_two_factors++;
        return 3;
    }
    unsigned long int retval = 3;

    //
    // Use the TIFA library to factor norm
    //
    mpz_array_t* const norm_factors = alloc_mpz_array(NFACTORS);

    MSGDEBUG("factor: tifa_factor() called with norm %Zd\n", norm);
    //
    // _TO_DO_: We could use TIFA's FIND_COMPLETE_FACTORIZATION mode
    //          but a) it is not fully tested / debugged yet and I'm not
    //          comfortable using it and b) it is overkill here.
    //
    ecode_t ecode = tifa_factor(norm_factors, NULL, norm, FIND_SOME_FACTORS);

    if (ecode != SOME_FACTORS_FOUND) {
        //
        // Should be rare but I cannot give any warranties here... So what
        // should we do? Rerun tifa_factor in another mode? Or just give up?
        // For the time being, we just discard the relation...
        //
        count_no_factor_found++;
        if (count_no_factor_found < MAXNMSG) {
            WARNING("TIFA found no factor!\n", norm);
        }
        goto clear_and_return;
    }
    //
    // We found non trivial factors of norm. If norm = p1 * p2 with pi prime,
    // it has obviously only 2 non trivial factors: p1 and p2.
    //
    if (norm_factors->length > 2) {
        //
        // We found more than 2 different non trivial factors of norm. They may
        // not be prime but that does not make any difference: norm != p1 * p2.
        //
        MSGDEBUG("Composite norm %Zd has more than 2 prime factors\n", norm);
        count_more_two_factors++;
        goto clear_and_return;
    }
    mpz_ptr x1 = norm_factors->data[0];

    if (norm_factors->length == 1) {
        //
        // If the found non trivial factor is not prime there is no way norm
        // could be a product of only two primes.
        //
        if (IS_PRIME(x1) && (BITSIZE(x1) <= lp)) {
            mpz_divexact(p2, norm, x1);
            if (IS_PRIME(p2) && (BITSIZE(p2) <= lp)) {
                //
                // norm = p1 * p2 with pi < 2^lp
                //
                mpz_set(p1, x1);
                retval = 2;
            } else {
                count_bad_factor++;
                MSGDEBUG("factor p2=%Zd not prime and/or too large!\n", p2);
            }
        } else {
            count_bad_factor++;
            MSGDEBUG("factor p1=%Zd not prime and/or too large!\n", p1);
        }
        goto clear_and_return;
    }
    mpz_ptr x2 = norm_factors->data[1];
    //
    // We now have exactly two different (since we took care of the square case
    // before) non trivial factors of norm (x1 and x2) but they may not be
    // prime.
    // If norm = p1 * p2 then: (x1, x2) = (p1, p2) or (x1, x2) = (p2, p1).
    // If one of the xi is not prime, then norm has more than 2 prime factors.
    // If both xi are prime we should check that x1 * x2 = norm.
    //
    if (IS_PRIME(x1) && IS_PRIME(x2)) {
        unsigned long int sx1 = BITSIZE(x1);
        unsigned long int sx2 = BITSIZE(x2);

        bool xok = (sx1 <= lp) && (sx2 <= lp);
        bool nok = (sx1 + sx2 == sn) || (sx1 + sx2 == sn + 1);

        if (xok && nok) {
            //
            // We have norm = p1 * p2 with pi < lp
            //
            // _WARNING_: Since we don't compute the actual product p1 * p2
            //            but just check the sizes of pi and norm, this may be
            //            wrong if norm is divisible by 2 (which I assume will
            //            never be the case since we divided it with small
            //            primes).
            //
            mpz_set(p1, x1);
            mpz_set(p2, x2);
            retval = 2;
        } else {
            count_bad_factor++;
            MSGDEBUG("factors %Zd and %Zd are prime but sizes mismatch!"
                     " (norm=%Zd)\n", x1, x2, norm);
            MSGDEBUG("size of %Zd is %lu\n", x1, sx1);
            MSGDEBUG("size of %Zd is %lu\n", x2, sx2);
            MSGDEBUG("size of %Zd is %lu\n", norm, sn);
            MSGDEBUG("lp is %lu\n", lp);
        }
    } else {
        //
        // See comment about the FIND_COMPLETE_FACTORIZATION mode above.
        // This mode is used by the function factor_complete which is used
        // when the '-factorall' option is passed on the command line...
        //
        count_more_two_factors++;
        MSGDEBUG("factor %Zd and / or %Zd is not prime\n", x1, x2);
    }

  clear_and_return:

    clear_mpz_array(norm_factors);
    free(norm_factors);

    return retval;

}
//-----------------------------------------------------------------------------
unsigned long factor_completely(mpz_t p1, mpz_t p2,
                                unsigned int* const m1,
                                unsigned int* const m2,
                                const mpz_t norm,
                                const size_t lp) {
    MSGDEBUG("factor_completely: entered\n");
    //
    // _WARNING_: This is NOT the 'standard' behaviour of checknorms as we
    //            knew it!
    //
    // Factor norm into at most 2 primes smaller than 2^lp:
    //
    //     0) if norm = 1, return 0
    //     1) if norm = p^m, with p a prime less than 2^lp, put it into p1,
    //        set *m1 to m and return 1
    //     2) if norm = pa^ma * pn^mb with pa, pb primes < 2^lp, set p1 to pa,
    //        p2 to pb, *m1 to ma, *m2 to mb and return 2
    //     3) otherwise if norm has 3 prime factors of more,
    //        or if one prime factor is >= 2^lp, return 3.
    //
    // _TO_DO_: Check that the rare case p1 == p2 doesn't destroy anything.
    //
    // _WARNING_: Uses and modifies the following global variables:
    //            - count_too_large_prime_norm
    //            - count_bad_factor
    //            - count_no_factor_found
    //            - count_perfect_power
    //            - count_more_two_factors
    //
    unsigned long sn = BITSIZE(norm);

    if (mpz_cmp_ui(norm, 1) == 0) {
        MSGDEBUG("norm is one!\n", norm);
        return 0;
    }
    if (IS_PRIME(norm)) {
        if (sn > lp) {
            //
            // This happens very frequently
            //
            count_too_large_prime_norm++;
            MSGDEBUG("Prime norm %Zd is too large\n", norm);
            return 3;
        }
        mpz_set(p1, norm);
        *m1 = 1;

        return 1;
    }
    unsigned long int retval = 3;

    //
    // Use the TIFA library to completely factor norm
    //
    mpz_array_t*    const factors = alloc_mpz_array(NFACTORS);
    uint32_array_t* const multis  = alloc_uint32_array(NFACTORS);

    MSGDEBUG("factor: tifa_factor() called with norm %Zd\n", norm);
    //
    // _TO_DO_: Testing for every possible problem than can arise is really a
    //          pain. This is only partially done now. Based on my numerous
    //          experiments I think that factoring failures are rare
    //          and not likely to arise in this special case. But I do not
    //          guarantee anything...
    //
    // _NOTE_:  The perfect square test is performed directly in tifa_factor.
    //
    ecode_t ecode = tifa_factor(
                        factors, multis, norm,
                        FIND_COMPLETE_FACTORIZATION
                    );

    switch (ecode) {

    case COMPLETE_FACTORIZATION_FOUND:
        //
        // We have found the complete factorization of the composite norm,
        // including the prime factors' multiplicities. This is the easy case...
        //
        switch (factors->length) {
        case 2: {
            mpz_ptr x1 = factors->data[0];
            mpz_ptr x2 = factors->data[1];

            if ((BITSIZE(x1) > lp) || (BITSIZE(x2) > lp)) {
                MSGDEBUG("norm has a too large prime factor!\n");
                count_bad_factor++;
                goto clear_and_return;
            }
            mpz_set(p1, x1);
            mpz_set(p2, x2);
            *m1 = multis->data[0];
            *m2 = multis->data[1];
            retval = 2;
            break;
        }
        case 1: {
            mpz_ptr x1 = factors->data[0];

            if (BITSIZE(x1) > lp) {
                count_bad_factor++;
                MSGDEBUG("norm has a too large prime factor!\n");
                goto clear_and_return;
            }
            count_perfect_power++;
            mpz_set(p1, x1);
            *m1 = multis->data[0];
            retval = 1;
            break;
        }
        default:
            count_more_two_factors++;
            MSGDEBUG("norm %Zd has more than 2 prime factors\n", norm);
            goto clear_and_return;
        }
        break;

    case PARTIAL_FACTORIZATION_FOUND:
        if (count_partial_fact_found < MAXNMSG) {
            WARNING("TIFA only partially factored norm %Zd\n", norm);
        }
        //
        // Theoretically, this is the difficult case. It should be very rare,
        // particularly in this context, but no there is no way to be 100% sure.
        // So what can happen here?
        //
        // First if norm was a power of prime, it would have been catched when
        // building a coprime base for it and the complete factorization would
        // have been found. So it is safe to assume than we have found _at
        // least_ two coprime factors (but not neccesarily prime!) of norm.
        //
        // If we have found more than 2 coprime factors, this is easy: there is
        // no way norm could then be written as pa^ma * pb^mb with ma, mb prime.
        // So discard the relation.
        //
        // The hard case is when we found 2 coprime factors, with at least one
        // of them not prime. In this case, TIFA was not able to break up the
        // composite factor so it's of no use to try again. The only case we're
        // interested in here in when these unfactored composite are powers of
        // prime. This check has yet to be implemented. For the time being,
        // we just throw out everything even if there's a baby in the bath
        // water...
        //
        // _TO_DO_: Check for prime powers if exactly two coprime factors are
        //          found.
        //
        switch (factors->length) {
        case 2: {
            //
            // _TO_DO_: Check for prime powers!
            //
            count_factorization_pb++;
            fprintf(stderr, "PARTIAL FACTORIZATION!\n");
            goto clear_and_return;
        }
        default:
            count_more_two_factors++;
            MSGDEBUG("norm %Zd has more than 2 prime factors\n", norm);
            goto clear_and_return;
        }
        break;

    case NO_FACTOR_FOUND:
    case FATAL_INTERNAL_ERROR:
    default:
        //
        // Should be rare but I cannot give any warranties here... So what
        // should we do? Try different TIFA algorithms by hand picking them?
        // Or just give up? For the time being, we just discard the relation...
        //
        count_no_factor_found++;
        if (count_no_factor_found < MAXNMSG) {
            WARNING("TIFA found no factor!\n", norm);
        }
        goto clear_and_return;
        break;
    }

  clear_and_return:

    clear_mpz_array(factors);
    clear_uint32_array(multis);
    free(factors);
    free(multis);

    return retval;
}
//-----------------------------------------------------------------------------
inline bool coprime(long a, unsigned long b) {
    //
    // Returns true iff a and b are coprime.
    //
    if (a < 0) {
        a = -a;
    }
    return (gcd_ulint(a, b) == 1);
}
//-----------------------------------------------------------------------------
unsigned long checkrels(char *f, cado_poly cpoly,
                        int verbose, int factorall,
                        size_t mfbr, size_t mfba,
                        unsigned long *primes, unsigned long nprimes) {
    //
    // primes[0] ... primes[nprimes-1] are the small primes omitted in relations
    //
    FILE *fp;

    unsigned long count_large_algprime = 0;
    unsigned long count_large_ratprime = 0;
    unsigned long not_coprime_reports  = 0;
    unsigned long mfba_reports = 0;
    unsigned long mfbr_reports = 0;
    unsigned long nfac;
    unsigned long b;
    long a;

    fbprime_t p;

    char *lineptr = NULL;           // position in buffer for read line
    char line[REL_MAXLENGTH];       // buffer for read line
    char outrel[REL_MAXLENGTH];     // buffer for relation to output

    unsigned long nrels_in  = 0;    // number of relations read
    unsigned long nrels_out = 0;    // number of relations output
    unsigned long outlength;        // length of current relation

    int c;      // number of variables set by sscanf
    int n;      // number of characters consumed by sscanf
    int npr;    // number of primes read by sscanf

    mpz_t norm;
    mpz_t p1;            // first prime factor of norm (if any)
    mpz_t p2;            // second prime factor of norm (if any)

    unsigned int m1 = 1; // multiplicity of first prime factor of norm
    unsigned int m2 = 1; // multiplicity of second prime factor of norm

    unsigned long (*factor_func) (mpz_t, mpz_t, unsigned int* const,
                                  unsigned int* const, const mpz_t,
                                  const size_t);

    mpz_init(norm);
    mpz_init(p1);
    mpz_init(p2);

    if (verbose) {
        MSG("Checking relations from file %s\n", f);
    }
    if (factorall) {
        MSG("Performing full factorization of cofactors' norms...\n");
        factor_func = factor_completely;
    } else {
        MSG("Performing partial factorization of cofactors' norms...\n");
        factor_func = factor;
    }

    fp = fopen(f, "r");
    if (fp == NULL) {
        ERROR("Unable to open file %s\n", f);
        exit(EXIT_FAILURE);
    }
#if DEBUG
    char cont = 0;
#endif

    while (true) {

      next_relation:

#if DEBUG
        MSGDEBUG("[Type enter to read new relation] ");
        cont = getchar();
#endif

        //
        // Read a whole line with a single call to fgets. Not that it makes
        // a huge difference however: most of the in / out time is spent in
        // parsing the line (be it with fscanf or sscanf) and preparing the
        // output (with sprintf).
        //
        // For what it's worth, for c20.rels and on a 867 MHz PowerPC laptop:
        //   sscanf  ~ 50 % of "IO time"
        //   sprintf ~ 40 % of "IO time"
        //
        // _WARNING_: Everything (even input parsing) will break if line is
        //            longer than REL_MAXLENGTH characters.
        //
        if (fgets(line, REL_MAXLENGTH, fp) == NULL) {
            break;
        }
        MSGDEBUG("Processing %s", line);

        if (line[0] == '#') {
            //
            // Skip commented lines...
            //
            printf("%s", line);
            goto next_relation;
        }

        lineptr = line;

        c = sscanf(lineptr, "%ld,%lu:%n", &a, &b, &n);
        if (c == 0 || c == EOF) {
            break;
        }
        lineptr += n;

        nrels_in++;

        outlength  = 0;
        outlength += sprintf(outrel, "%ld,%lu:", a, b);

        //
        // Check that a and b are coprime
        //
        if (!coprime(a, b)) {
            not_coprime_reports++;
            if (not_coprime_reports < MAXNMSG) {
                WARNING("Discarded relation (%ld, %lu) since (a, b) are "
                        "not coprime\n", a, b);
            }
            goto next_relation;
        }
        //
        // Evaluate norm on rational side
        //
        get_rational_norm(norm, cpoly->g, a, b);

        npr = 0;
        //
        // Read primes on rational side
        //
        while (sscanf(lineptr, "%x%n", &p, &n) != 0) {
            npr++;
            //
            // Check that p divides the norm. Otherwise discard the relation.
            //
            if (mpz_fdiv_q_ui(norm, norm, p) != 0) {
                WARNING("Prime %x=%u does not divide norm on rational side "
                      "for a=%ld, b=%lu\n", p, p, a, b);
                goto next_relation;
            }
            //
            // Check that p is smaller than the large prime bound. Otherwise
            // discard relation.
            //
            if (p > (1UL << cpoly->lpbr)) {
                count_large_ratprime++;
                if (count_large_ratprime <= MAXNMSG) {
                    WARNING("Rational prime %x exceeds large prime "
                            "bound 2^%d\n", p, cpoly->lpbr);
                }
                goto next_relation;
            }
            lineptr += n;
            if (lineptr[0] == ':') {
                outlength += sprintf(outrel + outlength, "%x", p);
                lineptr++;
                break;
            } else {
                outlength += sprintf(outrel + outlength, "%x,", p);
                lineptr++;
            }
        }
        //
        // Divide by small primes
        //
        for (unsigned long i = 0; i < nprimes; i++) {
            while (mpz_divisible_ui_p(norm, primes[i])) {
                if (npr == 0) {
                    outlength += sprintf(outrel + outlength, "%lx", primes[i]);
                } else {
                    outlength += sprintf(outrel + outlength, ",%lx", primes[i]);
                }
                mpz_divexact_ui(norm, norm, primes[i]);
                npr++;
            }
        }
        //
        // Check that the residue is smaller than 2^mfbr and factor it.
        // Otherwise discard relation.
        //
        if (BITSIZE(norm) > mfbr) {
            if (++mfbr_reports <= MAXNMSG) {
                WARNING("Rat. residue %Zd exceeds bound for a=%ld, b=%lu\n",
                         norm, a, b);
            }
            goto next_relation;
        }

        nfac = factor_func(p1, p2, &m1, &m2, norm, cpoly->lpbr);

        //
        // Write additional primes from the factorization of norm (if any).
        //
        switch (nfac) {
        case 1:
            MSGDEBUG("checkrels: rational side: nfac = 1\n");
            CHECK_PRIME(p1, cpoly->rlim, a, b);

            if (npr == 0) {
                outlength += gmp_sprintf(outrel + outlength, "%Zx", p1);
            } else {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", p1);
            }
            for (unsigned int i = 1; i < m1; i++) {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", p1);
            }
            outlength += gmp_sprintf(outrel + outlength, ":");
            break;

        case 2:
            MSGDEBUG("checkrels: rational side: nfac = 2\n");
            CHECK_PRIME(p1, cpoly->rlim, a, b);
            CHECK_PRIME(p2, cpoly->rlim, a, b);

            if (npr == 0) {
                outlength += gmp_sprintf(outrel + outlength, "%Zx", p1);
            } else {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", p1);
            }
            for (unsigned int i = 1; i < m1; i++) {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", p1);
            }
            for (unsigned int i = 0; i < m2; i++) {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", p2);
            }
            outlength += gmp_sprintf(outrel + outlength, ":");
            break;

        case 0:
            MSGDEBUG("checkrels: rational side: nfac = 0\n");
            outlength += gmp_sprintf(outrel + outlength, ":");
            break;

        default:
            MSGDEBUG("checkrels: rational side: nfac = 3 -> REL DISCARDED\n");
            goto next_relation;
        }

        //
        // Evaluate norm on algebraic side
        //
        get_algebraic_norm(norm, cpoly->f, cpoly->degree, a, b);

        npr = 0;
        //
        // Read primes on algebraic side
        //
        while (lineptr[0] != '\n' && sscanf(lineptr, "%x%n", &p, &n) != 0) {
            npr++;
            //
            // Check that p divides the norm. Otherwise terminate with error.
            //
            if (mpz_fdiv_q_ui(norm, norm, p) != 0) {
                WARNING("Prime %x=%u does not divide norm on algebraic side "
                      "for a=%ld, b=%lu\n", p, p, a, b);
                goto next_relation;
            }
            //
            // Check that p is smaller than the large prime bound. Otherwise
            // discard relation.
            //
            if (p > (1UL << cpoly->lpba)) {
                count_large_algprime++;
                if (count_large_algprime <= MAXNMSG) {
                    WARNING("Algebraic prime %x exceeds large prime "
                            "bound 2^%d\n", p, cpoly->lpba);
                }
                goto next_relation;
            }
            lineptr += n;
            if (lineptr[0] != ',') {
                outlength += sprintf(outrel + outlength, "%x", p);
                break;
            } else {
                outlength += sprintf(outrel + outlength, "%x,", p);
                lineptr++;
            }
        }
        //
        // Divide by small primes
        //
        for (unsigned long i = 0; i < nprimes; i++) {
            while (mpz_divisible_ui_p(norm, primes[i])) {
                if (npr == 0) {
                    outlength += sprintf(outrel + outlength, "%lx", primes[i]);
                } else {
                    outlength += sprintf(outrel + outlength, ",%lx", primes[i]);
                }
                mpz_divexact_ui(norm, norm, primes[i]);
                npr++;
            }
        }
        //
        // Check that the residue is smaller than 2^mfba and factor it.
        // Otherwise discard relation.
        //
        if (BITSIZE(norm) > (unsigned) mfba) {
            if (++mfba_reports <= MAXNMSG) {
                WARNING("Alg. residue %Zd exceeds bound for a=%ld, b=%lu\n",
                         norm, a, b);
            }
            goto next_relation;
        }

        nfac = factor_func(p1, p2, &m1, &m2, norm, cpoly->lpba);

        //
        // Write additional primes from the factorization of norm (if any).
        //
        switch (nfac) {
        case 1:
            MSGDEBUG("checkrels: algebraic side: nfac = 1\n");
            CHECK_PRIME(p1, cpoly->alim, a, b);

            if (npr == 0) {
                outlength += gmp_sprintf(outrel + outlength, "%Zx", p1);
            } else {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", p1);
            }
            for (unsigned int i = 1; i < m1; i++) {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", p1);
            }
            break;

        case 2:
            MSGDEBUG("checkrels: algebraic side: nfac = 2\n");
            CHECK_PRIME(p1, cpoly->alim, a, b);
            CHECK_PRIME(p2, cpoly->alim, a, b);

            if (npr == 0) {
                outlength += gmp_sprintf(outrel + outlength, "%Zx", p1);
            } else {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", p1);
            }
            for (unsigned int i = 1; i < m1; i++) {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", p1);
            }
            for (unsigned int i = 0; i < m2; i++) {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", p2);
            }
            break;

        case 0:
            MSGDEBUG("checkrels: algebraic side: nfac = 0\n");
            break;

        default:
            MSGDEBUG("checkrels: algebraic side: nfac = 3 -> REL DISCARDED\n");
            goto next_relation;
        }

        nrels_out++;
        if (verbose && !(nrels_out & 4095)) {
            MSG("%lu-th relation found\t\t(%lu read)\n", nrels_out, nrels_in);
        }
        outrel[outlength] = '\0';
        printf("%s\n", outrel);
    }
    fclose(fp);
    mpz_clear(norm);
    mpz_clear(p1);
    mpz_clear(p2);

#if PRINT_SUMMARY

    MSG("-----------------------------------------------------------------\n");
    MSG("File %s\n", f);
    MSG("-----------------------------------------------------------------\n");
    MSG("Miscellaneous information:\n");

    if (!factorall) {
        MSG("    > Found %lu norm(s) to be perfect prime square(s)\n",
            count_perfect_square);
    } else {
        MSG("    > Found %lu norms to be perfect prime power(s)\n",
            count_perfect_power);
        MSG("    > Found %lu partial factorization(s)\n",
            count_partial_fact_found);
    }
    MSG("    > check_prime issued %lu warning(s)\n", check_prime_reports);

    MSG("-----------------------------------------------------------------\n");
    MSG("Discarded relations:\n");
    MSG("    > %lu relations because (a, b) were not coprime\n",
        not_coprime_reports);
    MSG("    > %lu relations because of too large rational residue\n",
        mfbr_reports);
    MSG("    > %lu relations because of too large algebraic residue\n",
        mfba_reports);
    MSG("    > %lu relations because of too large rational prime\n",
        count_large_ratprime);
    MSG("    > %lu relations because of too large algebraic prime\n",
        count_large_algprime);
    MSG("    > %lu relations because of too large prime norm\n",
        count_too_large_prime_norm);
    MSG("    > %lu relations because of too large composite norm\n",
        count_too_large_norm);
    MSG("    > %lu relations because of failed factorization\n",
        count_no_factor_found + count_factorization_pb);
    MSG("    > %lu relations because of norm had more than 2 factors\n",
        count_more_two_factors);
    MSG("    > %lu relations because of non suitable factor of norm\n",
        count_bad_factor);
#endif
    MSG("-----------------------------------------------------------------\n");
    MSG("File %s: output %lu out of %lu relations\n", f, nrels_out, nrels_in);
    MSG("-----------------------------------------------------------------\n");
    
    return nrels_out;
}
//-----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    
    uint64_t t = ckn_mu_secs();
        
    char *polyfilename = NULL;
    cado_poly cpoly;

    unsigned long tot_rels = 0;
    unsigned long p;
    unsigned long tdmax = 1; // bound for trial division
    unsigned long *primes;   // small primes omitted in relations

    int nprimes;
    int mfbr = 0;
    int mfba = 0;
    int nb_files;
    int verbose   = 0;
    int factorall = 0;
    int i;

    MSG("%s.r%s", argv[0], REV);
    for (i = 1; i < argc; i++) {
        MSG(" %s", argv[i]);
    }
    MSG("\n");

    while (argc > 1 && argv[1][0] == '-') {
        if (argc > 1 && strcmp (argv[1], "-v") == 0) {
            verbose++;
            argc--;
            argv++;

        } else if (argc > 1 && strcmp (argv[1], "-factorall") == 0) {
            factorall++;
            argc--;
            argv++;
            tdmax = (tdmax < 20) ? 20 : tdmax;

        } else if (argc > 2 && strcmp (argv[1], "-poly") == 0) {
            polyfilename = argv[2];
            argc -= 2;
            argv += 2;

        } else if (argc > 2 && strcmp (argv[1], "-mfbr") == 0) {
            mfbr = atoi (argv[2]);
            argc -= 2;
            argv += 2;

        } else if (argc > 2 && strcmp (argv[1], "-mfba") == 0) {
            mfba = atoi (argv[2]);
            argc -= 2;
            argv += 2;

        } else if (argc > 2 && strcmp (argv[1], "-t") == 0) {
            tdmax = atoi (argv[2]);
            argc -= 2;
            argv += 2;

        } else  {
            MSG("Usage: %s [-v] [-factorall] [-mfbr <n>] [-mfba <n>] [-t <n>]"
                " -poly <file> rels1 ... relsn\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }
    if (polyfilename == NULL) {
        MSG("Please specify a polynomial file with -poly\n");
        exit(EXIT_FAILURE);
    }
    if (!read_polynomial(cpoly, polyfilename)) {
        MSG("Error reading polynomial file\n");
        exit(EXIT_FAILURE);
    }
    //
    // Check that n divides Res(f,g) (might be useful to factor n)
    //
    check_polynomials(cpoly);
    
    //
    // Default residue bounds are those from file, if not overridden by
    // command line parameters
    //
    if (mfbr == 0) {
        mfbr = cpoly->mfbr;
    }
    if (mfba == 0) {
        mfba = cpoly->mfba;
    }
    //
    // Compute small primes
    //
    primes = (unsigned long*) malloc (tdmax * sizeof(unsigned long));
    nprimes = 0;
    for (p = 2; p < tdmax; p += 1 + (p & 1)) {
        for (i = 0; (i < nprimes) && (primes[i] * primes[i] <= p); i++) {
            if (p % primes[i] == 0) {
                i = nprimes;
                break;
            }
        }
        if (nprimes == 0 || i < nprimes) {
            primes[nprimes++] = p;
        }
    }
    primes = realloc(primes, nprimes * sizeof(unsigned long));

    nb_files = argc - 1;

    while (argc > 1) {
        //
        // Reset reporting info...
        //
        count_more_two_factors     = 0;
        count_too_large_norm       = 0;
        count_bad_factor           = 0;
        count_too_large_prime_norm = 0;
        count_no_factor_found      = 0;
        count_partial_fact_found   = 0;
        count_factorization_pb     = 0;
        count_perfect_square       = 0;
        count_perfect_power        = 0;
        check_prime_reports        = 0;
        
        tot_rels += checkrels(argv[1], cpoly, verbose, factorall, mfbr, mfba,
                            primes, nprimes);
        argc--;
        argv++;
    }
    free (primes);

    if (nb_files > 1) {
        MSG("All input files: output %lu relations\n", tot_rels);
    }
    t = ckn_mu_secs() - t;
    MSG("   Completed in: %.3f seconds\n", t / 1000000.0);
    
    return 0;
}
//-----------------------------------------------------------------------------
