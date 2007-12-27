/**
 * \file    checknorms.c
 * \author  Jerome Milan (heavily based on Paul Zimmermann's checknorms program)
 * \date    Thu Dec 20 2007
 * \brief   Check relations produced by siever and factor residues.
 * \version 0.1
 *
 * This is a slight modification of Paul Zimmermann's checknorms program
 * to use TIFA's factoring primitives instead of ECM. Usage is unchanged:
 *
 * <tt>checknorms -poly c80.poly c80.rels1 [c80.rels2 ...]</tt>
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "cado.h"
#include "utils.h"
#include "tifa.h"

//------------------------------------------------------------------------------
#define VERBOSE 1
//
// _WARNING_: Setting DEBUG to one will prompt for user input before reading
//            a relation and print a truckload of messages.
//
#define DEBUG   0
#if DEBUG
    #define MSGDEBUG(...) gmp_printf(__VA_ARGS__);fflush(stdout);
#else
    #define MSGDEBUG(...) /* voluntarily left empty */
#endif
#ifdef VERBOSE
    #define WARNING(...) gmp_fprintf(stderr, "WARNING: " __VA_ARGS__)
#else
    #define WARNING(...) /* voluntarily left empty */
#endif
#define ERROR(...)   gmp_fprintf(stderr, "ERROR:   " __VA_ARGS__)
#define MSG(...)     gmp_fprintf(stderr, __VA_ARGS__)
//------------------------------------------------------------------------------
#define NFACTORS        8    /* (initial) size of array for norm's factors */
#define MAXNMSG         10   /* max number of warning messages             */
#define REL_MAXLENGTH   1024 /* max length in characters of a relation     */
#define NMILLER_RABIN   12   /* self explanatory                           */
//------------------------------------------------------------------------------
#define IS_PRIME(X)     (0 != mpz_probab_prime_p((X), NMILLER_RABIN))
#define IS_COMPOSITE(X) (0 == mpz_probab_prime_p((X), NMILLER_RABIN))
#define IS_SQUARE(X)    (0 != mpz_perfect_square_p(X))
#define BITSIZE(X)      (mpz_sizeinbase(X, 2))
//------------------------------------------------------------------------------

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
void check_prime(mpz_t p, unsigned long sb, long a, unsigned long b) {
    //
    // Prints a warning message if p <= sb
    //
    static unsigned long reports = 0;
    if (mpz_cmp_ui(p, sb) <= 0) {
        if (++reports <= MAXNMSG) {
            WARNING("(a=%ld, b=%lu): prime %Zd smaller "
                    "than factor base bound %lu\n", a, b, p, sb);
        }
    }
}
//-----------------------------------------------------------------------------
unsigned long factor(mpz_t p1, mpz_t p2, mpz_t norm, size_t lp) {

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
    static unsigned long count_too_large_norm  = 0;
    static unsigned long count_no_factor_found = 0;
    static unsigned long count_perfect_square  = 0;

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
        // than 2^lp or norm has more than 3 factors. Either way, we discard the
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
        count_perfect_square ++;
        if (count_perfect_square < MAXNMSG)
          WARNING("norm is a perfect square!\n");
        mpz_sqrt(p1, norm);
        if (IS_PRIME(p1)) {
            mpz_set(p2, p1);
            return 2;
        }
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
        goto clear_and_return;
    }
    mpz_ptr x1 = norm_factors->data[0];

    if (norm_factors->length == 1) {
        //
        // If the found factor is not prime there is no way norm could be a
        // product of only two primes.
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
                MSGDEBUG("factor p2=%Zd not prime and/or too large!\n", p2);
            }
        } else {
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

        if (   (sx1 <= lp) && (sx2 <= lp)
            && ((sx1 + sx2 == sn) || (sx1 + sx2 == sn + 1)) ) {
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
            MSGDEBUG("factors %Zd and %Zd are prime but sizes mismatch!"
                     " (norm=%Zd)\n", x1, x2, norm);
            MSGDEBUG("size of %Zd is %lu\n", x1, sx1);
            MSGDEBUG("size of %Zd is %lu\n", x2, sx2);
            MSGDEBUG("size of %Zd is %lu\n", norm, sn);
            MSGDEBUG("lp is %lu\n", lp);
        }
    } else {
        /* see comment about FIND_COMPLETE_FACTORIZATION above */
        MSGDEBUG("factor %Zd and / or %Zd is not prime\n", x1, x2);
    }

  clear_and_return:

    clear_mpz_array(norm_factors);
    free(norm_factors);

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
unsigned long checkrels(char *f, cado_poly cpoly, int verbose,
                        size_t mfbr, size_t mfba,
                        unsigned long *primes, unsigned long nprimes) {
    //
    // primes[0]..primes[nprimes-1] are the small primes omitted in relations
    //
    FILE *fp;

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

    mpz_t norm;
    mpz_t p1;
    mpz_t p2;

    mpz_init(norm);
    mpz_init(p1);
    mpz_init(p2);

    if (verbose) {
        MSG("Checking relations from file %s\n", f);
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

        MSGDEBUG("[Type enter to read new relation] ");
#if DEBUG
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
            continue;
        }
        nrels_in++;

        lineptr = line;

        c = sscanf(lineptr, "%ld,%lu:%n", &a, &b, &n);
        if (c == 0 || c == EOF) {
            break;
        }
        lineptr += n;

        outlength  = 0;
        outlength += sprintf(outrel, "%ld,%lu:", a, b);

        //
        // Check that a and b are coprime
        //
        if (!coprime(a, b)) {
            static unsigned long report = 0;
            if (report++ < MAXNMSG) {
                WARNING("Discarded relation (%ld, %lu) since (a, b) are "
                        "not coprime\n", a, b);
            }
            continue;
        }
        //
        // Evaluate norm on rational side
        //
        get_rational_norm(norm, cpoly->g, a, b);

        //
        // Read primes on rational side
        //
        while (sscanf(lineptr, "%x%n", &p, &n) != 0) {
            //
            // Check that p divides the norm. Otherwise terminate with error.
            //
            if (mpz_fdiv_q_ui(norm, norm, p) != 0) {
                ERROR("Prime %x=%u does not divide norm on rational side "
                      "for a=%ld, b=%lu\n", p, p, a, b);
                exit(EXIT_FAILURE);
            }
            //
            // Check that p is smaller than the large prime bound. Otherwise
            // discard relation.
            //
            if (p > (1UL << cpoly->lpbr)) {
                WARNING("Rational prime %x exceeds large prime bound 2^%d\n",
                         p, cpoly->lpbr);
                continue;
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
                outlength += sprintf(outrel + outlength, ",%lx", primes[i]);
                mpz_divexact_ui(norm, norm, primes[i]);
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
            continue;
        }

        nfac = factor(p1, p2, norm, cpoly->lpbr);

        //
        // Write additional primes
        //
        switch (nfac) {
        case 1:
            MSGDEBUG("checkrels: rational side: nfac = 1\n");
            check_prime(p1, cpoly->rlim, a, b);
            outlength += gmp_sprintf(outrel + outlength, ",%Zx:", p1);
            break;

        case 2:
            MSGDEBUG("checkrels: rational side: nfac = 2\n");
            check_prime(p1, cpoly->rlim, a, b);
            check_prime(p2, cpoly->rlim, a, b);
            outlength += gmp_sprintf(outrel + outlength, ",%Zx,%Zx:", p1, p2);
            break;

        case 0:
            MSGDEBUG("checkrels: rational side: nfac = 0\n");
            outlength += gmp_sprintf(outrel + outlength, ":");
            break;

        default:
            MSGDEBUG("checkrels: rational side: nfac = 3 -> REL DISCARDED\n");
            continue;
            break;
        }
        //
        // Evaluate norm on algebraic side
        //
        get_algebraic_norm(norm, cpoly->f, cpoly->degree, a, b);

        //
        // Read primes on algebraic side
        //
        while (lineptr[0] != '\n' && sscanf(lineptr, "%x%n", &p, &n) != 0) {
            //
            // Check that p divides the norm. Otherwise terminate with error.
            //
            if (mpz_fdiv_q_ui(norm, norm, p) != 0) {
                ERROR("Prime %x=%u does not divide norm on algebraic side "
                      "for a=%ld, b=%lu\n", p, p, a, b);
                exit(EXIT_FAILURE);
            }
            //
            // Check that p is smaller than the large prime bound. Otherwise
            // discard relation.
            //
            if (p > (1UL << cpoly->lpba)) {
                WARNING("Algebraic prime %x exceeds large prime bound 2^%d\n",
                        p, cpoly->lpba);
                continue;
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
                outlength += sprintf(outrel + outlength, ",%lx", primes[i]);
                mpz_divexact_ui(norm, norm, primes[i]);
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
            continue;
        }

        nfac = factor(p1, p2, norm, cpoly->lpba);

        //
        // Write additional primes
        //
        switch (nfac) {
        case 1:
            MSGDEBUG("checkrels: algebraic side: nfac = 1\n");
            check_prime(p1, cpoly->alim, a, b);
            outlength += gmp_sprintf(outrel + outlength, ",%Zx", p1);
            break;

        case 2:
            MSGDEBUG("checkrels: algebraic side: nfac = 2\n");
            check_prime(p1, cpoly->alim, a, b);
            check_prime(p2, cpoly->alim, a, b);
            outlength += gmp_sprintf(outrel + outlength, ",%Zx,%Zx", p1, p2);
            break;

        case 0:
            MSGDEBUG("checkrels: algebraic side: nfac = 0\n");
            break;

        default:
            MSGDEBUG("checkrels: algebraic side: nfac = 3 -> REL DISCARDED\n");
            continue;
            break;
        }
        nrels_out++;
        if(verbose && !(nrels_out & 1023)) {
            MSG("%lu-th relation found\t\t(%lu read)\n",
                nrels_out, nrels_in);
        }
        outrel[outlength] = '\0';
        printf("%s\n", outrel);
    }
    fclose(fp);
    mpz_clear(norm);
    mpz_clear(p1);
    mpz_clear(p2);

    MSG("File %s: output %lu out of %lu relations\n", f, nrels_out, nrels_in);

    return nrels_out;
}
//-----------------------------------------------------------------------------
int main(int argc, char *argv[]) {

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
    int verbose = 0;
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
            MSG("Usage: %s [-v] [-mfbr <n>] [-mfba <n>] [-t <n>] -poly <file> "
                "rels1 ... relsn\n", argv[0]);
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
      tot_rels += checkrels(argv[1], cpoly, verbose, mfbr, mfba,
                            primes, nprimes);
      argc--;
      argv++;
    }
    free (primes);

    if (nb_files > 1) {
        MSG("All input files: output %lu relations\n", tot_rels);
    }

    return 0;
}
//-----------------------------------------------------------------------------
