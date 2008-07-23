/**
 * \file    checknorms.c
 * \author  Jerome Milan (heavily based on Paul Zimmermann's checknorms program)
 * \date    Created on Thu Dec 20 2007
 * \date    Last updated on Thu Feb 14 2008
 * \brief   Check relations produced by siever and factor residues.
 * \version 0.3.1
 *
 * This is a slight modification of Paul Zimmermann's checknorms program
 * to use TIFA's factoring primitives instead of ECM. Usage is unchanged:
 *
 * <tt>checknorms -poly c80.poly c80.rels1 [c80.rels2 ...]</tt>
 *
 * Some options have been added. For more information, type
 * <tt>checknorms -h</tt>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>     // for getrusage
#include <sys/resource.h> // for getrusage

#include "cado.h"   // For cado_poly type
#include "utils.h"  // For read_polynomial(...)
#include "tifa.h"   // Use TIFA for residues factorization
#include "funcs.h"  // "Private" header file from the TIFA library to use
                    // gcd_ulint(...)
#include "checknorms.h"

//-----------------------------------------------------------------------------
#if !defined(MAX)
#define MAX(a, b) ( ((a) > (b)) ? (a) : (b) )
#endif
//-----------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    uint64_t t = ckn_mu_secs();

    char *polyfilename = NULL;
    cado_poly cpoly;

    unsigned long nrels_out = 0; // total number of accepted relations

    unsigned long tdmax   = DFLT_TDMAX; // bound for trial division
    unsigned long tdrmax  = 0;    // same bound for rational side (optional)
    unsigned long tdamax  = 0;    // same bound for algebraic side (optional)
    unsigned long *primes = NULL; // small primes omitted in relations
    unsigned long nprimes = 0;    // number of small primes to precompute
    unsigned long npr     = 0;    // number of small primes for rational side
    unsigned long npa     = 0;    // number of small primes for algebraic side

    unsigned int mfbr   = 0;  // bound for rational  residues is 2^mfbr
    unsigned int mfba   = 0;  // bound for algebraic residues is 2^mfba
    unsigned int nfiles = 0;  // number of relation files to check

    unsigned int maxnlp = DFLT_MAX_NLP; // max number of large primes allowed
    unsigned int cmult  = 0;            // count multiplicities in maxlp limit

    unsigned int verbose       = 0; // verbosity switch
    unsigned int keep_comments = 0; // should we keep comments in output file?

    //
    // Parse command line arguments and handle options.
    //
    char* progname = argv[0];
    while (argc > 1 && argv[1][0] == '-') {
        if (argc > 1 && strcmp(argv[1], "-h") == 0) {
            print_help(progname);
            return 0;

        } else if (argc > 1 && strcmp(argv[1], "-v") == 0) {
            verbose++;
            argc--;
            argv++;

        } else if (argc > 1 && strcmp(argv[1], "-com") == 0) {
            keep_comments++;
            argc--;
            argv++;
        
        } else if (argc > 1 && strcmp(argv[1], "-cmult") == 0) {
            cmult++;
            argc--;
            argv++;

        } else if (argc > 1 && strcmp(argv[1], "-maxnlp") == 0) {
            maxnlp = atoi(argv[2]);
            argc -= 2;
            argv += 2;

        } else if (argc > 2 && strcmp(argv[1], "-poly") == 0) {
            polyfilename = argv[2];
            argc -= 2;
            argv += 2;

        } else if (argc > 2 && strcmp(argv[1], "-mfbr") == 0) {
            mfbr = atoi(argv[2]);
            argc -= 2;
            argv += 2;

        } else if (argc > 2 && strcmp(argv[1], "-mfba") == 0) {
            mfba = atoi(argv[2]);
            argc -= 2;
            argv += 2;

        } else if (argc > 2 && strcmp(argv[1], "-t") == 0) {
            tdmax = atoi(argv[2]);
            if (tdmax > MAX_X_PIX) {
                tdmax = MAX_X_PIX;
            }
            argc -= 2;
            argv += 2;

        } else if (argc > 2 && strcmp(argv[1], "-tr") == 0) {
            tdrmax = atoi(argv[2]);
            if (tdrmax > MAX_X_PIX) {
                tdrmax = MAX_X_PIX;
            }
            argc -= 2;
            argv += 2;

        } else if (argc > 2 && strcmp(argv[1], "-ta") == 0) {
            tdamax = atoi(argv[2]);
            if (tdamax > MAX_X_PIX) {
                tdamax = MAX_X_PIX;
            }
            argc -= 2;
            argv += 2;
            
        } else  {
            print_usage(progname);
            return 0;
        }
    }
    if (tdrmax == 0) {
        tdrmax = tdmax;
    }
    if (tdamax == 0) {
        tdamax = tdmax;
    }
    tdmax = MAX(tdamax, tdrmax);
    
    nfiles = argc - 1;
    if (nfiles == 0) {
        MSG("Please specify at least one relation file.\n");
        exit(EXIT_FAILURE);
    }

    //
    // Keep a trace of the command line for debugging...
    //
    MSG("%s.r%s", argv[0], REV);
    for (int i = 1; i < argc; i++) {
        MSG(" %s", argv[i]);
    }
    MSG("\n");

    if (polyfilename == NULL) {
        MSG("Please specify a polynomial file with -poly.\n");
        exit(EXIT_FAILURE);
    }
    cado_poly_init(cpoly);
    if (!cado_poly_read(cpoly, polyfilename)) {
        ERROR("Could not read polynomial file %s.\n", polyfilename);
        exit(EXIT_FAILURE);
    }
    //
    // Check that n divides Res(f,g) (might be useful to factor n)
    //
    check_polynomials(cpoly);

    //
    // Default residue bounds are those from the polynomial file, if not
    // overridden by command line parameters.
    //
    if (mfbr == 0) {
        mfbr = cpoly->mfbr;
    }
    if (mfba == 0) {
        mfba = cpoly->mfba;
    }
    //
    // Compute table of small primes.
    //
    primes  = (unsigned long*) malloc(pi_approx(tdmax) * sizeof(unsigned long));
    nprimes = 0;
    unsigned long int i = 0;
    unsigned long int p = 0;
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
    //
    // Get number of small primes for rational and algebraic side
    //
    unsigned long min = 0;
    unsigned long max = nprimes;
    
    min = 0;
    max = nprimes;
    while ((max - min) > 1) {
        i = (max + min) / 2;
        if (tdamax <= primes[i]) {
            max = i;
        } else {
            min = i;
        }
    }
    npa = i;
    
    min = 0;
    max = nprimes;
    while ((max - min) > 1) {
        i = (max + min) / 2;
        if (tdrmax <= primes[i]) {
            max = i;
        } else {
            min = i;
        }
    }
    npr = i;    
    //
    // Check and complete relations from each relation file...
    //
    while (argc > 1) {
        //
        // Reset reporting info...
        //
        count_too_many_factors     = 0;
        count_too_large_prime_norm = 0;
        count_too_large_factor     = 0;
        count_no_factor_found      = 0;
        count_partial_fact_found   = 0;
        count_perfect_power        = 0;
        count_check_prime          = 0;
        count_large_algprime       = 0;
        count_large_ratprime       = 0;
        count_not_coprime          = 0;
        count_not_dividing         = 0;
        count_larger_mfba          = 0;
        count_larger_mfbr          = 0;
        
        nrels_out += checkrels(
                       argv[1], cpoly, verbose, keep_comments, maxnlp, cmult,
                       mfbr, mfba, primes, npr, npa
                     );
        argc--;
        argv++;
    }
    free(primes);
    
    if (nfiles > 1) {
        MSG("\n");
        MSG("--------------------------------------------------------------\n");
        MSG(" All input files: output %lu out of %lu relations\n",
            nrels_out, total_nrelations_read);
        MSG("--------------------------------------------------------------\n");
    }
    t = ckn_mu_secs() - t;
    MSG("\nCompleted in: %.3f seconds\n", t / 1000000.0);
    
    return 0;
}
//-----------------------------------------------------------------------------
unsigned long checkrels(char *f, cado_poly cpoly,
                        unsigned int verbose, unsigned int keep_comments,
                        unsigned int maxnlp, unsigned int cmult,
                        size_t mfbr, size_t mfba,
                        unsigned long *primes,
                        unsigned long nprimes_rat, unsigned long nprimes_alg) {
    //
    // primes[0], primes[1], primes[2] ... should be the small primes
    // omitted in relations.
    //
    // _WARNING_: This functions uses and/or modifies the count_* global
    //            variables declared in checknorms.h.
    //
    FILE *fp;

    unsigned long ecode;
    unsigned long b;
    long a;

    char *lineptr = NULL;           // position in buffer for read line
    char line[REL_MAXLENGTH];       // buffer for read line
    char outrel[REL_MAXLENGTH];     // buffer for relation to output

    unsigned long nrels_in  = 0;    // number of relations read
    unsigned long nrels_out = 0;    // number of relations output
    unsigned long outlength;        // length of current relation

    int c;        // number of variables set by sscanf
    int n;        // number of characters consumed by sscanf
    int npr;      // number of primes read by sscanf
    fbprime_t p;  // current prime read from relation file

    mpz_t norm;
    mpz_init(norm);

    // Factors of norm a priori not in the factor base and their
    // associated multiplicities
    mpz_array_t*    const factors = alloc_mpz_array(NFACTORS);
    uint32_array_t* const multis  = alloc_uint32_array(NFACTORS);

    if (verbose) {
        MSG("Checking relations from file %s...\n", f);
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
        //
        // Reset factors and multiplicities before getting a new relation.
        //
        reset_mpz_array(factors);
        multis->length = 0;

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
            if (keep_comments) {
                printf("%s", line);
            }
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
        // Check that a and b are coprime.
        //
        if (!coprime(a, b)) {
            count_not_coprime++;
            if (count_not_coprime < MAXNMSG) {
                WARNING("Discarded relation (%ld, %lu) since (a, b) are "
                        "not coprime\n", a, b);
            }
            goto next_relation;
        }
        //
        // Evaluate norm on rational side.
        //
        get_rational_norm(norm, cpoly->g, a, b);

        npr = 0;
        //
        // Read primes on rational side.
        //
        while (sscanf(lineptr, "%x%n", &p, &n) != 0) {          
            npr++;
            //
            // Check that p divides the norm. Otherwise discard the relation.
            //
            if (mpz_fdiv_q_ui(norm, norm, p) != 0) {
                count_not_dividing++;
                if (count_not_dividing <= MAXNMSG) {
                    WARNING("Prime %x=%u does not divide norm on rational side "
                            "for a=%ld, b=%lu\n", p, p, a, b);
                }
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
                //
                // Last factor on the rational side
                //
                outlength += sprintf(outrel + outlength, "%x", p);
                lineptr++;
                break;
            } else {
                outlength += sprintf(outrel + outlength, "%x,", p);
                lineptr++;
            }
        }
        if (npr == 0) {
            //
            // Handle special case: no factor on the rational side
            //
            lineptr++;
        }
        
        //
        // Divide by small primes ommited during sieving.
        //
        for (unsigned long i = 0; i < nprimes_rat; i++) {
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
            count_larger_mfbr++;
            if (count_larger_mfbr <= MAXNMSG) {
                WARNING("Rat. residue %Zd exceeds bound for a=%ld, b=%lu\n",
                         norm, a, b);
            }
            goto next_relation;
        }

        ecode = factor_completely(
                    factors, multis, norm, cpoly->lpbr, maxnlp, cmult
                );

        //
        // Write additional primes from the factorization of norm (if any).
        //
        switch (ecode) {
        case NORM_IS_ACCEPTED: {
            MSGDEBUG("checkrels: rational side: NORM_IS_ACCEPTED\n");

            mpz_ptr  fac = factors->data[0];
            uint32_t mul = multis->data[0];

            CHECK_PRIME(fac, cpoly->rlim, a, b);

            if (npr != 0) {
                outlength += sprintf(outrel + outlength, ",");
            }
            outlength += gmp_sprintf(outrel + outlength, "%Zx", fac);

            for (uint32_t j = 1; j < mul; j++) {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", fac);
            }
            for (uint32_t i = 1; i < factors->length; i++) {
                fac = factors->data[i];
                mul = multis->data[i];

                CHECK_PRIME(fac, cpoly->rlim, a, b);

                for (uint32_t j = 0; j < mul; j++) {
                    outlength += gmp_sprintf(outrel + outlength, ",%Zx", fac);
                }
            }
            break;
        }
        case NORM_IS_ONE:
            MSGDEBUG("checkrels: rational side: NORM_IS_ONE\n");
            break;

        case NORM_IS_REJECTED:
            MSGDEBUG("checkrels: rational side: NORM_IS_REJECTED\n");
            goto next_relation;

        case ERROR_FACTORING_NORM:
        default:
            MSGDEBUG("checkrels: rational side: ERROR_FACTORING_NORM\n");
            goto next_relation;
        }
        outlength += gmp_sprintf(outrel + outlength, ":");

        //
        // Evaluate norm on algebraic side.
        //
        get_algebraic_norm(norm, cpoly->f, cpoly->degree, a, b);

        npr = 0;
        //
        // Read primes on algebraic side.
        //
        while (lineptr[0] != '\n' && sscanf(lineptr, "%x%n", &p, &n) != 0) {
            npr++;
            //
            // Check that p divides the norm. Otherwise terminate with error.
            //
            if (mpz_fdiv_q_ui(norm, norm, p) != 0) {
                count_not_dividing++;
                if (count_not_dividing <= MAXNMSG) {
                    WARNING("Prime %x=%u does not divide norm on algebraic "
                            "side for a=%ld, b=%lu\n", p, p, a, b);
                }
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
                //
                // Last factor on the rational side
                //
                outlength += sprintf(outrel + outlength, "%x", p);
                break;
            } else {
                outlength += sprintf(outrel + outlength, "%x,", p);
                lineptr++;
            }
        }
        //
        // Divide by small primes ommited during sieving.
        //
        for (unsigned long i = 0; i < nprimes_alg; i++) {
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
            count_larger_mfba++;
            if (count_larger_mfba <= MAXNMSG) {
                WARNING("Alg. residue %Zd exceeds bound for a=%ld, b=%lu\n",
                         norm, a, b);
            }
            goto next_relation;
        }
        //
        // Reset factors and multiplicities before new factorization.
        //
        reset_mpz_array(factors);
        multis->length = 0;
        
        ecode = factor_completely(
                    factors, multis, norm, cpoly->lpba, maxnlp, cmult
                );

        //
        // Write additional primes from the factorization of norm (if any).
        //
        switch (ecode) {
        case NORM_IS_ACCEPTED: {
            MSGDEBUG("checkrels: algebraic side: NORM_IS_ACCEPTED\n");

            mpz_ptr  fac = factors->data[0];
            uint32_t mul = multis->data[0];

            CHECK_PRIME(fac, cpoly->alim, a, b);

            if (npr != 0) {
                outlength += sprintf(outrel + outlength, ",");
            }
            outlength += gmp_sprintf(outrel + outlength, "%Zx", fac);

            for (uint32_t j = 1; j < mul; j++) {
                outlength += gmp_sprintf(outrel + outlength, ",%Zx", fac);
            }

            for (uint32_t i = 1; i < factors->length; i++) {
                fac = factors->data[i];
                mul = multis->data[i];

                CHECK_PRIME(fac, cpoly->alim, a, b);

                for (uint32_t j = 0; j < mul; j++) {
                    outlength += gmp_sprintf(outrel + outlength, ",%Zx", fac);
                }
            }
            break;
        }
        case NORM_IS_ONE:
            MSGDEBUG("checkrels: algebraic side: NORM_IS_ONE\n");
            break;

        case NORM_IS_REJECTED:
            MSGDEBUG("checkrels: algebraic side: NORM_IS_REJECTED\n");
            goto next_relation;

        case ERROR_FACTORING_NORM:
        default:
            MSGDEBUG("checkrels: algebraic side: ERROR_FACTORING_NORM\n");
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
    clear_mpz_array(factors);
    clear_uint32_array(multis);

#if PRINT_SUMMARY
    MSG("-----------------------------------------------------------------\n");
    MSG("File %s\n", f);
    MSG("-----------------------------------------------------------------\n");
    MSG("Miscellaneous information:\n");
    MSG("    > Found %lu norms to be perfect prime power(s)\n",
        count_perfect_power);
    MSG("    > Found %lu partial factorization(s)\n",
        count_partial_fact_found);
    MSG("    > check_prime issued %lu warning(s)\n", count_check_prime);
    MSG("-----------------------------------------------------------------\n");
    MSG("Discarded relations:\n");
    if (count_not_coprime != 0) {
        MSG("    > %lu relations because (a, b) were not coprime\n",
            count_not_coprime);
    }
    if (count_not_dividing != 0) {
        MSG("    > %lu relations because prime not dividing norm\n",
            count_not_dividing);
    }
    if (count_larger_mfbr != 0) {
        MSG("    > %lu relations because of too large rational residue\n",
            count_larger_mfbr);
    }
    if (count_larger_mfba != 0) {
        MSG("    > %lu relations because of too large algebraic residue\n",
            count_larger_mfba);
    }
    if (count_large_ratprime != 0) {
        MSG("    > %lu relations because of too large rational prime\n",
            count_large_ratprime);
    }
    if (count_large_algprime != 0) {
        MSG("    > %lu relations because of too large algebraic prime\n",
            count_large_algprime);
    }
    if (count_too_large_prime_norm != 0) {
        MSG("    > %lu relations because of too large prime norm\n",
            count_too_large_prime_norm);
    }
    if (count_no_factor_found != 0) {
        MSG("    > %lu relations because of failed factorization\n",
            count_no_factor_found);
    }
    if (count_too_many_factors != 0) {
        MSG("    > %lu relations because norm had too many prime factors\n",
            count_too_many_factors);
    }
    if (count_too_large_factor != 0) {
        MSG("    > %lu relations because norm had a too large prime factor\n",
            count_too_large_factor);
    }
#endif

    MSG("-----------------------------------------------------------------\n");
    MSG("File %s: output %lu out of %lu relations\n", f, nrels_out, nrels_in);
    MSG("-----------------------------------------------------------------\n");
    
    total_nrelations_read += nrels_in;
    
    return nrels_out;
}
//-----------------------------------------------------------------------------
unsigned long factor_completely(mpz_array_t*    const factors,
                                uint32_array_t* const multis,
                                const mpz_t norm,
                                const size_t lp,
                                const unsigned int maxnlp,
                                const unsigned int cmult) {
    MSGDEBUG("factor_completely: entered\n");

    //
    // _WARNING_: Uses and modifies the following global variables:
    //              - count_too_large_prime_norm
    //              - count_too_many_factors
    //              - count_no_factor_found
    //              - count_perfect_power
    //
    if (mpz_cmp_ui(norm, 1) == 0) {
        MSGDEBUG("norm is one!\n", norm);
        return NORM_IS_ONE;
    }
    unsigned long sn = BITSIZE(norm);

    if (IS_PRIME(norm)) {
        if (sn > lp) {
            //
            // This happens very frequently
            //
            count_too_large_prime_norm++;
            MSGDEBUG("Prime norm %Zd is too large\n", norm);
            return NORM_IS_REJECTED;
        }
        append_mpz_to_array(factors, norm);
        append_uint32_to_array(multis, 1);

        return NORM_IS_ACCEPTED;
    }
    unsigned long int retval = ERROR_FACTORING_NORM;

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

    case COMPLETE_FACTORIZATION_FOUND: {
        //
        // We have found the complete factorization of the composite norm,
        // including the prime factors' multiplicities. This is the easy case...
        //
        if (cmult != 0) {
            uint32_t mtot = 0;
            for (uint32_t i = 0; i < multis->length; i++) {
                mtot += multis->data[i];
            }
            if (mtot > maxnlp) {
                count_too_many_factors++;
                if (count_too_many_factors <= MAXNMSG) {
                    WARNING("Residue %Zd has too many prime factors\n", norm);
                }
                retval = NORM_IS_REJECTED;
                goto clean_and_return;
            }
        }
        if (factors->length == 1) {
            count_perfect_power++;
        } else if (factors->length > maxnlp) {
            count_too_many_factors++;
            if (count_too_many_factors <= MAXNMSG) {
                WARNING("Residue %Zd has too many prime factors\n", norm);
            }
            retval = NORM_IS_REJECTED;
            goto clean_and_return;
        }
        for (uint32_t i = 0; i < factors->length; i++) {
            if (BITSIZE(factors->data[i]) > lp) {
                count_too_large_factor++;
                if (count_too_large_factor <= MAXNMSG) {
                    WARNING(
                        "Residue %Zd has too large prime factor %Zd "
                        "(%lu bits)\n", norm, factors->data[i],
                        BITSIZE(factors->data[i])
                    );
                }
                retval = NORM_IS_REJECTED;
                goto clean_and_return;
            }
        }
        retval = NORM_IS_ACCEPTED;
        break;
    }
    case PARTIAL_FACTORIZATION_FOUND:
        //
        // _NOTE_: This is the more tedious case. Try to do something
        //         smarter than just returning with an error...
        //
        count_partial_fact_found++;
        if (count_partial_fact_found <= MAXNMSG) {
            WARNING("TIFA only partially factored norm %Zd\n", norm);
        }
        retval = ERROR_FACTORING_NORM;
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
        retval = ERROR_FACTORING_NORM;
        break;
    }

  clean_and_return:

    return retval;
}
//-----------------------------------------------------------------------------
void get_rational_norm(mpz_t norm, mpz_t *g, long a, unsigned long b) {
    mpz_mul_si(norm, g[1], a);
    mpz_addmul_ui(norm, g[0], b);
    mpz_abs(norm, norm);
}
//-----------------------------------------------------------------------------
void get_algebraic_norm(mpz_t norm, mpz_t *f, int d, long a, unsigned long b) {
    
    mpz_t bpow;
    mpz_init_set_ui(bpow, 1);
    mpz_set(norm, f[d]);
    
    d--;
    while (d >= 0) {
        mpz_mul_si(norm, norm, a);
        mpz_mul_ui(bpow, bpow, b);
        mpz_addmul(norm, f[d], bpow);
        d--;
    }
    mpz_clear(bpow);
    mpz_abs(norm, norm);
}
//-----------------------------------------------------------------------------
static inline void
check_prime(mpz_t p, unsigned long sb, long a, unsigned long b) {
    //
    // _WARNING_: Uses and modifies the following global variable:
    //            - count_check_prime
    //
    if (mpz_cmp_ui(p, sb) <= 0) {
        count_check_prime++;
        if (count_check_prime <= MAXNMSG) {
            WARNING("(a=%ld, b=%lu): prime %Zd smaller "
                    "than factor base bound %lu\n", a, b, p, sb);
        }
    }
}
//-----------------------------------------------------------------------------
inline bool coprime(long a, unsigned long b) {
    if (a < 0) {
        a = -a;
    }
    return (gcd_ulint(a, b) == 1);
}
//-----------------------------------------------------------------------------
uint32_t pi_approx(uint64_t x) {
    size_t i = 0;
    while (x  > pi_table_x[i]) {
        i++;
    }
    return pi_table_pi[i];
}
//-----------------------------------------------------------------------------
void print_usage(char* progname) {
    printf("Usage:\n");
    printf("------\n");
    printf("%14s [-h] [-v] [-cmult]\n", progname);
    printf("%14s [-maxnlp <num>] [-mfbr <num>] [-mfba <num>]\n", "");
    printf("%14s [-t <num>] [-tr <num>] [-ta <num>]\n", "");
    printf("%14s -poly <file> <relfile_1> [<relfile_2> ... <relfile_n>]\n\n",
           "");
    printf("Type %s -h for more information.\n", progname);
}
//-----------------------------------------------------------------------------
void print_help(char* progname) {
    printf("Usage:\n");
    printf("------\n");
    printf("%14s [-h] [-v] [-com] [-cmult]\n", progname);
    printf("%14s [-maxnlp <num>] [-mfbr <num>] [-mfba <num>]\n", "");
    printf("%14s [-t <num>] [-tr <num>] [-ta <num>]\n", "");
    printf("%14s -poly <file> <relfile_1> [<relfile_2> ... <relfile_n>]\n\n",
           "");
    printf("Mandatory arguments:\n");
    printf("--------------------\n");
    printf("  -poly FILE\n");
    printf("      CADO polynomial file.\n");
    printf("  <relfile_i> FILES\n");
    printf("      Space-separated list of relation files obtained from CADO "
           "siever.\n\n");
    printf("Options:\n");
    printf("--------\n");
    printf("  -v\n");
    printf("      Turn verbose mode on.\n");
    printf("  -com\n");
    printf("      Copy comments from input files to output files.\n");
    printf("      By default, comments are not copied to output files.\n");
    printf("  -h\n");
    printf("      Print this help.\n");
    printf("  -cmult\n");
    printf("      Take into account multiplicities of residues' factors.\n");
    printf("  -maxnlp NUM\n");
    printf("      Maximum number of large primes to allow.\n");
    printf("      Multiplicities are taken into account if the option -cmult ");
    printf("is given.\n");
    printf("      Default: %i\n", DFLT_MAX_NLP);
    printf("  -mfbr NUM\n");
    printf("      Bound for rational residues (in bits).\n");
    printf("      Default: use value from polynomial file\n");
    printf("  -mfba NUM\n");
    printf("      Bound for algebraic residues (in bits).\n");
    printf("      Default: use value from polynomial file\n");
    printf("  -t NUM\n");
    printf("      Bound for largest prime for trial division.\n");
    printf("      Default: %i\n", DFLT_TDMAX);
    printf("  -tr NUM\n");
    printf("      Bound for largest prime for trial division on rational "
           "side.\n");
    printf("      Overrides generic value given by -t option.\n");
    printf("      Default: %i or value given by -t option\n", DFLT_TDMAX);
    printf("  -ta NUM\n");
    printf("      Bound for largest prime for trial division on algebraic "
           "side.\n");
    printf("      Overrides generic value given by -t option.\n");
    printf("      Default: %i or value given by -t option\n", DFLT_TDMAX);
}
//-----------------------------------------------------------------------------
uint64_t ckn_mu_secs() {
    struct rusage res[1];
    getrusage(RUSAGE_SELF, res);
    
    uint64_t r;
    r  = (uint64_t) res->ru_utime.tv_sec + res->ru_stime.tv_sec;
    r *= (uint64_t) 1000000UL;
    r += (uint64_t) res->ru_utime.tv_usec + res->ru_stime.tv_usec;
    return r;
}
//-----------------------------------------------------------------------------
