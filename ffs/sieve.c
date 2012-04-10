#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "macros.h"

#include "types.h"
#include "fppol.h"
#include "ffspol.h"
#include "qlat.h"
#include "fb.h"
#include "latsieve.h"
#include "norm.h"
#include "timing.h"
#include "ijvec.h"
#include "params.h"
#include "sublat.h"
#include "smoothness.h"
#include "polyfactor.h"
#include "buckets.h"

int my_factor_survivor(fppol_t a, fppol_t b, ffspol_t* F, int *B, sq_t q,
        int side) 
{
    fppol_t Nab;
    fppol_init(Nab);
    for (int twice = 0; twice < 2; twice++) {
        ffspol_norm(Nab, F[twice], a, b);
        if (side == twice) {
            fppol_t qq;
            fppol_init(qq);
            fppol_set_sq(qq, q);
            fppol_div(Nab, Nab, qq);
            fppol_clear(qq);
        }
        if (!fppol_is_smooth(Nab, B[twice])) {
            fppol_clear(Nab);
            return 0;
        }
    }

    // OK. Now we know that we have a relation.
    // Let's factor it completely and print it.
    fppol_fact_t factors;
    fppol_fact_init(factors);
    fppol_out(stdout, a); printf(",");
    fppol_out(stdout, b); printf(":");
    for (int twice = 0; twice < 2; twice++) {
        ffspol_norm(Nab, F[twice], a, b);
        fppol_factor(factors, Nab);
        fppol_fact_out(stdout, factors);
        if (!twice)
            printf(":");
    }
    printf("\n");
    fppol_fact_clear(factors);
    fppol_clear(Nab);
    return 1;
}


void usage(const char *argv0, const char * missing)
{
    fprintf(stderr, "Usage: %s [options | optionfile] \n", argv0);
    fprintf(stderr, "  Command line options have the form '-optname optval'\n");
    fprintf(stderr, "  and take the form 'optname=optval' in the optionfile\n");
    fprintf(stderr, "List of options (a * means mandatory, default value in []):\n");
    fprintf(stderr, "  pol0 *           function field polynomial on side 0\n");
    fprintf(stderr, "  pol1 *           function field polynomial on side 1\n");
    fprintf(stderr, "  fb0  *           factor base file on side 0\n");
    fprintf(stderr, "  fb1  *           factor base file on side 1\n");
    fprintf(stderr, "  fbb0             factor base bound on side 0\n");
    fprintf(stderr, "  fbb1             factor base bound on side 1\n");
    fprintf(stderr, "  I    *           degree bound for i\n");
    fprintf(stderr, "  J    *           degree bound for j\n");
    fprintf(stderr, "  lpb0 *           large prime bound on side 0\n");
    fprintf(stderr, "  lpb1 *           large prime bound on side 1\n");
    fprintf(stderr, "  thresh0 [2*lpb0] survivor threshold on side 0\n");
    fprintf(stderr, "  thresh1 [2*lpb1] survivor threshold on side 1\n");
    fprintf(stderr, "  q    *           q-poly of the special-q\n");
    fprintf(stderr, "  rho  *           rho-poly of the special-q\n");
    fprintf(stderr, "  sqfile *         file of special-q's. This option, disable the two previous\n");
    fprintf(stderr, "  sqside           side (0 or 1) of the special-q\n");
    fprintf(stderr, "  firstsieve       side (0 or 1) to sieve first\n");
    fprintf(stderr, "  sublat           toggle the sublattice sieving\n");
 
    if (missing != NULL)
        fprintf(stderr, "Missing parameter: %s\n", missing);
    exit(1);
}

// usage: ./a.out q rho
int main(int argc, char **argv)
{
    ffspol_t ffspol[2]; 
    qlat_t qlat;
    factor_base_t FB[2];
    int fbb[2] = {0, 0};
    int I=0, J=0;  
    int lpb[2] = {0, 0};  
    unsigned int threshold[2] = {0, 0};  
    char *argv0 = argv[0];
    int want_sublat = 0;
    const char *sqfile = NULL;
    FILE *sqFile = NULL;
    int sqside = 0;
    int firstsieve = 0;

    param_list pl;
    param_list_init(pl);
    param_list_configure_knob(pl, "-sublat", &want_sublat);
    argv++, argc--;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        /* Could also be a parameter file */
        FILE *f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f);
            fclose(f);
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage(argv0, NULL);
    }
    param_list_print_command_line(stdout, pl);

    // read function field polynomials
    {
        const char * polstr;
        ffspol_init(ffspol[0]);
        ffspol_init(ffspol[1]);
        polstr = param_list_lookup_string(pl, "pol0");
        if (polstr == NULL) usage(argv0, "pol0");
        ffspol_set_str(ffspol[0], polstr);
        polstr = param_list_lookup_string(pl, "pol1");
        if (polstr == NULL) usage(argv0, "pol1");
        ffspol_set_str(ffspol[1], polstr);
    }
    // read various bounds
    param_list_parse_int(pl, "I", &I); 
    param_list_parse_int(pl, "J", &J); 
    param_list_parse_int(pl, "lpb0", &lpb[0]);
    param_list_parse_int(pl, "lpb1", &lpb[1]);
    param_list_parse_int(pl, "fbb0", &fbb[0]);
    param_list_parse_int(pl, "fbb1", &fbb[1]);
    param_list_parse_uint(pl, "thresh0", &threshold[0]);
    param_list_parse_uint(pl, "thresh1", &threshold[1]);
    if (I == 0) usage(argv0, "I");
    if (J == 0) usage(argv0, "J");
    if (lpb[0] == 0) usage(argv0, "lpb0");
    if (lpb[1] == 0) usage(argv0, "lpb1");
    if (threshold[0] == 0) 
        threshold[0] = 2*lpb[0];
    if (threshold[1] == 0) 
        threshold[1] = 2*lpb[1];
    // check if there is an sqfile
    {
        sqfile = param_list_lookup_string(pl, "sqfile");
        if (sqfile != NULL) {
            sqFile = fopen(sqfile, "r");
            if (sqFile == NULL) {
                fprintf(stderr, "Could not open sqfile %s for reading\n",
                        sqfile);
                exit (EXIT_FAILURE);
            }
        }
    }
    // read q, rho
    if (!sqfile) {
        const char *sqstr;
        int noerr;
        sqstr = param_list_lookup_string(pl, "q");
        if (sqstr == NULL) usage(argv0, "q");
        noerr = sq_set_str(qlat->q, sqstr);
        if (!noerr) {
            fprintf(stderr, "Could not parse q: %s\n", sqstr);
            exit(EXIT_FAILURE);
        }
        sqstr = param_list_lookup_string(pl, "rho");
        if (sqstr == NULL) usage(argv0, "rho");
        noerr = sq_set_str(qlat->rho, sqstr);
        if (!noerr) {
            fprintf(stderr, "Could not parse rho: %s\n", sqstr);
            exit(EXIT_FAILURE);
        }
        if (sq_deg(qlat->q) > lpb[0]) {
            fprintf(stderr, "WARNING: not a good idea to have a special-q beyond the large prime bound!\n");
        }
    }
    param_list_parse_int(pl, "sqside", &sqside);
    if (sqside != 0 && sqside != 1) {
        fprintf(stderr, "sqside must be 0 or 1\n");
        exit(EXIT_FAILURE);
    }
    param_list_parse_int(pl, "firstsieve", &firstsieve);
    if (firstsieve != 0 && firstsieve != 1) {
        fprintf(stderr, "firstsieve must be 0 or 1\n");
        exit(EXIT_FAILURE);
    }
   
    // Read the factor bases
    {
        const char *filename;
        int noerr;
        for (int i = 0; i < 2; ++i) {
            char param[4] = {'f', 'b', '0', '\0'};
            if (i == 1) 
                param[2] = '1';
            filename = param_list_lookup_string(pl, param);
            if (filename == NULL) usage(argv0, param);
            double tm = seconds();
            noerr = factor_base_init(FB[i], filename, fbb[i]);
            fprintf(stdout, "# Reading factor base %d took %1.1f s\n", 
                    i, seconds()-tm);
            if (!noerr) {
                fprintf(stderr, "Could not read %s: %s\n", param, filename);
                exit(EXIT_FAILURE);
            }
        }
    }

    param_list_clear(pl);

    if (want_sublat) {
        fprintf(stderr, "# WARNING: sublattices seem to be broken, and won't be fixed soon.\n");
    }

#ifdef USE_F2
    sublat_ptr sublat;
    if (want_sublat)
        sublat = &nine_sublat[0];
    else
        sublat = &no_sublat[0];
#else
    if (want_sublat)
        fprintf(stderr, "# Sorry, no sublattices in characteristic > 2. Ignoring the 'sublat' option\n");
    sublat_ptr sublat = &no_sublat[0];
#endif


    // Most of what we do is at the sublattice level. 
    // So we fix I and J accordingly.
    I -= sublat->deg;
    J -= sublat->deg;

    // Allocated size of the sieve array.
    unsigned IIJJ = ijvec_get_max_pos(I, J);

    qlat->side = sqside; 

    double tot_time = seconds();
    int tot_nrels = 0;
    int tot_sq = 0;
    
    //
    // Begin of loop over special-q's
    //
    int end_list_of_sq = 0;
    do {
        // Select next special-q
        if (sqfile) {
            int noerr;
            noerr = sq_inp(qlat->q, sqFile);
            if (!noerr) break;
            noerr = sq_inp(qlat->rho, sqFile);
            if (!noerr) break;
            printf("#######################################\n");
            printf("# Have read next special q from sqfile.\n");
        } else // in the case where there is no file, we take the one that is
               // given on command line.
            end_list_of_sq = 1;
        tot_sq++;

        double t_tot = seconds();

        // Check the given special-q
        if (!is_valid_sq(qlat, ffspol[sqside])) {
            fprintf(stderr, "Error: the rho = ");
            sq_out(stderr, qlat->rho);
            fprintf(stderr, " is not a root modulo ");
            sq_out(stderr, qlat->q);
            fprintf(stderr, " of the polynomial %d\n", sqside);
            exit(EXIT_FAILURE);
        }

        // Reduce the q-lattice
        int noerr = skewGauss(qlat, 0);
        ASSERT_ALWAYS(noerr);
        print_qlat_info(qlat);

        // Precompute lambda for each element of the factor bases.
        for (int i = 0; i < 2; ++i) {
            double tm = seconds();
            factor_base_precomp_lambda(FB[i], qlat, sublat);
            fprintf(stdout, "# Precomputing lambda on side %d took %1.1f s\n",
                    i, seconds()-tm);
        }

        double t_norms = 0;
        double t_sieve = 0;
        double t_buckets = 0;
        double t_cofact = 0;
        double t_initS = 0;
        int nrels = 0;

        // Loop on all sublattices
        // In the no_sublat case, this loops degenerates into one pass, since
        // nb = 1.
        for (sublat->n = 0; sublat->n < sublat->nb; sublat->n++) {
            if (use_sublat(sublat)) {
                fprintf(stdout, "# Sublattice (");
                fppol16_out(stdout, sublat->lat[sublat->n][0]);
                fprintf(stdout, ", ");
                fppol16_out(stdout, sublat->lat[sublat->n][1]);
                fprintf(stdout, ") :\n");
            }

            t_initS -= seconds();
            // Allocate and init the sieve space
            uint8_t *S;
            S = (uint8_t *) malloc(IIJJ*sizeof(uint8_t));
            ASSERT_ALWAYS(S != NULL);
            memset(S, 0, IIJJ*sizeof(uint8_t));

            // Kill trivial positions.
            // When there are no sublattices:
            //   (i,0) for i != 1
            //   (0,j) for j != 1
            // When using sublattices, just the position (0,0)
            {
                S[0] = 255;  // that's (0,0)
                if (!use_sublat(sublat)) {
                    ij_t i, j;
                    for (ij_set_zero(j), ij_monic_set_next(j, j, J);
                            ij_monic_set_next(j, j, J); )
                        S[ijvec_get_start_pos(j, I, J)] = 255;
                    for (ij_set_zero(i), ij_set_next(i, i, I);
                            ij_set_next(i, i, I); )
                        S[ijvec_get_offset(i, I)] = 255;
                }
            }
// If the norm computation is too expensive, it might pay to activate
// this block.
#if 0
            // Kill positions where gcd(hat i, hat j) != 1
            {
                ij_t i, j;
                ij_t hati, hatj, g;
                int rci, rcj = 1;
                for (ij_set_zero(j); rcj; rcj = ij_monic_set_next(j, j, J)) {
                    ijpos_t start = ijvec_get_start_pos(j, I, J);
                    rci = 1;
                    for (ij_set_zero(i); rci; rci = ij_set_next(i, i, I)) {
                        ijpos_t pos = start + ijvec_get_offset(i, I);
                        ij_convert_sublat(hati, hatj, i, j, sublat);
                        ij_gcd(g, hati, hatj);
                        if (ij_deg(g) != 0 && ij_deg(hati)>0  && ij_deg(hatj)>0)
                            S[pos] = 255;
                    }
                }
            }
#endif
            t_initS += seconds();

            buckets_t buckets;
            // FIXME: The bucket capacity is hardcoded for the moment.
            buckets_init(buckets, I, J, 1<<20);
            if (sublat->n == 0)
                print_bucket_info(buckets);
            for (int twice = 0; twice < 2; twice++) {
                // Select the side to be sieved
                int side = (firstsieve)?(1-twice):twice;

                // Norm initialization.
                // convention: if a position contains 255, it must stay like
                // this. It means that the other side is hopeless.
                t_norms -= seconds();
                init_norms(S, ffspol[side], I, J, qlat, qlat->side == side,
                        sublat);
                t_norms += seconds();

                // sieve
                t_sieve -= seconds();
                sieveFB(S, FB[side], I, J, sublat);
                t_sieve += seconds();

                t_buckets -= seconds();
                buckets_reset(buckets);
                buckets_fill(buckets, FB[side], sublat, I, J);
                uint8_t *Sptr = S;
                for (unsigned k = 0; k < buckets->n; ++k) {
                    bucket_apply(Sptr, buckets, k);
                    Sptr += bucket_region_size();
                }
                t_buckets += seconds();

                // mark survivors
                // no need to check if this is a valid position
                Sptr = S;
                for (unsigned int k = 0; k < IIJJ; ++k, ++Sptr) {
                    if (*Sptr > threshold[side])
                        *Sptr = 255; 
                }
                // since (0,0) is divisible by everyone, its position might
                // have been clobbered.
                S[0] = 255;
            }
            buckets_clear(buckets);

            t_cofact -= seconds();
            // survivors cofactorization
            {
                fppol_t a, b;
                ij_t i, j, g;
                ij_t hati, hatj;
                fppol_init(a);
                fppol_init(b);

                int rci, rcj = 1;
                for (ij_set_zero(j); rcj; rcj = ij_monic_set_next(j, j, J)) {
                    ijpos_t start = ijvec_get_start_pos(j, I, J);
                    rci = 1;
                    for (ij_set_zero(i); rci; rci = ij_set_next(i, i, I)) {
                        ijpos_t pos = start + ijvec_get_offset(i, I);

                        if (S[pos] != 255) {
#ifdef TRACE_POS
                            if (pos == TRACE_POS) {
                                fprintf(stderr, "TRACE_POS(%d): ", pos);
                                fprintf(stderr,
                                   "entering cofactorization, S[pos] = %d\n",
                                   S[pos]);
                            }
#endif 
                            ij_convert_sublat(hati, hatj, i, j, sublat);
                            ij_gcd(g, hati, hatj);
                            if (ij_deg(g) != 0 && ij_deg(hati)>0  && ij_deg(hatj)>0)
                                continue;
                            ij2ab(a, b, hati, hatj, qlat);
                            nrels += my_factor_survivor(a, b, ffspol, lpb,
                                    qlat->q, qlat->side);
                        }
                    }
                }
                fppol_clear(a);
                fppol_clear(b);
            }
            t_cofact += seconds();

            free(S);
        }  // End of loop on sublattices.

        t_tot = seconds()-t_tot;
        fprintf(stdout, "# Total for this special-q: %d relations found "
                "in %1.1f s\n", nrels, t_tot);
        fprintf(stdout, "# Time of main steps: %1.1f s (initS); "
                "%1.1f s (norms); "
                "%1.1f s (sieve); "
                "%1.1f s (buckets); "
                "%1.1f s (cofact)\n",
                t_initS, t_norms, t_sieve, t_buckets, t_cofact);
        fprintf(stdout, "# Yield: %1.5f s/rel\n", t_tot/nrels);
        tot_nrels += nrels;

    } while (!end_list_of_sq); // End of loop over special-q's
    
    if (sqFile)
        fclose(sqFile);

    tot_time = seconds()-tot_time;
    fprintf(stdout, "###### General statistics ######\n");
    fprintf(stdout, "#   Computed %d special-q\n", tot_sq);
    fprintf(stdout, "#   %d relations found (%1.1f rel/sq)\n",
            tot_nrels, (double)tot_nrels / (double)tot_sq);
    fprintf(stdout, "#   Yield: %1.5f s/rel\n", tot_time/tot_nrels);

    ffspol_clear(ffspol[0]);
    ffspol_clear(ffspol[1]);

    return EXIT_SUCCESS;
}
