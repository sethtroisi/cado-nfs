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
#include "fq.h"
#include "fqpol.h"

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
    // But first, ensure that b is monic (destructively, but who
    // cares...)

    if (!fppol_is_monic(b)) {
        fp_t lc;
        fppol_get_coeff(lc, b, fppol_deg(b));
        fppol_sdiv(b, b, lc);
        fppol_sdiv(a, a, lc);
    }
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


int sq_roots(sq_t * roots, sq_srcptr q, ffspol_srcptr F) {
    fq_info_t Fq;
    fq_info_init(Fq, q);

    fqpol_t f;
    fqpol_init(f);
    fqpol_set_ffspol(f, F, Fq);

    int nr = fqpol_roots(roots, f, Fq);
    fqpol_clear(f);
    fq_info_clear(Fq);
    return nr;
}

int sq_is_irreducible(sq_srcptr p) {
    fppol_t P;
    fppol_init(P);
    fppol_set_sq(P, p);
    int ret = fppol_is_irreducible(P);
    fppol_clear(P);
    return ret;
}

int sq_cmp(sq_srcptr p, sq_srcptr q) {
    if (sq_eq(p, q))
        return 0;
    uint64_t ip, iq;
    ip = sq_get_ui(p, 64, 0);
    iq = sq_get_ui(q, 64, 0);
    if (ip < iq)
        return -1;
    else
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
    fprintf(stderr, "  q0   *           lower bound for special-q range\n");
    fprintf(stderr, "  q1   *           lower bound for special-q range\n");
    fprintf(stderr, "    Note: giving (q0,q1) is exclusive to giving (q,rho). In the latter case,\n" "    rho is optional.\n");
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
    int sqside = 0;
    int firstsieve = 0;
    sq_t q0, q1;
    int rho_given = 0;

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

    // read q0, q1
    {
        const char *sqstr;
        int noerr;
        sqstr = param_list_lookup_string(pl, "q0");
        if (sqstr != NULL) {
            noerr = sq_set_str(q0, sqstr);
            if (!noerr) {
                fprintf(stderr, "Could not parse q0: %s\n", sqstr);
                exit(EXIT_FAILURE);
            }
            if (!sq_is_monic(q0)) {
                fprintf(stderr, "Error: given q0 is not monic: %s\n", sqstr);
                exit(EXIT_FAILURE);
            }
            sqstr = param_list_lookup_string(pl, "q1");
            if (sqstr == NULL) usage(argv0, "q1");
            noerr = sq_set_str(q1, sqstr);
            if (!noerr) {
                fprintf(stderr, "Could not parse q1: %s\n", sqstr);
                exit(EXIT_FAILURE);
            }
            if (sq_deg(q1) > lpb[0]) {
                fprintf(stderr, "WARNING: not a good idea to have special-q beyond the large prime bound!\n");
            }
        } else {
            sq_set_zero(q0);
            sq_set_zero(q1);
        }
    }

    // read q, rho 
    // We store them in qlat for convenience, but these will be moved
    // away before the main loop.
    {
        const char *sqstr;
        int noerr;
        sqstr = param_list_lookup_string(pl, "q");
        if (sq_is_zero(q0) && (sqstr == NULL)) usage(argv0, "q");
        if (sqstr != NULL) {
            if (!sq_is_zero(q0)) {
                fprintf(stderr, "You can not provide both (q0,q1) and q\n");
                exit(EXIT_FAILURE);
            }
            noerr = sq_set_str(qlat->q, sqstr);
            if (!noerr) {
                fprintf(stderr, "Could not parse q: %s\n", sqstr);
                exit(EXIT_FAILURE);
            }
            if (!sq_is_irreducible(qlat->q)) {
                fprintf(stderr, "Error, q is not irreducible: %s\n", sqstr);
                exit(EXIT_FAILURE);
            }
            if (!sq_is_monic(qlat->q)) {
                fprintf(stderr, "Error, q is not monic: %s\n", sqstr);
                exit(EXIT_FAILURE);
            }
            sqstr = param_list_lookup_string(pl, "rho");
            if (sqstr != NULL) {
                rho_given = 1;
                noerr = sq_set_str(qlat->rho, sqstr);
                if (!noerr) {
                    fprintf(stderr, "Could not parse rho: %s\n", sqstr);
                    exit(EXIT_FAILURE);
                }
                if (sq_deg(qlat->q) > lpb[0]) {
                    fprintf(stderr, "WARNING: not a good idea to have a special-q beyond the large prime bound!\n");
                }
            }
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
            noerr = factor_base_init(FB[i], filename, I, fbb[i]);
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

    // Allocate storage space for the buckets.
    // FIXME: The bucket capacity is hardcoded for the moment.
    buckets_t buckets;
    buckets_init(buckets, I, J, 1<<16, I,
                 1+MAX(factor_base_max_degp(FB[0]),
                       factor_base_max_degp(FB[1])));
    print_bucket_info(buckets);

    qlat->side = sqside; 

    double tot_time = seconds();
    int tot_nrels = 0;
    int tot_sq = 0;
    
    sq_t * roots;
    roots = (sq_t *) malloc (ffspol[sqside]->deg * sizeof(sq_t));
    ASSERT_ALWAYS(roots != NULL);
    int nroots = 0; // number of roots still to work on for current q.

    if (sq_is_zero(q0)) {
        sq_set(q0, qlat->q);
        sq_set(q1, qlat->q);
    }
    if (rho_given) {
        nroots = 1;
        sq_set(roots[0], qlat->rho);
    } else {
        if (sq_is_irreducible(q0)) {
            nroots = sq_roots(roots, q0, ffspol[sqside]);
            sq_set(qlat->q, q0);
            printf("############################################\n");
            printf("# Roots for q = "); 
            sq_out(stdout, q0);
            printf(":");
            for (int i = 0; i < nroots; ++i) {
                printf(" ");
                sq_out(stdout, roots[i]);
            }
            printf("\n");
        }
    }


    //
    // Begin of loop over special-q's
    //
    do {
        // Select next special-q
        if (nroots == 0) { // find next q
            do {
                do {
                    sq_monic_set_next(q0, q0, 64);
                } while (!sq_is_irreducible(q0));
                nroots = sq_roots(roots, q0, ffspol[sqside]);
            } while (nroots == 0);

            if (sq_cmp(q0, q1) >= 0)
                break;

            printf("############################################\n");
            printf("# Roots for q = "); 
            sq_out(stdout, q0);
            printf(":");
            for (int i = 0; i < nroots; ++i) {
                printf(" ");
                sq_out(stdout, roots[i]);
            }
            printf("\n");

            sq_set(qlat->q, q0);
            sq_set(qlat->rho, roots[nroots-1]);
            nroots--;
        } else {
            sq_set(qlat->rho, roots[nroots-1]);
            nroots--;
        }

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
        printf("############################################\n");
        print_qlat_info(qlat);

        double t_norms = 0;
        double t_sieve = 0;
        double t_buckets = 0;
        double t_cofact = 0;
        double t_initS = 0;
        double t_lambda = 0;
        int nrels = 0;

        // Precompute lambda for each element of the factor bases.
        t_lambda -= seconds();
        for (int i = 0; i < 2; ++i) 
            factor_base_precomp_lambda(FB[i], qlat, sublat);
        t_lambda += seconds();


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

            for (int twice = 0; twice < 2; twice++) {
                // Select the side to be sieved
                int side = (firstsieve)?(1-twice):twice;

                // Fill the buckets.
                t_buckets -= seconds();
                buckets_fill(buckets, FB[side], sublat, I, J);
                t_buckets += seconds();

                uint8_t *Sptr = S;
                ij_t j;
                ij_set_zero(j);
                for (unsigned k = 0, pos = 0; k < buckets->n; ++k) {
                  // Norm initialization.
                  // convention: if a position contains 255, it must stay like
                  // this. It means that the other side is hopeless.
                  t_norms -= seconds();
                  init_norms(Sptr, ffspol[side], I, J, j, pos,
                             bucket_region_size(),
                             qlat, qlat->side == side, sublat);
                  t_norms += seconds();

                  Sptr += bucket_region_size();
                  pos  += bucket_region_size();
                  while (ij_monic_set_next(j, j, J) &&
                         ijvec_get_start_pos(j, I, J) < pos);
                }

                // sieve
                t_sieve -= seconds();
                sieveFB(S, FB[side], I, J, sublat);
                t_sieve += seconds();

                Sptr = S;
                for (unsigned k = 0; k < buckets->n; ++k) {
                  // Apply the updates from the corresponding bucket.
                  t_buckets -= seconds();
                  bucket_apply(Sptr, buckets, k);
                  t_buckets += seconds();

                  // since (0,0) is divisible by everyone, its position might
                  // have been clobbered.
                  if (!k) *Sptr = 255;

                  // mark survivors
                  // no need to check if this is a valid position
                  for (unsigned i = 0; i < bucket_region_size(); ++i, ++Sptr) {
                    if (*Sptr > threshold[side])
                      *Sptr = 255; 
                  }
                }
            }

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
        fprintf(stdout,
                "# Time of main steps: "
                "%1.1f s (lambda); "
                "%1.1f s (initS);   "
                "%1.1f s (norms);\n"
                "#                     "
                "%1.1f s (sieve);  "
                "%1.1f s (buckets); "
                "%1.1f s (cofact).\n",
                t_lambda, t_initS, t_norms, t_sieve, t_buckets, t_cofact);
        fprintf(stdout, "# Yield: %1.5f s/rel\n", t_tot/nrels);
        tot_nrels += nrels;

    } while (1); // End of loop over special-q's

    factor_base_clear(FB[0]);
    factor_base_clear(FB[1]);
    buckets_clear(buckets);
    free(roots);

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
