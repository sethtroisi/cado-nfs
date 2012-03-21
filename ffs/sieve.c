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
#include "cofactor.hh"
#include "timing.h"
#include "ijvec.h"
#include "params.h"
#include "sublat.h"


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
    // read q, rho
    {
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
            noerr = factor_base_init(FB[i], filename, fbb[i]);
            if (!noerr) {
                fprintf(stderr, "Could not read %s: %s\n", param, filename);
                exit(EXIT_FAILURE);
            }
        }
    }

    param_list_print_command_line(stderr, pl);
    param_list_clear(pl);

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

    qlat->side = 0; // assume that the special-q is on the algebraic side.

    // Reduce the q-lattice
    int noerr = skewGauss(qlat, 0);
    assert (noerr);
    print_qlat_info(qlat);

    // Precompute lambda for each element of the factor bases.
    for (int i = 0; i < 2; ++i)
        factor_base_precomp_lambda(FB[i], qlat);

    double t_norms = 0;
    double t_sieve = 0;
    double t_cofact = 0;
    int nrels = 0;

    // Loop on all sublattices
    // In the no_sublat case, this loops degenerates into one pass, since
    // nb = 1.
    for (sublat->n = 0; sublat->n < sublat->nb; sublat->n++) {
        if (use_sublat(sublat)) {
            fprintf(stderr, "# Sublattice (");
            fppol16_out(stderr, sublat->lat[sublat->n][0]);
            fprintf(stderr, ", ");
            fppol16_out(stderr, sublat->lat[sublat->n][1]);
            fprintf(stderr, ") :\n");
        }

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

        for (int twice = 0; twice < 2; twice++) {
            // This variable allows to change in which order we handle
            // the polynomials.
            // Put side = 1 - twice, to do the "rational" side first.
            int side = twice;

            // Norm initialization.
            // convention: if a position contains 255, it must stay like
            // this. It means that the other side is hopeless.
            t_norms -= seconds();
            init_norms(S, ffspol[side], I, J, qlat, qlat->side == side, sublat);
            t_norms += seconds();

            // sieve
            t_sieve -= seconds();
            sieveFB(S, FB[side], I, J, sublat);
            t_sieve += seconds();

            // mark survivors
            // no need to check if this is a valid position
            uint8_t *Sptr = S;
            for (unsigned int k = 0; k < IIJJ; ++k, ++Sptr) {
                if (*Sptr > threshold[side])
                    *Sptr = 255; 
            }
        }

        t_cofact -= seconds();
        // survivors cofactorization
        {
            fppol_t a, b;
            ij_t i, j, g;
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
                        ij_gcd(g, i, j);
                        if (ij_deg(g) != 0 && ij_deg(i)>0  && ij_deg(j)>0)
                            continue;
                        ij2ab(a, b, i, j, qlat);
                        nrels += factor_survivor(a, b, ffspol, lpb);
                    }
                }
            }
            fppol_clear(a);
            fppol_clear(b);
        }
        t_cofact += seconds();
        
        free(S);
    }  // End of loop on sublattices.

    fprintf(stdout, "# Total: %d relations found\n", nrels);
    fprintf(stdout, "# Time spent: %1.1f s (norms); %1.1f s"
            " (sieve); %1.1f s (cofact)\n", t_norms, t_sieve, t_cofact);
    fprintf(stdout, "# Rate: %1.5f s/rel\n", (t_norms+t_sieve+t_cofact)/nrels);

    ffspol_clear(ffspol[0]);
    ffspol_clear(ffspol[1]);

    return EXIT_SUCCESS;
}
