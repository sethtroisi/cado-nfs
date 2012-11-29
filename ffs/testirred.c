#include <stdio.h>
#include <stdlib.h>

#include "types.h"
#include "ffspol.h"
#include "fqpol.h"
#include "params.h"
#include "polyfactor.h"

// This function is duplicated from sieve.c
int sq_is_irreducible(sq_srcptr p) {
    fppol_t P;
    fppol_init(P);
    fppol_set_sq(P, p);
    int ret = fppol_is_irreducible(P);
    fppol_clear(P);
    return ret;
}

void usage(const char *argv0, const char * missing)
{
    fprintf(stderr, "Usage: %s [options | optionfile] < poly_file\n", argv0);
    fprintf(stderr, "  Command line options have the form '-optname optval'\n");
    fprintf(stderr, "  and take the form 'optname=optval' in the optionfile\n");
    fprintf(stderr, "Polynomials to be tested are given on stdin in sieve4ffs format\n");
    fprintf(stderr, "List of options:\n");
    fprintf(stderr, "  deg_ell          degree of ell used for testing\n");
    fprintf(stderr, "  ntest            number of ell tested per polynomial\n");
    if (missing != NULL)
        fprintf(stderr, "Missing parameter: %s\n", missing);
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
    // Default parameters:
    int DEG_ELL = 11;
    int NTEST = 10;

    // Command-line parameters:
    param_list pl;
    param_list_init(pl);
    char *argv0 = argv[0];
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

    param_list_parse_int(pl, "deg_ell", &DEG_ELL);
    param_list_parse_int(pl, "ntest", &NTEST);

    param_list_update_cmdline(pl, &argc, &argv);
    
    // read function field polynomials from stdin, one per line
    ffspol_t ffspol; 
    ffspol_init(ffspol);
    char line[1000];
    while (1) {
        if (fgets(line, 1000, stdin) == NULL) 
            break;
        if (!ffspol_set_str(ffspol, line))
            break;
        ffspol_out(stdout, ffspol);
        sq_t ell;
        sq_set_ti(ell, DEG_ELL); // degree of polynomials we use for testing
        int irred = 0;
        for (int i = 0; i < NTEST; ++i) {
            // select next irreducible polynomial ell
            do {
                sq_monic_set_next(ell, ell, 64);
            } while (!sq_is_irreducible(ell));

            // build the corresponding finite field and reduce the ffspol
            fq_info_t Fq;
            fq_info_init(Fq, ell);
            fqpol_t Fbar;
            fqpol_init(Fbar);
            fqpol_set_ffspol(Fbar, ffspol, Fq);

            // is it irred mod ell?
            if ((Fbar->deg == ffspol->deg) && (fqpol_is_irreducible(Fbar,Fq))) {
                irred = 1;
                fqpol_clear(Fbar);
                fq_info_clear(Fq);
                break;
            }
            fqpol_clear(Fbar);
            fq_info_clear(Fq);
        }
        if (irred) {
            printf(" irreducible; witness = ");
            sq_out(stdout, ell);
            printf("\n");
        } else {
            printf(" probably reducible\n");
        }
    }
    ffspol_clear(ffspol);
    param_list_clear(pl);
    return EXIT_SUCCESS;
}


