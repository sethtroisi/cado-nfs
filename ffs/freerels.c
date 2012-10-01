#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "macros.h"

#include "types.h"
#include "fppol.h"
#include "ffspol.h"
#include "timing.h"
#include "params.h"
#include "smoothness.h"
#include "polyfactor.h"
#include "fq.h"
#include "fqpol.h"


int sq_is_split(sq_t * roots, sq_srcptr q, ffspol_srcptr F) {
    fq_info_t Fq;
    fq_info_init(Fq, q);

    fqpol_t f;
    fqpol_init(f);
    fqpol_set_ffspol(f, F, Fq);

    int ret = fqpol_is_split(roots, f, Fq);
    fqpol_clear(f);
    fq_info_clear(Fq);
    return ret;
}


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
    fprintf(stderr, "Usage: %s [options | optionfile] \n", argv0);
    fprintf(stderr, "  Command line options have the form '-optname optval'\n");
    fprintf(stderr, "  and take the form 'optname=optval' in the optionfile\n");
    fprintf(stderr, "List of options (a * means mandatory, default value in []):\n");
    fprintf(stderr, "  pol0 *           function field polynomial on side 0\n");
    fprintf(stderr, "  pol1 *           function field polynomial on side 1\n");
    fprintf(stderr, "  lpb0 *           large prime bound on side 0\n");
    fprintf(stderr, "  lpb1 *           large prime bound on side 1\n");
    fprintf(stderr, "  gf               indicate the base field for sanity check\n");
 
    if (missing != NULL)
        fprintf(stderr, "Missing parameter: %s\n", missing);
    exit(1);
}

int main(int argc, char **argv)
{
    ffspol_t ffspol[2]; 
    int lpb[2] = {0, 0};  
    char *argv0 = argv[0];
    int gf = 0;
    
    param_list pl;
    param_list_init(pl);
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
    // Update parameter list at least once to register argc/argv pointers.
    param_list_update_cmdline(pl, &argc, &argv);
    param_list_print_command_line(stdout, pl);

    param_list_parse_int(pl, "gf", &gf);
    if (gf) {
        if (gf != FP_SIZE) {
            fprintf(stderr, "Error: base field mismatch.\n");
            fprintf(stderr, "  The binary is compiled for GF(%d)\n", FP_SIZE);
            fprintf(stderr, "  The parameters are for GF(%d)\n", gf);
            exit(EXIT_FAILURE);
        }
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
    param_list_parse_int(pl, "lpb0", &lpb[0]);
    param_list_parse_int(pl, "lpb1", &lpb[1]);
    if (lpb[0] == 0) usage(argv0, "lpb0");
    if (lpb[1] == 0) usage(argv0, "lpb1");

    param_list_clear(pl);

    uint64_t nrel = 0;
    sq_t p;
    sq_set_ti(p, 1);   // start with the polynomial t, always irreducible
#define MAX_FFS_DEG 20
    fq_t roots0[MAX_FFS_DEG], roots1[MAX_FFS_DEG];
    while (sq_deg(p) <= MIN(lpb[0], lpb[1])) {
        if (sq_is_split(roots0, p, ffspol[0])
                && sq_is_split(roots1, p, ffspol[1])) {
            sq_out(stdout, p);
            printf(",0:");
            for (int i = 0; i < ffspol[0]->deg - 1; ++i) {
                sq_out(stdout, roots0[i]);
                printf(",");
            }
            sq_out(stdout, roots0[ffspol[0]->deg-1]);
            printf(":");
            for (int i = 0; i < ffspol[1]->deg - 1; ++i) {
                sq_out(stdout, roots1[i]);
                printf(",");
            }
            sq_out(stdout, roots1[ffspol[1]->deg-1]);
            printf("\n");
            nrel++;
        }
        do {
            sq_monic_set_next(p, p, 64);
        } while (!sq_is_irreducible(p));
    } 

    fprintf(stdin, "# Computed %" PRIu64 " free relations\n", nrel);
    ffspol_clear(ffspol[0]);
    ffspol_clear(ffspol[1]);

    return EXIT_SUCCESS;
}
