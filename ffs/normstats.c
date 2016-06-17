#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <inttypes.h>
#include "macros.h"

#include "types.h"
#include "fppol.h"
#include "ffspol.h"
#include "qlat.h"
#include "norm.h"
#include "timing.h"
#include "ijvec.h"
#include "params.h"
#include "fq.h"
#include "fqpol.h"



void usage(const char *argv0, const char * missing)
{
    fprintf(stderr, "Usage: %s options\n", argv0);
    fprintf(stderr, "  Command line options have the form '-optname optval'\n");
    fprintf(stderr, "List of options (a * means mandatory, default value in []):\n");
    fprintf(stderr, "  pol *            function field polynomial\n");
    fprintf(stderr, "  A      [8]       strict degree bound for a\n");
    fprintf(stderr, "  B      [8]       strict degree bound for b\n");
    fprintf(stderr, "  maxgap [20]      accumulate stats for gaps>maxgap on a single line\n");
    fprintf(stderr, "  gf               indicate the base field for sanity check\n");
 
    if (missing != NULL)
        fprintf(stderr, "Missing parameter: %s\n", missing);
    exit(1);
}

int main(int argc, char **argv)
{
    ffspol_t ffspol[1]; 
    int I=8, J=8;  
    char *argv0 = argv[0];
    int gf = 0;
    int maxgap = 20;

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
        polstr = param_list_lookup_string(pl, "pol");
        if (polstr == NULL) usage(argv0, "pol");
        ffspol_set_str(ffspol[0], polstr);
    }
    // read various bounds
    param_list_parse_int(pl, "A", &I); 
    param_list_parse_int(pl, "B", &J); 
    param_list_parse_int(pl, "maxgap", &maxgap); 
    if (I == 0) usage(argv0, "A");
    if (J == 0) usage(argv0, "B");

    param_list_clear(pl);

    maxgap++;
    {
        double gaps[maxgap];
        for (int i = 0; i < maxgap; ++i)
            gaps[i] = 0;
        double cpt = 0;

        // Compute the norm and the gap everywhere in the (i,j)-plane.
        ij_t i, j, one;
        ij_set_one(one);
        int rci, rcj = 1;
        for(ij_set_zero(j); rcj; rcj = ij_monic_set_next(j, j, J)) {
            rci = 1;
            for(ij_set_zero(i); rci; rci = ij_set_next(i, i, I)) {
                if (ij_is_zero(i) && !ij_eq(j, one))
                    continue;
                if (ij_is_zero(j) && !ij_eq(i, one))
                    continue;

                int gap;
                deg_norm_ij(ffspol[0], i, j, &gap);
                if (gap < 0)
                    gap = 0;
                if (gap >= maxgap)
                    gap = maxgap-1;
                gaps[gap]++;
                cpt++;
            }
        }

        printf("Stats on the gap in the norm for (A, B) = (%d, %d):\n", I, J);
        for(int i = 0; i < maxgap-1; ++i)
            printf("   %3d: %f\n", i, gaps[i]/cpt);
        printf("  >=%2d: %f\n", maxgap-1, gaps[maxgap-1]/cpt);
    } 

    ffspol_clear(ffspol[0]);

    return EXIT_SUCCESS;
}
