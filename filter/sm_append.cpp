/* Shirokauer maps 

   This is largely copied from sm_simple.cpp ; however the purpose of
   sm_simple was to be kept simple and stupid, so let's keep it like
   this.

   This program is to be used as an mpi accelerator for reconstructlog-dl

   Given a relation file where each line begins with an a,b pair IN
   HEXADECIMAL FORMAT, output the same line, appending the SM values in
   the end.

   SM computation is offloaded to the (many) MPI jobs.

   This program works exclusively as a filter (stdin --> stdout).

*/

#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>

#include "macros.h"
#include "utils.h"
#include "relation.h"

static void my_sm(sm_side_info *sm_info, int nb_polys)
{
    char buf[1024];
    mpz_poly_t pol, smpol;
    int maxdeg = sm_info[0]->f->deg;
    for(int side = 1; side < nb_polys; side++)
        maxdeg = MAX(maxdeg, sm_info[side]->f->deg);
    mpz_poly_init(pol, maxdeg);
    mpz_poly_init(smpol, maxdeg);
    while (fgets(buf, 1024, stdin)) {
        int n = strlen(buf);
        if (!n) break;
        buf[n-1]='\0';

        if (buf[0] == '#') {
            puts(buf);
            continue;
        }

        char * p = buf;
        int64_t a;		/* only a is allowed to be negative */
        uint64_t b;
        int64_t sign = 1;
        if (*p=='-') {
            sign=-1;
            p++;
        }
        if (sscanf(p, "%" SCNx64 ",%" SCNx64 ":", &a, &b) < 2) {
            fprintf(stderr, "Parse error at line: %s\n", buf);
            exit(EXIT_FAILURE);
        }

        mpz_poly_init_set_ab(pol, a*sign, b);

        fputs(buf, stdout);
        putchar(':');
        for (int side = 0; side < nb_polys; ++side) {
            compute_sm_piecewise(smpol, pol, sm_info[side]);
            print_sm2(stdout, smpol, sm_info[side]->nsm, sm_info[side]->f->deg, ",");
            if (side == 0 && sm_info[0]->nsm > 0 && sm_info[1]->nsm > 0)
                putchar(',');
        }
        putchar('\n');
    }
}

static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "(required) poly file");
    param_list_decl_usage(pl, "ell", "(required) group order");
    param_list_decl_usage(pl, "nsm", "number of SMs to use per side");
    verbose_decl_usage(pl);
}

static void usage (const char *argv, const char * missing, param_list pl)
{
    if (missing) {
        fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
                missing);
    }
    param_list_print_usage(pl, argv, stderr);
    exit (EXIT_FAILURE);
}

/* -------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
    char *argv0 = argv[0];

    const char *polyfile = NULL;

    param_list pl;
    cado_poly pol;
    mpz_poly_ptr F[NB_POLYS_MAX];

    mpz_t ell;
    double t0;

    /* read params */
    param_list_init(pl);
    declare_usage(pl);

    if (argc == 1)
        usage (argv[0], NULL, pl);

    argc--,argv++;
    for ( ; argc ; ) {
        if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
        usage (argv0, NULL, pl);
    }

    /* Read poly filename from command line */
    if ((polyfile = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: parameter -poly is mandatory\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    /* Read ell from command line (assuming radix 10) */
    mpz_init (ell);
    if (!param_list_parse_mpz(pl, "ell", ell)) {
        fprintf(stderr, "Error: parameter -ell is mandatory\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    /* Init polynomial */
    cado_poly_init (pol);
    cado_poly_read(pol, polyfile);
    for(int side = 0; side < pol->nb_polys; side++)
        F[side] = pol->pols[side];

    int nsm_arg[NB_POLYS_MAX];
    for(int side = 0; side < pol->nb_polys; side++)
        nsm_arg[side]=-1;

    param_list_parse_int_list (pl, "nsm", nsm_arg, pol->nb_polys, ",");

    if (param_list_warn_unused(pl))
        usage (argv0, NULL, pl);

    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);

    sm_side_info sm_info[NB_POLYS_MAX];

    for(int side = 0 ; side < pol->nb_polys; side++) {
        sm_side_info_init(sm_info[side], F[side], ell);
        if (nsm_arg[side] >= 0)
            sm_info[side]->nsm = nsm_arg[side]; /* command line wins */
        printf("# Will use %d SMs on side %d\n", sm_info[side]->nsm, side);
    }

    /*
       for (int side = 0; side < pol->nb_polys; side++) {
       printf("\n# Polynomial on side %d:\nF[%d] = ", side, side);
       mpz_poly_fprintf(stdout, F[side]);

       printf("# SM info on side %d:\n", side);
       sm_side_info_print(stdout, sm_info[side]);

       fflush(stdout);
       }
       */

    t0 = seconds();

    my_sm(sm_info, pol->nb_polys);

    printf("\n# sm completed in %2.2lf seconds\n", seconds() - t0);
    fflush(stdout);

    for(int side = 0 ; side < pol->nb_polys ; side++) {
        sm_side_info_clear(sm_info[side]);
    }

    mpz_clear(ell);
    cado_poly_clear(pol);
    param_list_clear(pl);

    return 0;
}
