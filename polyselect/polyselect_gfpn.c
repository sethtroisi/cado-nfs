#include "cado.h"
#include <stdio.h>
#include "portability.h"
#include "utils.h"
#include "gfpkdlpolyselect.h"

static void
declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "p", "defines the prime p in GF(p^n)");
    param_list_decl_usage(pl, "n", "defines the degree n in GF(p^n), atm n=2 is mandatory");
    param_list_decl_usage(pl, "ell", "defines the prime order subgroup of GF(p^n), ell divides Phi_n(p)");
    param_list_decl_usage(pl, "mnfs", "defines the number of number fields on second side (\'g\' side or \'linear\' side), should be >= 1, and, for the moment, <= 2");
    param_list_decl_usage(pl, "out",
            "filename where to write selected polynomials");
    verbose_decl_usage(pl);
}

static void
usage (const char *argv, const char * missing, param_list pl)
{
    if (missing) {
        fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
                missing);
    }
    param_list_print_usage(pl, argv, stderr);
    param_list_clear (pl);
    exit (EXIT_FAILURE);
}

int main (int argc, char *argv[])
{
    char **argv0 = argv;
    mpz_t p, ell;
    unsigned int n = 0;
    unsigned int mnfs = 0;

    mpz_init(p);
    mpz_init(ell);

    /* read params */
    param_list pl;
    param_list_init (pl);
    declare_usage(pl);

    if (argc == 1)
        usage (argv0[0], NULL, pl);
    argv++, argc--;
    for ( ; argc; ) {
        if (param_list_update_cmdline (pl, &argc, &argv)) continue;
        fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
        usage (argv0[0], NULL, pl);
    }

    if (!param_list_parse_mpz(pl, "p", p)) // fill mpz_t p with p if a right value is given
        usage(argv0[0], "p", pl);
    if (!param_list_parse_mpz(pl, "ell", ell)) // fill mpz_t ell with ell if a right value is given
        usage(argv0[0], "ell", pl);

    if (!param_list_parse_uint(pl, "n", &n) || (n != 2))
        usage(argv0[0], "n", pl);

    if (!param_list_parse_uint(pl, "mnfs", &mnfs) || (mnfs < 1) || (mnfs > 2))
        usage(argv0[0], "mnfs", pl);

    char *out = param_list_lookup_string (pl, "out");
    if (out == NULL){
      // usage(argv0[0], "out", pl);
      // use a default .poly name: p<n>dd<log_10p>-f<deg f>g<deg g>-<polyselect method>.poly
      //e.g. p2dd31f4g2-conj.poly
      // log_10 (p) i.e. size of p in dd (decimal digits)
      int dd_p = mpz_sizeinbase(p, 10);
      out = (char*)malloc(48);
      sprintf(out, "p%ddd%d-f4g2-conj.poly", n, dd_p);
    }
    verbose_interpret_parameters(pl);

    /* check unused and print command line */
    if (param_list_warn_unused(pl))
        usage (argv0[0], NULL, pl);
    param_list_print_command_line (stdout, pl);

    gfpk_print_params(n,p,ell);

    int return_code = gfpkdlpolyselect( n, p, ell, mnfs, out);

    if (!return_code){
      printf("no polynomial found.\n");
      /* test 1 */
      // build/cormoran/polyselect/polyselect_gfpn -n 2 -p 1000000000000000000000000000099 -ell 127522001020150503761 -mnfs 1 -out superpoly.poly
      // gdb --args build/cormoran/polyselect/polyselect_gfpn -n 2 -p 1000000000000000000000000000099 -ell 127522001020150503761 -mnfs 1 -out superpoly.poly

      /* Start to do something */
      printf("WARNING: we just write a hardcoded polynomial to %s\n", out);
      FILE * outfile = fopen(out, "w");
      if (outfile == NULL) {
        perror("Error while trying to open output file for writing");
        return EXIT_FAILURE;
      }
      fprintf(outfile, "n: 100000000000000000039\n");
      fprintf(outfile, "skew: 1\n");
      fprintf(outfile, "poly0: 1,-35,-48,-35,1\n");
      fprintf(outfile, "poly1: 8385962152,-2002641011,8385962152\n");
      fclose(outfile);
    }

    mpz_clear(p);
    mpz_clear(ell);
    param_list_clear (pl);
    return EXIT_SUCCESS;
}
