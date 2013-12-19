#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "portability.h"
#include "utils.h"
#include "auxiliary.h" /* for common routines with polyselect.c */

extern int MAX_k;

static void
usage_and_die (char *argv0)
{
  fprintf (stderr, "usage: %s [-v] poly kmax\n", argv0);
  fprintf (stderr, "  apply rotation f += (j*x+k)*g to poly.\n");
  fprintf (stderr, "  poly: filename of polynomial\n");
  fprintf (stderr, "  j,k : integers\n");
  exit (1);
}

int
main (int argc, char **argv)
{
    cado_poly poly;
    long kmax, jmin, kmin;
    unsigned long alim = 2000;
    mpz_t b, m;
    int argc0 = argc, verbose = 0;
    char **argv0 = argv;

    while (argc >= 2 && strcmp (argv[1], "-v") == 0)
      {
        argv ++;
        argc --;
        verbose ++;
      }

    mpz_init(b);
    mpz_init(m);
    if (argc != 3)
        usage_and_die (argv0[0]);
    cado_poly_init (poly);
    if (!cado_poly_read(poly, argv[1])) {
        fprintf(stderr, "Problem when reading file %s\n", argv[1]);
        usage_and_die (argv0[0]);
    }
    kmax = strtol(argv[2], NULL, 10);
    MAX_k = kmax;

    mpz_poly_t F;
    F->coeff = poly->alg->coeff;
    F->deg = poly->alg->deg;
    poly->skew = L2_skewness (F, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);

    printf ("Initial polynomial:\n");
    if (verbose)
      print_cadopoly_extra (stdout, poly, argc0, argv0, 0);
    else
      printf ("skewness=%1.2f, alpha=%1.2f\n", poly->skew,
              get_alpha (F, ALPHA_BOUND));
    optimize (F, poly->rat->coeff, verbose - 1, 1);
    poly->skew = L2_skewness (F, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    
    printf ("After norm optimization:\n");
    if (verbose)
      print_cadopoly_extra (stdout, poly, argc0, argv0, 0);
    else
      printf ("skewness=%1.2f, alpha=%1.2f\n",
              poly->skew, get_alpha (F, ALPHA_BOUND));

    mpz_set (b, poly->rat->coeff[1]);
    mpz_neg (m, poly->rat->coeff[0]);
    rotate (F, alim, m, b, &jmin, &kmin, 0, verbose - 1, DEFAULT_L2_METHOD);
    mpz_set (poly->rat->coeff[1], b);
    mpz_neg (poly->rat->coeff[0], m);
    /* optimize again, but only translation */
    fprint_polynomial (stdout, poly->rat->coeff, poly->rat->deg);
    optimize (F, poly->rat->coeff, verbose - 1, 0);
    poly->skew = L2_skewness (F, SKEWNESS_DEFAULT_PREC, DEFAULT_L2_METHOD);
    mpz_clear(b);
    mpz_clear(m);

    print_cadopoly_extra (stdout, poly, argc0, argv0, 0);
    return 0;
} 
