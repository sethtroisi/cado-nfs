#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "portability.h"
#include "utils.h"
#include "auxiliary.h" /* for common routines with polyselect.c */
#include "size_optimization.h"
#include "area.h"

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

    poly->skew = L2_skewness (poly->alg, SKEWNESS_DEFAULT_PREC);

    printf ("Initial polynomial:\n");
    if (verbose)
      print_cadopoly_extra (stdout, poly, argc0, argv0, 0);
    else
      printf ("skewness=%1.2f, alpha=%1.2f\n", poly->skew,
              get_alpha (poly->alg, ALPHA_BOUND));
    size_optimization (poly->alg, poly->rat, poly->alg, poly->rat,
                       SOPT_DEFAULT_EFFORT, verbose - 1);
    poly->skew = L2_skewness (poly->alg, SKEWNESS_DEFAULT_PREC);
    
    printf ("After norm optimization:\n");
    if (verbose)
      print_cadopoly_extra (stdout, poly, argc0, argv0, 0);
    else
      printf ("skewness=%1.2f, alpha=%1.2f\n",
              poly->skew, get_alpha (poly->alg, ALPHA_BOUND));

    mpz_set (b, poly->rat->coeff[1]);
    mpz_neg (m, poly->rat->coeff[0]);
    rotate (poly->alg, alim, m, b, &jmin, &kmin, 0, verbose - 1);
    mpz_set (poly->rat->coeff[1], b);
    mpz_neg (poly->rat->coeff[0], m);
    /* optimize again, but only translation */
    mpz_poly_fprintf (stdout, poly->rat);
    sopt_local_descent (poly->alg, poly->rat, poly->alg, poly->rat, 1, -1,
                                          SOPT_DEFAULT_MAX_STEPS, verbose - 1);
    poly->skew = L2_skewness (poly->alg, SKEWNESS_DEFAULT_PREC);
    mpz_clear(b);
    mpz_clear(m);

    print_cadopoly_extra (stdout, poly, argc0, argv0, 0);
    return 0;
} 
