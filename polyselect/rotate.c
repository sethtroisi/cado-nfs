#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "portability.h"
#include "utils.h"
#include "auxiliary.h" /* for common routines with polyselect.c */

double bound_f = BOUND_F, bound_g = BOUND_G, area = AREA;

static void usage_and_die(char *argv0) {
    fprintf(stderr, "usage: %s [-area a] [-Bf b] [-Bg c] poly j k\n", argv0);
    fprintf(stderr, "  apply rotation f += (j*x+k)*g to poly.\n");
    fprintf(stderr, "  poly: filename of polynomial\n");
    fprintf(stderr, "  j,k : integers\n");
    exit(1);
}

int main(int argc, char **argv) {
    cado_poly poly;
    long j, k;
    mpz_t b, m;

    mpz_init(b);
    mpz_init(m);
    while (argc >= 2 && argv[1][0] == '-')
      {
        if (strcmp (argv[1], "-area") == 0)
          {
            area = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-Bf") == 0)
          {
            bound_f = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-Bg") == 0)
          {
            bound_g = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else
          break;
      }
    if (argc != 4)
        usage_and_die(argv[0]);
    cado_poly_init (poly);
    if (!cado_poly_read(poly, argv[1])) {
        fprintf(stderr, "Problem when reading file %s\n", argv[1]);
        usage_and_die(argv[0]);
    }
    j = strtol(argv[2], NULL, 10);
    k = strtol(argv[3], NULL, 10);

    mpz_set (b, poly->rat->coeff[1]);
    mpz_neg (m, poly->rat->coeff[0]);
    rotate_aux (poly->alg->coeff, b, m, 0, k, 0);
    rotate_aux (poly->alg->coeff, b, m, 0, j, 1);
    mpz_clear(b);
    mpz_clear(m);

    print_cadopoly_extra (stdout, poly, argc, argv, 0);
    return 0;
} 
