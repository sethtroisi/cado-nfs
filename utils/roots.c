/* Usage: roots -poly xxx.poly -q <q> [-side <side>]
   prints the roots of the polynomial mod q
 */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "modul_poly.h"
#include <gmp.h>
#include "utils.h"
#include "macros.h"
#include "portability.h"

void usage(const char *argv0)
{
    fprintf(stderr, "Usage: %s -poly xxx.poly -q <q> [-side <side>]\n", argv0);
    fprintf(stderr, "prints the roots of polynomial[side] mod q");
    fprintf(stderr, " (default side is ALGEBRAIC)\n");
    exit(1);
}

int main(int argc0, char *argv0[])
{
    const char *polyfilename = NULL;
    cado_poly pol;
    param_list pl;
    uint64_t q = 0;
    uint64_t q0 = 0;
    uint64_t q1 = 0;
    int argc = argc0;
    char **argv = argv0;
    int side = ALGEBRAIC_SIDE;

    param_list_init(pl);
    argv++, argc--;
    for (; argc;) {
	if (param_list_update_cmdline(pl, &argc, &argv)) {
	    continue;
	}
	/* Could also be a file */
	FILE *f;
	if ((f = fopen(argv[0], "r")) != NULL) {
	    param_list_read_stream(pl, f);
	    fclose(f);
	    argv++, argc--;
	    continue;
	}
	fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
	usage(argv0[0]);
    }

    param_list_parse_int(pl, "side", &side);
    param_list_parse_uint64(pl, "q", &q);
    param_list_parse_uint64(pl, "q0", &q0);
    param_list_parse_uint64(pl, "q1", &q1);

    if (!((q != 0) ^ (q0 != 0 && q1 != 0))) {
	usage(argv0[0]);
    }

    cado_poly_init(pol);
    ASSERT_ALWAYS((polyfilename =
		   param_list_lookup_string(pl, "poly")) != NULL);
    if (!cado_poly_read(pol, polyfilename)) {
	fprintf(stderr, "Error reading polynomial file\n");
	usage(argv[0]);
	exit(EXIT_FAILURE);
    }
    if (pol->skew <= 0.0) {
	fprintf(stderr, "Error, please provide a positive skewness\n");
	exit(EXIT_FAILURE);
    }

    mpz_poly_ptr ps = pol->pols[side];

    if (q) {
        /* check that q fits in an unsigned long */
        if (q > (uint64_t) ULONG_MAX) {
            fprintf(stderr, "Error, q=%" PRIu64 " exceeds ULONG_MAX\n", q);
            exit(EXIT_FAILURE);
        }

        modulusul_t qq;
        modul_initmod_ul(qq, q);
        residueul_t * r = malloc(ps->deg * sizeof(residueul_t));
        for(int i = 0 ; i < ps->deg ; i++)
            modul_init(r[i], qq);
        int nr = modul_poly_roots(r, ps, qq);
        for(int i = 0 ; i < nr ; i++) {
            printf("%lu\n", modul_get_ul(r[i], qq));
        }
        for(int i = 0 ; i < ps->deg ; i++)
            modul_clear(r[i], qq);
    } else {
        if (q0 > q1)
            usage(argv[0]);
        /* check that q1 fits in an unsigned long */
        if (q1 > (uint64_t) ULONG_MAX) {
            fprintf(stderr, "Error, q=%" PRIu64 " exceeds ULONG_MAX\n", q);
            exit(EXIT_FAILURE);
        }

        mpz_t qz;
        mpz_init(qz);
        mpz_set_ui(qz, q0);

        if (!mpz_probab_prime_p(qz, 10)) {
            mpz_nextprime(qz,qz);
        }
        for( ; mpz_get_ui(qz) < q1 ; mpz_nextprime(qz,qz) ) {
            modulusul_t qq;
            modul_initmod_ul(qq, mpz_get_ui(qz));
            residueul_t * r = malloc(ps->deg * sizeof(residueul_t));
            for(int i = 0 ; i < ps->deg ; i++)
                modul_init(r[i], qq);
            int nr = modul_poly_roots(r, ps, qq);
            if (nr) {
                printf("%lu", mpz_get_ui(qz));
                for(int i = 0 ; i < nr ; i++) {
                    printf(" %lu", modul_get_ul(r[i], qq));
                }
                printf("\n");
            }
            for(int i = 0 ; i < ps->deg ; i++)
                modul_clear(r[i], qq);
        }
    }

    cado_poly_clear(pol);
    param_list_clear(pl);
    return 0;
}
