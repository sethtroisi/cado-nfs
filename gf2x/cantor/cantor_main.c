#include <inttypes.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>

#include "cantor128.h"

/* Be ugly */
extern char             __gmp_rands_initialized;
extern gmp_randstate_t  __gmp_rands;

int main(int argc, char **argv)
{
    int i, N;
    uint64_t *f, *g, *h;

    if (argc != 2) {
	fprintf(stderr, "usage: %s N\n", argv[0]);
	fprintf(stderr,
		"  where N is the number of 64-bit limbs of operands\n");
	exit(1);
    }
    N = atoi(argv[1]);

    f = (uint64_t *) malloc(N * sizeof(uint64_t));
    g = (uint64_t *) malloc(N * sizeof(uint64_t));
    h = (uint64_t *) malloc(2 * N * sizeof(uint64_t));

    __gmp_rands_initialized = 1;
    gmp_randinit_mt (__gmp_rands);
    gmp_randseed_ui(__gmp_rands, time(NULL));

    mpn_random((mp_limb_t *) f, (sizeof(uint64_t) / sizeof(mp_limb_t)) * N);
    mpn_random((mp_limb_t *) g, (sizeof(uint64_t) / sizeof(mp_limb_t)) * N);

    printf("f := [\n");
    for (i = 0; i < N - 1; ++i)
	printf("%" PRIu64 ", ", f[i]);
    printf("%" PRIu64 "\n];\n", f[N - 1]);
    printf("\n");
    printf("g := [\n");
    for (i = 0; i < N - 1; ++i)
	printf("%" PRIu64 ", ", g[i]);
    printf("%" PRIu64 "\n];\n", g[N - 1]);
    printf("\n");

    //for (i = 0; i < 10; ++i) 
    mulCantor128(h, f, N, g, N);

    printf("fg := [\n");
    for (i = 0; i < 2 * N - 1; ++i)
	printf("%" PRIu64 ", ", h[i]);
    printf("%" PRIu64 "\n];\n", h[2 * N - 1]);

    free(h);
    free(g);
    free(f);
    return 0;
}
