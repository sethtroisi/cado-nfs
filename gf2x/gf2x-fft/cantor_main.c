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
    unsigned long *f, *g, *h;

    if (argc != 2) {
	fprintf(stderr, "usage: %s N\n", argv[0]);
	fprintf(stderr,
		"  where N is the number of 64-bit limbs of operands\n");
	exit(1);
    }
    N = atoi(argv[1]);

    f = (unsigned long *) malloc(N * sizeof(unsigned long));
    g = (unsigned long *) malloc(N * sizeof(unsigned long));
    h = (unsigned long *) malloc(2 * N * sizeof(unsigned long));

    __gmp_rands_initialized = 1;
    gmp_randinit_default (__gmp_rands);
    // gmp_randseed_ui(__gmp_rands, time(NULL));
    gmp_randseed_ui(__gmp_rands, 0); //time(NULL));

    mpn_random((mp_limb_t *) f, (sizeof(unsigned long) / sizeof(mp_limb_t)) * N);
    mpn_random((mp_limb_t *) g, (sizeof(unsigned long) / sizeof(mp_limb_t)) * N);

    printf("w := %u;\n", GMP_LIMB_BITS);
    printf("f := [\n");
    for (i = 0; i < N - 1; ++i)
	printf("%lu, ", f[i]);
    printf("%lu\n];\n", f[N - 1]);
    printf("\n");
    printf("g := [\n");
    for (i = 0; i < N - 1; ++i)
	printf("%lu, ", g[i]);
    printf("%lu\n];\n", g[N - 1]);
    printf("\n");

    //for (i = 0; i < 10; ++i) 
    mulCantor128(h,f,N,g,N);

    printf("fg := [\n");
    for (i = 0; i < 2 * N - 1; ++i)
	printf("%lu, ", h[i]);
    printf("%lu\n];\n", h[2 * N - 1]);

    free(h);
    free(g);
    free(f);
    return 0;
}
