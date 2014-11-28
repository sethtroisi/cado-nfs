#include "cado.h"
#include <inttypes.h>
#include <stdlib.h>
#include <gmp.h>
#include "gf2x.h"
#include "gf2x-fft.h"

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "usage: %s N\n", argv[0]);
        fprintf(stderr,
                "  where N is the number of limbs of operands\n");
        exit(1);
    }
    int N = atoi(argv[1]);
    
    unsigned long *f, *g, *h1, *h2;
    f = (unsigned long *) malloc(N * sizeof(unsigned long));
    g = (unsigned long *) malloc(N * sizeof(unsigned long));
    h1 = (unsigned long *) malloc(2 * N * sizeof(unsigned long));
    h2 = (unsigned long *) malloc(2 * N * sizeof(unsigned long));
    ASSERT_ALWAYS(f != NULL);
    ASSERT_ALWAYS(g != NULL);
    ASSERT_ALWAYS(h1 != NULL);
    ASSERT_ALWAYS(h2 != NULL);

    gmp_randstate_t state;
    gmp_randinit_default(state);
    for (int i = 0; i < N; ++i) {
        f[i] = gmp_urandomb_ui(state, 64);
        g[i] = gmp_urandomb_ui(state, 64);
    }

    mulCantor128(h1,f,N,g,N);
    gf2x_mul(h2,f,N,g,N);
    for (int i = 0; i < 2*N; ++i) {
        ASSERT_ALWAYS(h1[i] == h2[i]);
    }
    return EXIT_SUCCESS;
}
