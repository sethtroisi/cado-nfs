#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <ctime>
#include <gmp.h>
#include "utils.h"
#include "gf2x.h"
#include "cantor128.h"
#include "fake_fft.h"
#include "fft_adapter.hpp"

#include "lingen_mat_types.hpp"

void usage()
{
    fprintf(stderr, "Usage: ./bench_polmatmul [--nrep <k>] <deg> <dim>\n");
    exit(1);
}

/* handy globals */
unsigned int nrep = 10;

void bench_polmul(unsigned int d1, unsigned int d2)
{
    unsigned int nw1 = BITS_TO_WORDS(d1 + 1, ULONG_BITS);
    unsigned int nw2 = BITS_TO_WORDS(d2 + 1, ULONG_BITS);
    unsigned long * f = new unsigned long[nw1];
    unsigned long * g = new unsigned long[nw2];
    unsigned long * h = new unsigned long[nw1 + nw2];
    int tt;

    mpn_random((mp_limb_t *) f, nw1);
    mpn_random((mp_limb_t *) g, nw2);

    printf("gf2x_mul (%u x %u):", d1, d2);

    tt = cputime();
    for(unsigned int r = 0 ; r < nrep ; r++) {
        gf2x_mul(h, f, nw1, g, nw2);
    }
    printf(" %f", (double) (cputime()-tt)/nrep);
    printf("\n");

    delete[] f;
    delete[] g;
    delete[] h;
}
void bench_c128(unsigned int d1, unsigned int d2)
{
    unsigned int nw1 = BITS_TO_WORDS(d1 + 1, ULONG_BITS);
    unsigned int nw2 = BITS_TO_WORDS(d2 + 1, ULONG_BITS);
    unsigned long * f = new unsigned long[nw1];
    unsigned long * g = new unsigned long[nw2];
    unsigned long * h = new unsigned long[nw1 + nw2];
    int tt;
    int tt0;

    mpn_random((mp_limb_t *) f, nw1);
    mpn_random((mp_limb_t *) g, nw2);

    c128_info_t o;

    c128_setup(o, d1, d2);

    c128_t tf = c128_alloc(o, 1);
    c128_t tg = c128_alloc(o, 1);
    c128_t th = c128_alloc(o, 1);

    printf("c128 (%u x %u):", d1, d2);

    tt0 = tt = cputime();
    for(unsigned int r = 0 ; r < nrep ; r++) {
        c128_dft(o, tf, f, d1);
    }
    printf(" %f", (double) (cputime()-tt)/nrep);

    tt = cputime();
    for(unsigned int r = 0 ; r < nrep ; r++) {
        c128_dft(o, tg, g, d2);
    }
    printf(" %f", (double) (cputime()-tt)/nrep);

    tt = cputime();
    for(unsigned int r = 0 ; r < nrep ; r++) {
        c128_compose(o, th, tf, tg);
    }
    printf(" %f", (double) (cputime()-tt)/nrep);

    tt = cputime();
    for(unsigned int r = 0 ; r < nrep ; r++) {
        c128_ift(o, h, d1 + d2, th);
    }
    printf(" %f", (double) (cputime()-tt)/nrep);

    printf(" tot %f", (double) (cputime()-tt0)/nrep);
    printf("\n");
    delete[] f;
    delete[] g;
    delete[] h;
}

template<typename T>
void bench_polmatmul(const char * s,
        unsigned int n1, unsigned int n2,
        unsigned int d1, unsigned int d2)
{
    unsigned int nc1 = d1 + 1;
    unsigned int nc2 = d2 + 1;
    unsigned int nw1 = BITS_TO_WORDS(d1 + 1, ULONG_BITS);
    unsigned int nw2 = BITS_TO_WORDS(d2 + 1, ULONG_BITS);

    printf("polmatmul<%s> (dim %u x %u deg %u x %u)", s, n1, n2, d1, d2);

    polmat f(n1, n2, nc1);
    polmat g(n2, n2, nc2);
    polmat h(n1, n2, nc1 + nc2 - 1);

    int tt, tt0;

    for(unsigned int i = 0 ; i < n1 ; i++) {
        for(unsigned int j = 0 ; j < n2 ; j++) {
            mpn_random((mp_limb_t *) f.poly(i,j), nw1);
        }
    }

    for(unsigned int i = 0 ; i < n2 ; i++) {
        for(unsigned int j = 0 ; j < n2 ; j++) {
            mpn_random((mp_limb_t *) g.poly(i,j), nw2);
        }
    }

    T o(d1, d2);

    tpolmat<T> tf(n1,n2,o);
    tpolmat<T> tg(n2,n2,o);
    tpolmat<T> th(n1,n2,o);

    tt0 = tt = cputime();
    for(unsigned int r = 0 ; r < nrep ; r++) {
        transform(tf, f, o, d1);
    }
    printf(" %f", (double) (cputime()-tt)/nrep);

    tt = cputime();
    for(unsigned int r = 0 ; r < nrep ; r++) {
        transform(tg, g, o, d2);
    }
    printf(" %f", (double) (cputime()-tt)/nrep);

    tt = cputime();
    for(unsigned int r = 0 ; r < nrep ; r++) {
        compose(th, tf, tg, o);
    }
    printf(" %f", (double) (cputime()-tt)/nrep);

    tt = cputime();
    for(unsigned int r = 0 ; r < nrep ; r++) {
        itransform(h, th, o, d1 + d2);
    }
    printf(" %f", (double) (cputime()-tt)/nrep);

    printf(" tot %f", (double) (cputime()-tt0)/nrep);
    printf("\n");
}

int main(int argc, char * argv[])
{
    unsigned int d = 0;
    unsigned int n = 0;
    argc--,argv++;

    for( ; argc ; ) {
        if (strcmp(argv[0],"--nrep") == 0) {
            argc--,argv++;
            if (!argc) usage();
            nrep = atoi(argv[0]);
            argc--,argv++;
            continue;
        }
        if (d == 0) {
            d = atoi(argv[0]);
            argc--,argv++;
            continue;
        }
        if (n == 0) {
            n = atoi(argv[0]);
            argc--,argv++;
            continue;
        }
        usage();
    }
    if (d == 0 || n == 0)
        usage();



    /* Typical situation when both blocking parameters are equal
     *
     * E * pi_left:
     *
     * n x 2n, degree d times 2n x 2n, degree d/4
     *
     * pi_left * pi_right
     *
     * 2n x 2n, degree d/4 times the same.
     *
     * We're timing instead:
     * nxn deg d, homogeneous.
     * n/2xn, d x d/4
     */

    /* This is a bench program. We want the output quick.  */

    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    /* Seed the random state. Ugly at will. */

    extern char             __gmp_rands_initialized;
    extern gmp_randstate_t  __gmp_rands;

    __gmp_rands_initialized = 1;
    gmp_randinit_default (__gmp_rands);
    gmp_randseed_ui(__gmp_rands, time(NULL));


    printf("Degree d is %d (%d words)\n", d, BITS_TO_WORDS(d+1,ULONG_BITS));
    printf("Dimension n is %d\n", n);
    printf("Timings are in milliseconds\n");


    /* Start by benching degree d polynomial multiplication */
    bench_polmul(d, d);

    /* Then time for cantor, separating the different steps */
    bench_c128(d, d);

    /* Now the matrix versions */
    bench_polmatmul<fake_fft>("fake_fft", n, n, d, d);
    bench_polmatmul<cantor_fft>("cantor_fft", n, n, d, d);

    
    /* Now do these again, but for unbalanced computations */

    bench_polmul(d, d/4);
    bench_c128(d, d/4);
    bench_polmatmul<fake_fft>("fake_fft", n/2, n, d, d/4);
    bench_polmatmul<cantor_fft>("cantor_fft", n/2, n, d, d/4);
}


