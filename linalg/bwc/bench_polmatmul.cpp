#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <ctime>
#include <gmp.h>
#include "utils.h"
#include "gf2x.h"
#include "gf2x-fft.h"

#include "lingen_mat_types.hpp"

void usage()
{
    fprintf(stderr, "Usage: ./bench_polmatmul [--nrep <k>] <deg> <m> <n>\n");
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
#define do_and_noprinttime(v, line, nrep, tlim) do {			\
    clock_t tt = clock();						\
    unsigned int r;							\
    clock_t clocklim = tt + tlim * CLOCKS_PER_SEC;			\
    for(r = 0 ; r < nrep && clock() < clocklim ; r++) {			\
        line;								\
    }									\
    v=(double) (clock()-tt)/r/CLOCKS_PER_SEC;           		\
} while (0)

#define do_and_printtime(v, line, nrep, tlim) do {			\
    do_and_noprinttime(v, line, nrep, tlim);                            \
    printf(" %f", v);                                                   \
} while (0)


void bench_c128(unsigned int d1, unsigned int d2)
{
    unsigned int nw1 = BITS_TO_WORDS(d1 + 1, ULONG_BITS);
    unsigned int nw2 = BITS_TO_WORDS(d2 + 1, ULONG_BITS);
    unsigned long * f = new unsigned long[nw1];
    unsigned long * g = new unsigned long[nw2];
    unsigned long * h = new unsigned long[nw1 + nw2];
    int tt0;

    mpn_random((mp_limb_t *) f, nw1);
    mpn_random((mp_limb_t *) g, nw2);

    c128_info_t o;

    c128_init(o, d1, d2);

    c128_t * tf = c128_alloc(o, 1);
    c128_t * tg = c128_alloc(o, 1);
    c128_t * th = c128_alloc(o, 1);

    printf("c128 (%u x %u):", d1, d2);

    tt0 = cputime();
    double v;

    do_and_printtime(v, c128_dft(o, tf, f, d1), nrep, 1);
    do_and_printtime(v, c128_dft(o, tg, g, d2), nrep, 1);
    do_and_printtime(v, c128_compose(o, th, tf, tg), nrep, 1);
    do_and_printtime(v, c128_ift(o, h, d1 + d2, th), nrep, 1);

    c128_clear(o);
    printf(" tot %f", (double) (cputime()-tt0)/nrep);
    printf("\n");
    delete[] f;
    delete[] g;
    delete[] h;
}

/* Since strassen is a recursive algorithm that falls back on matrices
 * whose sizes are halved, and since for our applications, the matrices
 * do have a large power of 2 in their dimension, we consider different
 * possible dimensions of the form (d1*2^i,d2*2^i) times (d2*2^i,d3*2^i),
 * for all values of i from 0 to some bound, and for all value of d from
 * 1 to some power of 2 (of course 0 for d makes no sense).
 */

#define BITS_IN_DIM_D   3
#define BITS_IN_DIM_I   4

#define MAX_I_FOR_TUNING        256    /* < (1 << BITS_IN_DIM_I) */
#define MAX_CONSIDERED_THRESHOLD        1000000
#define MIN_CONSIDERED_THRESHOLD        0
// dimension must be <= d * 2^i
// d <= (1 << BITS_IN_DIM_D)
// i < (1 << BITS_IN_DIM_I)
//
// one tuning table entry contains
// 3*BITS_IN_DIM_D+BITS_IN_DIM_I bits
//
// representing a max dimension (1 << BITS_IN_DIM_D) * (1 << ((1 <<
// BITS_IN_DIM_I) - 1)) ; however not all dimensions can be represented
// up to this size, for obvious reasons.
//
// 1,3 means at most 256        ;  6 bits per entry index
// 2,3 means at most 512        ;  9 bits per entry index
// 1,4 means at most 65536      ;  7 bits per entry index
// 2,4 means at most 131072     ;  10 bits per entry index
// 3,4 means at most 262144     ;  13 bits per entry index

// contains the value n0 such that for nbits >= n0, we'd better use strassen

#define SELECTOR_NB_INDICES (1 << (3*BITS_IN_DIM_D+BITS_IN_DIM_I))

struct my_strassen_selector {

    unsigned int strassen_threshold[SELECTOR_NB_INDICES];

    unsigned int threshold_index(unsigned int m, unsigned int n, unsigned int p) const {
        unsigned int combined = m|n|p;
        unsigned int b = ctzl(combined);
        if (b >= (1UL << BITS_IN_DIM_I)) {
            // arrange so that we 
            b = (1UL << BITS_IN_DIM_I) - 1;
        }
        m >>= b; m-=1; ASSERT(m < (1 << BITS_IN_DIM_D));
        n >>= b; n-=1; ASSERT(n < (1 << BITS_IN_DIM_D));
        p >>= b; p-=1; ASSERT(p < (1 << BITS_IN_DIM_D));

        unsigned int index = 0;
        index = m;
        index = (index << BITS_IN_DIM_D) | n;
        index = (index << BITS_IN_DIM_D) | p;
        index = (index << BITS_IN_DIM_I) | b;
        return index;
    }

    unsigned int const & threshold(unsigned int m, unsigned int n, unsigned int p) const {
        return strassen_threshold[threshold_index(m,n,p)];
    }

    unsigned int& threshold(unsigned int m, unsigned int n, unsigned int p) {
        return strassen_threshold[threshold_index(m,n,p)];
    }

    int operator()(unsigned int m, unsigned int n, unsigned int p, unsigned int nbits) const {
        int answer = nbits >= threshold(m,n,p);
        // printf("([%u,%u,%u,%u]=>%d)",m,n,p,nbits,answer);
        return answer;
    }
};


struct logbook {
    double t1;
    double t2;
    double c;
    double i;
};

template<typename T>
struct foo {
    static logbook l;
    static my_strassen_selector s;
};

template<> logbook foo<fake_fft>::l = logbook();
template<> logbook foo<cantor_fft>::l = logbook();
template<> my_strassen_selector foo<fake_fft>::s = my_strassen_selector();
template<> my_strassen_selector foo<cantor_fft>::s = my_strassen_selector();

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

    do_and_printtime(foo<T>::l.t1, transform(tf, f, o, d1), nrep, 1);
    do_and_printtime(foo<T>::l.t2, transform(tg, g, o, d2), nrep, 1);
    do_and_printtime(foo<T>::l.c, compose(th, tf, tg, o), nrep, 1);
    do_and_printtime(foo<T>::l.i, itransform(h, th, o, d1 + d2), nrep, 1);

    double t = foo<T>::l.t1 + foo<T>::l.t2 + foo<T>::l.c + foo<T>::l.i;

    printf(" tot %f", t);
    printf("\n");
}

template<typename fft_type>
void randomize(tpolmat<fft_type> & t)
{
    for(unsigned int i = 0 ; i < t.nrows ; i++) {
        for(unsigned int j = 0 ; j < t.ncols ; j++) {
            mpn_random((mp_limb_t *) t.poly(i,j),
                    (t.o.size() * sizeof(typename fft_type::t)) / sizeof(mp_limb_t));
        }
    }
}

template<typename fft_type>
void tune_strassen1(const char * name,
        unsigned int n1, unsigned int n2, unsigned int n3)
{
    my_strassen_selector& s(foo<fake_fft>::s);

    printf("[%s] Tuning strassen for %ux%u * %ux%u\n",name,n1,n2,n2,n3);
    unsigned int bp = ctzl(n1 | n2 | n3);
    printf("%u is %u*2^%u\n", n1, n1 >> bp, bp);
    printf("%u is %u*2^%u\n", n2, n2 >> bp, bp);
    printf("%u is %u*2^%u\n", n3, n3 >> bp, bp);
    unsigned int nn1 = n1 >> bp;
    unsigned int nn2 = n2 >> bp;
    unsigned int nn3 = n3 >> bp;
    ASSERT((nn1-1) < (1 << BITS_IN_DIM_D));
    ASSERT((nn2-1) < (1 << BITS_IN_DIM_D));
    ASSERT((nn3-1) < (1 << BITS_IN_DIM_D));
    // for the very small matrices, we just _cannot_ use strassen because
    // there is an odd dimension hanging around.
    s.threshold(nn1,nn2,nn3) = UINT_MAX;
        printf("Threshold for %ux%u * %ux%u is %u bits [index %u]\n",
                nn1, nn2, nn2, nn3, UINT_MAX,
                s.threshold_index(nn1,nn2,nn3));
    nn1 <<= 1;
    nn2 <<= 1;
    nn3 <<= 1;
    // If we know that for smaller matrices, strassen was paying off from
    // nbits on, then of course it is going to pay off at worst at this
    // size.
    // old_wt stands for old threshold in words.
    unsigned int max_wt = iceildiv(MAX_CONSIDERED_THRESHOLD, ULONG_BITS);
    unsigned int min_wt = iceildiv(MIN_CONSIDERED_THRESHOLD, ULONG_BITS);
    unsigned int old_wt = max_wt;
    for( ; nn1 <= n1 ; nn1 <<= 1, nn2 <<= 1, nn3 <<= 1) {
        printf("Tuning for %ux%u * %ux%u\n",nn1, nn2, nn2, nn3);
        if (old_wt == min_wt) {
            printf("[keeping low threshold]\n");
            s.threshold(nn1,nn2,nn3) = 0;
            continue;
        }
        // in order to avoid accumulating errors, we introduce some
        // looseness here. Setting wt1 to twice the old wt value will
        // force re-checking at the previous cutoff value. If ever this
        // cutoff value was wrong, we'll settle one one which sits above
        // in the interval.
        unsigned int wt1 = 2 * old_wt;
        unsigned int wt0 = min_wt;
        for( ; wt1 - wt0 > 1  && wt0 < max_wt; ) {
            // assuming cubic is better for 64 * t0, ans strassen better
            // for 64 * t1
            unsigned int wt = (wt1 + wt0) / 2;
            unsigned int t1 = wt1 * ULONG_BITS;
            unsigned int t0 = wt0 * ULONG_BITS;
            unsigned int t = wt * ULONG_BITS;
            if (t == 0) { t = 2; }
            fft_type o(t / 2, t / 2);
            typedef typename fft_type::t ft;
            unsigned int tr = (o.size() * sizeof(ft)) * CHAR_BIT;
            tpolmat<fft_type> tf(nn1, nn2, o);
            tpolmat<fft_type> tg(nn2, nn3, o);
            tpolmat<fft_type> th(nn1, nn3, o);

            randomize(tf);
            randomize(tg);
            printf("%u [%u]", t, tr);
            double t_cubic;
            double t_strassen;
            // force cubic
            s.threshold(nn1,nn2,nn3) = t1;
            do_and_printtime(t_cubic, compose_inner(th, tf, tg, o, s), nrep, 1);
            // force strassen
            s.threshold(nn1,nn2,nn3) = t0;
            do_and_printtime(t_strassen, compose_inner(th, tf, tg, o, s), nrep, 1);
            if (t_cubic <= t_strassen) {
                printf(" cubic [%.2f]\n", t_cubic / t_strassen);
                wt0 = wt;
            } else {
                printf(" strassen [%.2f]\n", t_cubic / t_strassen);
                wt1 = wt;
            }
        }
        if (wt1 > max_wt)
            wt1 = max_wt;
        old_wt = wt1;
        if (wt1 == max_wt) wt1 = UINT_MAX;
        if (wt1 == min_wt) wt1 = 0;
        printf("Threshold for %ux%u * %ux%u is %u bits [index %u]\n",
                nn1, nn2, nn2, nn3, wt1 * ULONG_BITS,
                s.threshold_index(nn1,nn2,nn3));
        s.threshold(nn1,nn2,nn3) = wt1 * ULONG_BITS;
    }
}

template<typename fft_type>
void tune_strassen(const char * name)
{
    my_strassen_selector& s(foo<fake_fft>::s);
    for(unsigned int i = 1 ; i <= (1 << BITS_IN_DIM_D) ; i++) {
    for(unsigned int j = 1 ; j <= (1 << BITS_IN_DIM_D) ; j++) {
    for(unsigned int k = 1 ; k <= (1 << BITS_IN_DIM_D) ; k++) {
        tune_strassen1<fft_type>(name,
                i * MAX_I_FOR_TUNING,
                j * MAX_I_FOR_TUNING,
                k * MAX_I_FOR_TUNING);
    }
    }
    }
    for(int i = 0 ; i < SELECTOR_NB_INDICES ; i++) {
        unsigned int t = s.strassen_threshold[i];
        if (t) printf("%u %u\n", i, t);
    }
}

int main(int argc, char * argv[])
{
    unsigned int m = 0, n = 0, d = 0;
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
        if (m == 0) {
            m = atoi(argv[0]);
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
    // if (d == 0 || n == 0) usage();


    tune_strassen<cantor_fft>("cantor_fft");
    tune_strassen<fake_fft>("fake_fft");

#if 0
    /* Typical situation
     *
     * E * pi_left:
     *
     * n x (m + n), degree d times (m + n) x (m + n), degree d/4
     *
     * pi_left * pi_right
     *
     * (m + n) x (m + n), degree d/4 times the same.
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
#endif
}


