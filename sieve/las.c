#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>   // for ceiling, floor in cfrac
#include "cado.h"
#include "../utils/mod_ul.h"
#include "fb.h"
#include "../utils/utils.h"


// General information about the siever
typedef struct {
    uint32_t I;
    uint32_t J;
    uint64_t q;
    uint64_t rho;
    int32_t a0, b0, a1, b1;
    int logI; // such that I = 1<<logI
    int bucket_size;    // should be around L1 cache size
    int nb_buckets;
    int bucket_limit;   // maximal number of bucket_reports allowed in one bucket.
} sieve_info_t;

/*
 * Some bucket stuff, for later...
 */

typedef struct {
    char logp;
    uint8_t p_low;  // to keep alignment and help trial div.
    uint16_t x;
} bucket_report_t;

typedef struct {
    bucket_report_t * reports;
    int nbr;
} bucket_t;

typedef struct {
    uint32_t r;
    int32_t alpha, beta, gamma, delta;
} prime_info_t;

void
push_bucket_report(bucket_t * B, char logp, uint8_t p_low, uint32_t x)
{
    int nbr = B->nbr;
    B->reports[nbr].logp = logp;
    B->reports[nbr].p_low = p_low;
    B->reports[nbr].x = x;
    B->nbr++;
}

// Compute the root r describing the lattice inside the q-lattice
// corresponding to the factor base prime (p,R).
// Formula: r = - (a1-R*b1)/(a0-Rb0) mod p
// In the case where denominator is zero, returns p.
// Otherwise r in [0,p-1]
fbprime_t
fb_root_in_qlattice(const fbprime_t p, const fbprime_t R, const sieve_info_t * si)
{
    modulus_t m;
    residueul_t RR, aa, bb, x, y;
    mod_initmod_ul(m, p);
    modul_initmod_ul(RR, R);  // already reduced
    // nuemrator
    if (si->a1 < 0) {
        modul_initmod_ul(aa, ((unsigned long) (-si->a1)) % p);
        modul_neg(aa, aa, m);
    } else {
        modul_initmod_ul(aa, ((unsigned long) si->a1) % p);
    }
    if (si->b1 < 0) {
        modul_initmod_ul(bb, ((unsigned long) (-si->b1)) % p);
        modul_neg(bb, bb, m);
    } else {
        modul_initmod_ul(bb, ((unsigned long) si->b1) % p);
    }

    modul_mul(x, bb, RR, m);
    modul_sub(x, aa, x, m);
    if (modul_is0(x, m))
        return 0;
    
    // denominator
    if (si->a0 < 0) {
        modul_initmod_ul(aa, ((unsigned long) (-si->a0)) % p);
        modul_neg(aa, aa, m);
    } else {
        modul_initmod_ul(aa, ((unsigned long) si->a0) % p);
    }
    if (si->b0 < 0) {
        modul_initmod_ul(bb, ((unsigned long) (-si->b0)) % p);
        modul_neg(bb, bb, m);
    } else {
        modul_initmod_ul(bb, ((unsigned long) si->b0) % p);
    }

    modul_mul(y, bb, RR, m);
    modul_sub(y, y, aa, m);
    if (modul_is0(y, m))
        return p;

    // divide
    modul_inv(y, y, m);
    modul_mul(x, x, y, m);
    return modul_get_ul(x, m);
}


static void
special_case_0() 
{
#ifndef NO_WARNING
    fprintf(stderr, "Warning: special_case_0() not implemented\n");
#endif
}

static void
special_case_p() 
{
#ifndef NO_WARNING
    fprintf(stderr, "Warning: special_case_p() not implemented\n");
#endif
}

static void
sieve_slow (unsigned char *S, const factorbase_degn_t *fb,
        const sieve_info_t *si)
{
    fbprime_t p;
    fbprime_t r, R;
    unsigned int d;
    unsigned char logp;
    unsigned char *S_ptr;

    // Loop over all primes in the factor base.
    while (fb->p > 0) {
        unsigned char nr;
        p = fb->p;
        logp = fb->plog;
        unsigned long Is2modp = (unsigned long)(si->I>>1) % p;

        for (nr = 0; nr < fb->nr_roots; ++nr) {
            R = fb->roots[nr];
            r = fb_root_in_qlattice(p, R, si);
            if (r == 0) {
                special_case_0();
                continue;
            } 
            if (r == p) {
                special_case_p();
                continue;
            } 

            unsigned long j;
            for (j = 0; j < si->J; ++j) {
                const uint32_t I = si->I;
                long i0;
                // init i0 = -I/2 + ( (rj+I/2) mod p)
                {
                    modulus_t m;
                    residueul_t x, rr, jj;

                    // TODO:
                    // Assuming that p < 2^32 and we are on 64bit
                    // machine, this can be improved a lot! No need to
                    // use mod_ul lib.
                    mod_initmod_ul(m, p);
                    modul_initmod_ul(rr, r);
                    if (j >= p) 
                        modul_initmod_ul(jj, j % p);
                    else
                        modul_initmod_ul(jj, j);
                    modul_mul(x, rr, jj, m);            // here there is a modulo p
                    modul_add_ul(x, x, Is2modp, m);
                    i0 = (long)modul_get_ul(x, m) - (long)(I>>1) ;
                }

                S_ptr = S + (j*I + (I>>1));
                // sieve one j-line
                for (; i0 < (I>>1); i0 += p)
                    S_ptr[i0] -= logp;
            }
        }
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
}




/*
 * Algorithm by Franke and Kleinjung for lattice sieving of largish
 * primes.
 */

typedef struct {
    int32_t alpha, beta, gamma, delta;  // coordinates of the basis
    uint32_t a, c;       // sieving offsets in x-coordinate
    uint32_t b0, b1;     // thresholds for the branch inside siever
} plattice_info_t;

// Proposition 1 of [FrKl05]:
// Compute a basis <(alpha, beta), (gamma, delta)> of the p-lattice
// inside the q-lattice, such that
//    beta, delta > 0
//    -I < alpha <= 0 <= gamma < I
//    gamma-alpha >= I
//
// Sizes:
//    p is less than 32 bits and I fits easily in 32 bits.
//    So, alpha and beta fit easily in 32 bits, since they are less than I
//    Now, gamma and delta are also bounded by p, so 32 bits is enough
//    However: a and c can be as large as p*I (not both ?).
//    We still store them in 32 bits, since if they are larger, it means
//    that as soon as they are added to the offset for S, the index will
//    be out of range for S and the loop stops. Hence, this is safe to
//    replace a and c by a large value within 32 bits, when they are
//    larger than 32 bits.
//
void
reduce_plattice(plattice_info_t *pli, const fbprime_t p, const fbprime_t r, const sieve_info_t * si)
{
    int64_t a0, a1, b0, b1, c0, c1, I, J, ak;
    int k = 1;
    I = si->I;
    J = si->J;
    a0 = -((int64_t)p); a1 = 0;
    b0 = r;  b1 = 1;
#if 0
    /* subtractive variant of Euclid's algorithm */
    while ( b0 >= I )
      {
        /* a0 < 0, b0 > 0 with |a0| > |b0|: this loop is executed at least
           once */
        do
          {
            a0 += b0;
            a1 += b1;
          }
        while (a0 + b0 <= 0);
        /* b0 > 0, a0 < 0 with |b0| > |a0|: since b0 >= I here, this loop is
           also executed at least once */
        do
          {
            b0 += a0;
            b1 += a1;
          }
        while (b0 + a0 >= 0 && b0 >= I);
    }
    while (a0 <= -I)
      {
        a0 += b0;
        a1 += b1;
      }
    pli->alpha = a0;
    pli->beta = a1;
    pli->gamma = b0;
    pli->delta = b1;
#else 
    /* another subtractive variant of Euclid's algorithm */
    /* Seems to be faster... */
    while ( b0 >= I )
    {
      /* a0 < 0, b0 > 0 with |a0| > |b0| */
        while (a0 + b0 <= 0)
          {
            a0 += b0;
            a1 += b1;
          }
        k++;
        if (-a0 < I)
          {
            /* b0 > 0, a0 < 0, with |b0| > |a0| */
            while (b0 >= I)
              {
                b0 += a0;
                b1 += a1;
              }
            goto case_k_even;
          }
        /* b0 > 0, a0 < 0 with |b0| > |a0| */
        while (b0 + a0 >= 0)
          {
            b0 += a0;
            b1 += a1;
          }
        k++;
    }
    /* k is odd here */
    while (a0 <= -I)
      {
        a0 += b0;
        a1 += b1;
      }
 case_k_even:
    pli->alpha = (int32_t) a0;
    pli->beta = (int32_t) a1;
    pli->gamma = (int32_t) b0;
    pli->delta = (int32_t) b1;
#endif
    assert (pli->beta > 0);
    assert (pli->delta > 0);
    assert ((pli->alpha <= 0) && (pli->alpha > -I));
    assert ((pli->gamma >= 0) && (pli->gamma < I));
    assert (pli->gamma-pli->alpha >= I);

    // WARNING: Here, we assume a lot on a bound on I,J
    // TODO: clean these bound problems
    int64_t aa = ((int64_t)pli->beta)*I + (int64_t)(pli->alpha);
    if (aa > I*J)
        pli->a = (uint32_t)(INT32_MAX/2);
    else
        pli->a = (uint32_t)aa;
    int64_t cc = ((int64_t)pli->delta)*I + (int64_t)(pli->gamma);
    if (cc > I*J)
        pli->c = (uint32_t)(INT32_MAX/2);
    else
        pli->c = (uint32_t)cc;
    pli->b0 = -pli->alpha;
    pli->b1 = I - pli->gamma;
}

static void line_sieve(unsigned char *S, factorbase_degn_t **fb_ptr, 
        const sieve_info_t *si)
{
    const uint32_t I = si->I;
    double tm = seconds();
    unsigned char *S_ptr;
    fbprime_t p, r, R;
    unsigned char logp;
    while ((*fb_ptr)->p > 0 && (*fb_ptr)->p <= I) {
        unsigned char nr;
        p = (*fb_ptr)->p;
        logp = (*fb_ptr)->plog;
        unsigned long Is2modp = (unsigned long)(I>>1) % p;

        for (nr = 0; nr < (*fb_ptr)->nr_roots; ++nr) {
            R = (*fb_ptr)->roots[nr];
            r = fb_root_in_qlattice(p, R, si);
            if (r == 0) {
                special_case_0();
                continue;
            } 
            if (r == p) {
                special_case_p();
                continue;
            } 
            
            { 
                // line sieving
                unsigned long j;
                long ii0 = Is2modp; // this will be (rj+I/2) mod p
                for (j = 0; j < si->J; ++j) {
                    long i0;
                    // init i0 = -I/2 + ( (rj+I/2) mod p)
                    i0 = ii0 - (long)(I>>1);
                    S_ptr = S + (j*I + (I>>1));
                    // sieve one j-line
                    for (; i0 < (I>>1); i0 += p)
                        S_ptr[i0] -= logp;
                    // update starting point
                    ii0 += r;
                    if (ii0 >= p)
                        ii0 -= p;

                }
            }
        }
        (*fb_ptr) = fb_next ((*fb_ptr)); // cannot do fb++, due to variable size !
    }
    printf("small primes sieved in %f sec\n", seconds()-tm);
}


typedef uint32_t v4si __attribute__ ((vector_size (16)));

// r = (x < y) ? a : b
static inline void
xmm_cmovlt(v4si *r, v4si x, v4si y, v4si a, v4si b) {
    v4si tmp;
    tmp = __builtin_ia32_pcmpgtd128(y, x); // tmp = x<y ? 0xffff..fff. : 0
    a = __builtin_ia32_pand128(a, tmp);  // a = x<y ? a : 0
    tmp = __builtin_ia32_pandn128(tmp, b); // tmp = x<y ? 0 : b
    *r = __builtin_ia32_por128(a, tmp); 
}

// r = (x >= y) ? a : b
static inline void
xmm_cmovge(v4si *r, v4si x, v4si y, v4si a, v4si b) {
    v4si tmp0, tmp1;
    tmp0 = __builtin_ia32_pcmpgtd128(x, y); // tmp0 = x > y ? 0xffff..fff. : 0
    tmp1 = __builtin_ia32_pcmpeqd128(x, y); // tmp1 = x == y ? 0xffff..fff. : 0
    tmp0 |= tmp1;  // tmp0 = x >= y ? 0xffff..fff. : 0
    a = __builtin_ia32_pand128(a, tmp0);  
    tmp0 = __builtin_ia32_pandn128(tmp0, b); 
    *r = __builtin_ia32_por128(a, tmp0); 
}

// This version of the sieve implements the following:
//   - naive line-sieving for small p
//   - Franke-Kleinjung lattice-sieving for large p
//   - write directly the reports in the sieve-array (no buckets)
// WARNING: STILL AT WORK!
static void
sieve_random_access (unsigned char *S, factorbase_degn_t *fb,
        const sieve_info_t *si)
{
    // Sieve primes <= I.
    line_sieve(S, &fb, si);

    double tm = seconds();
#ifdef TRY_SSE
    if (fb->p == 0)
        return;
    int finished = 0;
    int cur_root = 0;
    fbprime_t p[4];
    fbprime_t r[4];
    fbprime_t R[4];
    unsigned char logp[4];
    int how_many = 0;
    // when we enter this loop, we assume that 
    //   fb->p > 0
    //   cur_root < fb->nb_roots
    // describe the next prime ideal to process.
    // (except if finished is set, of course)
    while (!finished) {
        // Get next 4 primes in fb.
        how_many = 0;
        while (how_many < 4) {
            //printf("Dealing with p = %d, r = %d\n", fb->p, fb->roots[cur_root]);
            p[how_many] = fb->p;
            R[how_many] = fb->roots[cur_root];
            logp[how_many] = fb->plog;
            r[how_many] = fb_root_in_qlattice(p[how_many], R[how_many], si);
            how_many++;
            if (r[how_many-1] == 0) {
                special_case_0();
                how_many--;
            } else if (r[how_many-1] == p[how_many-1]) {
                special_case_p();
                how_many--;
            }
            if (cur_root < fb->nr_roots-1)
                cur_root++;
            else {
                cur_root = 0;
                fb = fb_next (fb);
                if (fb->p == 0) {
                    finished = 1;
                    break;
                }
            }
        }
        if (how_many != 4)
            break;

        // Here we have 4 consecutives prime ideals to sieve.

        const uint32_t I = si->I;
        const uint32_t maskI = I-1;
        plattice_info_t pli[4];
        // precomupte appropriate basis for p-lattices.
        {
            int i;
            for (i = 0; i < 4; ++i) {
                reduce_plattice(pli + i, p[i], r[i], si);
            }
        }

        // Start sieving from (0,0) which is I/2 in x-coordinate
        uint32_t x0, x1, x2, x3;
        x0 = (I>>1);
        x1 = (I>>1);
        x2 = (I>>1);
        x3 = (I>>1);
        const uint32_t IJ = I*si->J;

        v4si vec_b0 = { pli[0].b0, pli[1].b0, pli[2].b0, pli[3].b0 };
        v4si vec_b1 = { pli[0].b1, pli[1].b1, pli[2].b1, pli[3].b1 };
        v4si vec_a = { pli[0].a, pli[1].a, pli[2].a, pli[3].a };
        v4si vec_c = { pli[0].c, pli[1].c, pli[2].c, pli[3].c };
        v4si vec_x = { x0, x1, x2, x3 };
        v4si vec_maskI = { maskI, maskI, maskI, maskI };
        v4si vec_IJ = { IJ, IJ, IJ, IJ };
        v4si vec_zero = { 0, 0, 0, 0};
        asm("## Inner sieving routine starts here!!!\n");
        while (x0 < IJ) {
            uint32_t i;
            v4si vec_i, vec_tmp;
            vec_i = vec_x & vec_maskI;
            S[x0] -= logp[0];
            S[x1] -= logp[1];
            S[x2] -= logp[2];
            S[x3] -= logp[3];
            xmm_cmovge(&vec_tmp, vec_i, vec_b1, vec_a, vec_zero);
            vec_x += vec_tmp;
            xmm_cmovlt(&vec_tmp, vec_i, vec_b0, vec_c, vec_zero);
            vec_x += vec_tmp;
            // x = ( IJ < x ) ? IJ : x
            xmm_cmovlt(&vec_x, vec_IJ, vec_x, vec_IJ, vec_x);
            x0 = ((uint32_t *)(&vec_x))[0];
            x1 = ((uint32_t *)(&vec_x))[1];
            x2 = ((uint32_t *)(&vec_x))[2];
            x3 = ((uint32_t *)(&vec_x))[3];
        }
        asm("## Inner sieving routine stops here!!!\n");
    }
#endif

    // Loop over all primes in the factor base > I
    while (fb->p > 0) {
        unsigned char nr;
        fbprime_t p = fb->p;
        unsigned char logp = fb->plog;

        for (nr = 0; nr < fb->nr_roots; ++nr) {
            fbprime_t r, R;
            R = fb->roots[nr];
            r = fb_root_in_qlattice(p, R, si);
            if (r == 0) {
                special_case_0();
                continue;
            } 
            if (r == p) {
                special_case_p();
                continue;
            } 
            
            const uint32_t I = si->I;
            const uint32_t maskI = I-1;
 
            plattice_info_t pli;
            reduce_plattice(&pli, p, r, si);

            // Start sieving from (0,0) which is I/2 in x-coordinate
            uint32_t x;
            x = (I>>1);
            // TODO: to gain speed, by aligning the start of the
            // loop, in assembly.
            // Besides this small trick, gcc does a good job with
            // this loop: no branching for the if(), and a dozen of
            // instructions inside the while().
            asm("## Inner sieving routine starts here!!!\n");
            while (x < I*si->J) {
                uint32_t i;
                i = x & maskI;   // x mod I
                S[x] -= logp;
                if (i >= pli.b1)
                    x += pli.a;
                if (i < pli.b0)
                    x += pli.c;
            }
            asm("## Inner sieving routine stops here!!!\n");
        }
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
    printf("large primes sieved in %f sec\n", seconds()-tm);
}

// Conversions between different representations for sieve locations:
//   x          is the index in the sieving array. x in [0,I*J[
//   (i,j)      is the coordinates in the q-lattice. i in [-I/2,I/2[
//                                                   j in [0,J[
//   (a,b)      is the original coordinates. a is signed, b is unsigned.

void 
xToIJ(int *i, int *j, int x, sieve_info_t * si)
{
    *i = (x % (si->I)) - (si->I >> 1);
    *j = x / si->I;
}

void
IJTox(int *x, int i, int j, sieve_info_t * si)
{
    *x = i + (si->I)*j + (si->I>>1);
}

void
IJToAB(int64_t *a, uint64_t *b, int i, int j, sieve_info_t * si)
{
    *a = i*si->a0 + j*si->a1;
    *b = i*si->b0 + j*si->b1;
}

void
xToAB(int64_t *a, uint64_t *b, int x, sieve_info_t * si)
{
    int i, j;
    i = (x % (si->I)) - (si->I >> 1);
    j = x / si->I;
    *a = i*si->a0 + j*si->a1;
    *b = i*si->b0 + j*si->b1;
}


// this main is just to play with rsa155
// compile with:
//   gcc -O4 -o las -I.. las.c fb.o ../utils/libutils.a -lgmp -lm
// and run with:
//   ./las rsa155.roots
int main(int argc, char ** argv) {
    sieve_info_t si;
    const double log_scale = 1.4426950408889634073599246810018921374;


    si.I = 1<<14;
    si.logI = 14;
    si.J = 1<<13;
    si.q = 16777291;
    si.rho = 4078255;
    si.a0 = 464271;
    si.b0 = -4;
    si.a1 = -100184;
    si.b1 = 37;

    unsigned char * S;

    // One more is for the garbage during SSE
    S = (unsigned char *) malloc ((1+si.I*si.J)*sizeof(unsigned char));
    ASSERT_ALWAYS(S != NULL);

    factorbase_degn_t * fb;
    fb = fb_read (argv[1], log_scale, 1);
    ASSERT_ALWAYS(fb != NULL);

    // Put 255 everywhere, since we don't want to compute norms for the
    // moment.
    memset(S, 255, si.I*si.J);

    double tm = seconds();
    //sieve_slow(S, fb, &si);
    sieve_random_access(S, fb, &si);
    printf("Done sieving in %f sec\n", seconds()-tm);

    unsigned char min=255;
    int kmin = 0;
    int i, j, k;
    for (k = 0; k < si.I*si.J; ++k)
        if ((S[k] < min) && (S[k]>50) ) {
            kmin = k;
            min = S[k];
        }
    
    xToIJ(&i, &j, kmin, &si);
    printf("Min of %d is obtained for (i,j) = (%d,%d)\n", min, i, j);

    int stat[256];
    for (i = 0; i < 256; ++i)
        stat[i]=0;
    for (i = 0; i < si.I*si.J; ++i)
        stat[S[i]]++;
    for (i = 0; i < 256; ++i)
        printf("[%d]%d ", i, stat[i]);
    printf("\n");

}

