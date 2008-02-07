#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
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
    fprintf(stderr, "Warning: special_case_0() not implemented\n");
}

static void
special_case_p() 
{
    fprintf(stderr, "Warning: special_case_p() not implemented\n");
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
//        fprintf(stderr, "done %lu\n", p);
        fb = fb_next (fb); // cannot do fb++, due to variable size !
    }
}

// Conversions between different representations for sieve locations:
//   x          is the index in the sieving array. x in [0,I*J[
//   (i,j)      is the coordinates in the q-lattice. i in [-I/2,I/2[
//                                                   j in [0,J[
//   (a,b)      is the original coordinates. a is signed, b is unsigned.

void 
xToIJ(int *i, int *j, int x, sieve_info_t * si)
{
    *i = (x % (si.I)) - (si.I >> 1);
    *j = x / si.I;
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
    i = (x % (si.I)) - (si.I >> 1);
    j = x / si.I;
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


    si.I = 1<<13;
    si.J = 1<<12;
    si.q = 16777291;
    si.rho = 4078255;
    si.a0 = 464271;
    si.b0 = -4;
    si.a1 = -100184;
    si.b1 = 37;

    unsigned char * S;

    S = (unsigned char *) malloc (si.I*si.J*sizeof(unsigned char));
    ASSERT_ALWAYS(S != NULL);
    // Put 255 everywhere, since we don't want to compute norms for the
    // moment.
    memset(S, 255, si.I*si.J);

    factorbase_degn_t * fb;
    fb = fb_read (argv[1], log_scale, 1);
    ASSERT_ALWAYS(fb != NULL);

    sieve_slow(S, fb, &si);

    unsigned char min=255;
    int kmin = 0;
    int i, j, k;
    for (k = 0; k < si.I*si.J; ++k)
        if ((S[k] < min) && (S[k]>50) ) {
            kmin = k;
            min = S[k];
        }
    
    xToIJ(&i, &j, kmin, si);
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

