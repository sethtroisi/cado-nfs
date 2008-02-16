#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>   // for ceiling, floor in cfrac
#include "cado.h"
#include "../utils/mod_ul.h"
#include "fb.h"
#include "../utils/utils.h"
#include "basicnt.h" /* for gcd_ul */

#define LOG_SCALE 1.4426950408889634 /* 1/log(2) to 17 digits, rounded to
                                        nearest. This is enough to uniquely
                                        identify the corresponding IEEE 754
                                        double precision number */

#define LOG_MAX 255.9 /* this should be as near as possible from 256, to
                         get more accuracy in the norm computations, but not
                         too much, otherwise a norm might be rounded to zero */

// General information about the siever
typedef struct {
    // sieving area
    uint32_t I;
    uint32_t J;
    int logI; // such that I = 1<<logI
    // description of the q-lattice
    uint64_t q;
    uint64_t rho;
    int32_t a0, b0, a1, b1;
    // parameters for bucket sieving
    int bucket_size;    // should be around L1 cache size, and a multiple of I.
    int nb_buckets;
    int bucket_limit;   // maximal number of bucket_reports allowed in one bucket.
    unsigned int degree;   /* polynomial degree */
    double scale;     /* LOG_MAX / log (max |F(a0*i+a1*j, b0*i+b1*j)|) */
    mpz_t *fij;       /* coefficients of F(a0*i+a1*j, b0*i+b1*j) */
} sieve_info_t;

/*
 * Some bucket stuff, for later...
 */

/*
 * For the moment, we will keep the bucket reports aligned by adding an
 * 8-bit field that can contain, for instance, the low bits of p.
 *
 * TODO:
 * If the memory pressure becomes too high with this, we can remove this
 * p_low field and pack the reports as follows:
 *    [ x0 ] [ x1 ] [ x2 ] [ x3 ] [logp0] [logp1] [logp2] [logp3] 
 * in the bucket. 
 * 
 * This will be slightly more tricky to store/load bucket reports, but
 * the additional cost should be negligible.
 */

typedef struct {
    uint16_t x;
    char logp;
    uint8_t p_low;  // to keep alignment and help trial div.
} bucket_report_t;

typedef struct {
    bucket_report_t * reports;
    int nbr;
} bucket_t;

/************************** sieve info stuff *********************************/

void
sieve_info_init (sieve_info_t *si, unsigned int d)
{
  unsigned int k;

  si->degree = d;
  si->fij = malloc ((d + 1) * sizeof (mpz_t));
  for (k = 0; k <= d; k++)
    mpz_init (si->fij[k]);
}

void
sieve_info_clear (sieve_info_t *si)
{
  unsigned int d = si->degree;
  unsigned int k;

  for (k = 0; k <= d; k++)
    mpz_clear (si->fij[k]);
  free (si->fij);
}

/*****************************************************************************/

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
special_case_0 (unsigned long p)
{
#ifndef NO_WARNING
  fprintf(stderr, "Warning: special_case_0(%lu) not implemented\n", p);
#endif
}

static void
special_case_p (unsigned long p)
{
#ifndef NO_WARNING
  fprintf(stderr, "Warning: special_case_p(%lu) not implemented\n", p);
#endif
}

#if 0
static void
sieve_slow (unsigned char *S, const factorbase_degn_t *fb,
        const sieve_info_t *si)
{
    fbprime_t p;
    fbprime_t r, R;
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
                special_case_0 (p);
                continue;
            } 
            if (r == p) {
                special_case_p (p);
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
#endif


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
    int64_t a0, a1, b0, b1, I, J;
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
        do {
            a0 += b0;
            a1 += b1;
        } while (a0 + b0 <= 0);
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
        do {
            b0 += a0;
            b1 += a1;
        } while (b0 + a0 >= 0);
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
                special_case_0 (p);
                continue;
            } 
            if (r == p) {
                special_case_p (p);
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
    fprintf (stderr, "small primes sieved in %f sec\n", seconds()-tm);
}

#ifdef TRY_SSE
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
#endif

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
                special_case_0 (fb->p);
                how_many--;
            } else if (r[how_many-1] == p[how_many-1]) {
                special_case_p (fb->p);
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
                special_case_0 (p);
                continue;
            } 
            if (r == p) {
                special_case_p (p);
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
    fprintf (stderr, "large primes sieved in %f sec\n", seconds()-tm);
}

// Conversions between different representations for sieve locations:
//   x          is the index in the sieving array. x in [0,I*J[
//   (i,j)      is the coordinates in the q-lattice. i in [-I/2,I/2[
//                                                   j in [0,J[
//   (a,b)      is the original coordinates. a is signed, b is unsigned.

void 
xToIJ(int *i, unsigned int *j, int x, sieve_info_t * si)
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

/* Warning: b might be negative, in which case we return (-a,-b) */
void
xToAB(int64_t *a, uint64_t *b, int x, sieve_info_t * si)
{
    int i, j;
    int64_t c;

    i = (x % (si->I)) - (si->I >> 1);
    j = x / si->I;
    *a = (int64_t) i * (int64_t) si->a0 + (int64_t) j * (int64_t) si->a1;
    c =  (int64_t) i * (int64_t) si->b0 + (int64_t) j * (int64_t) si->b1;
    if (c >= 0)
      *b = c;
    else
      {
        *a = -*a;
        *b = -c;
      }
}

/*********************** norm computation ************************************/

/* Puts in fij[] the coefficients of f'(i) = F(a0*i+a1, b0*i+b1).
   Assumes the coefficients of fij[] are initialized. */
static void
fij_from_f (mpz_t *fij, mpz_t *f, int d, int32_t a0, int32_t a1,
            int32_t b0, int32_t b1)
{
  int k, l;
  mpz_t *g; /* will contain the coefficients of (b0*i+b1)^l */
  mpz_t f0;

  for (k = 0; k <= d; k++)
    mpz_set (fij[k], f[k]);

  g = malloc ((d + 1) * sizeof (mpz_t));
  for (k = 0; k <= d; k++)
    mpz_init (g[k]);
  mpz_init (f0);

  /* Let h(x) = quo(f(x), x), then F(x,y) = H(x,y)*x + f0*y^d, thus
     F(a0*i+a1, b0*i+b1) = H(a0*i+a1, b0*i+b1)*(a0*i+a1) + f0*(b0*i+b1)^d.
     We use that formula recursively. */

  mpz_set_ui (g[0], 1); /* g = 1 */

  for (k = d - 1; k >= 0; k--)
    {
      /* invariant: we have already translated coefficients of degree > k,
         in f[k+1..d], and g = (b0*i+b1)^(d - (k+1)), with coefficients in
         g[0..d - (k+1)]:
         f[k]   <- f[k] + a1*f[k+1]
         ...
         f[l] <- a0*f[l]+a1*f[l+1] for k < l < d
         ...
         f[d] <- a0*f[d] */
      mpz_swap (f0, fij[k]); /* save the new constant coefficient */
      mpz_mul_si (fij[k], fij[k + 1], a1);
      for (l = k + 1; l < d; l++)
        {
          mpz_mul_si (fij[l], fij[l], a0);
          mpz_addmul_si (fij[l], fij[l + 1], a1);
        }
      mpz_mul_si (fij[d], fij[d], a0);

      /* now compute (b0*i+b1)^(d-k) from the previous (b0*i+b1)^(d-k-1):
         g[d-k] = b0*g[d-k-1]
         ...
         g[l] = b1*g[l]+b0*g[l-1] for 0 < l < d-k
         ...
         g[0] = b1*g[0]
      */
      mpz_mul_si (g[d - k], g[d - k - 1], b0);
      for (l = d - k - 1; l > 0; l--)
        {
          mpz_mul_si (g[l], g[l], b1);
          mpz_addmul_si (g[l], g[l-1], b0);
        }
      mpz_mul_si (g[0], g[0], b1);

      /* now g has degree d-k, and we add f0*g */
      for (l = k; l <= d; l++)
        mpz_addmul (fij[l], g[l - k], f0);
    }

  mpz_clear (f0);
  for (k = 0; k <= d; k++)
    mpz_clear (g[k]);
  free (g);
}

#if 0
/* f(x) -> f(x+1) */
static void
translate_right (mpz_t *f, int d)
{
  int k, l;
  
  for (k = d - 1; k >= 0; k--)
    {
      /* if f(x) = h(x)*x + f0, then f(x+1) = h(x+1)*(x+1) + f0:
         f[d] is unchanged
         ...
         f[l] <- f[l] + f[l+1] */
      for (l = k; l < d; l++)
        mpz_add (f[l], f[l], f[l + 1]);
    }
}
#endif

/* f(x) -> f(x-1) */
static void
translate_left (mpz_t *f, int d)
{
  int k, l;
  
  for (k = d - 1; k >= 0; k--)
    {
      /* if f(x) = h(x)*x + f0, then f(x-1) = h(x-1)*(x-1) + f0:
         f[d] is unchanged
         ...
         f[l] <- f[l] - f[l+1] */
      for (l = k; l < d; l++)
        mpz_sub (f[l], f[l], f[l + 1]);
    }
}

/* f(x) -> f(x+2^i) */
static void
translate_right_2exp (mpz_t *f, unsigned int d, unsigned long i)
{
  unsigned int k;

  for (k = 1; k <= d; k++)
    mpz_mul_2exp (f[k], f[k], i * k);
}

/* f(x) -> f(x-2^i) */
static void
translate_left_2exp (mpz_t *f, unsigned int d, unsigned long i)
{
  unsigned int k;

  for (k = 1; k <= d; k++)
    {
      mpz_mul_2exp (f[k], f[k], i * k);
      if (k & 1)
        mpz_neg (f[k], f[k]);
    }
}

/* f(x) -> f(x+a) */
static void
translate_hori_si (mpz_t *f, int d, long a)
{
  int k, l;
  
  for (k = d - 1; k >= 0; k--)
    {
      for (l = k; l < d; l++)
        mpz_addmul_si (f[l], f[l + 1], a);
    }
}

/* f(x,y) -> f(x,y+1) */
static void
translate_up (mpz_t *f, int d)
{
  int k, l;
  
  for (k = 1; k <= d; k++)
    {
      /* if f(y) = h(y)*y + c, then f(y+1) = h(y+1)*(y+1) + c:
         f[0] is unchanged
         ...
         f[l] <- f[l] + f[l-1] */
      for (l = k; l > 0; l--)
        mpz_add (f[l], f[l], f[l - 1]);
    }
}

/* f(x,y) -> f(x,y+a), a positive: up, a negative: down */
static void
translate_vert_si (mpz_t *f, int d, long a)
{
  int k, l;
  
  for (k = 1; k <= d; k++)
    {
      for (l = k; l > 0; l--)
        mpz_addmul_si (f[l], f[l - 1], a);
    }
}

/* returns the maximal value of log |F(a,b)| for 
   a = a0 * i + a1 * j, b = b0 * i + b1 * j,
   -I/2 <= i <= I/2, 0 <= j <= J */
static double
get_maxnorm (cado_poly cpoly, sieve_info_t si, mpz_t *fij)
{
  mpz_t max_norm, qq;
  double ret;
  unsigned int d, k, j0;
  int i0;
  uint64_t time = cputime ();

  d = cpoly->degree;

  mpz_init_set_ui (max_norm, 0);
  mpz_init (qq);

  fij_from_f (fij, cpoly->f, d, si.a0, si.a1, si.b0, si.b1);

  /* do not divide by q, since at the end we must take q into account in the
     remaining norm */
#if 0
  /* divide the coefficients of fij by q */
  mpz_set_uint64 (qq, si.q);
  for (k = 0; k <= d; k++)
    {
      ASSERT_ALWAYS (mpz_divisible_p (fij[k], qq));
      mpz_divexact (fij[k], fij[k], qq);
    }
#endif

#ifdef VERBOSE
  fprintf (stderr, "F_q(x,1) = ");
  fprint_polynomial (stderr, fij, d);
#endif

  /* for 0 <= i <= I/2, we have F_q(i,0) = fij[d]*i^d, thus the maximum is
     attained at i=I/2, which is on the right border treated below */

  /* translate to i=I/2 */
  translate_right_2exp (fij, d, si.logI - 1);

  /* F_q(I/2,0) is the leading coefficient */
  mpz_abs (max_norm, fij[d]);
  i0 = si.I / 2;
  j0 = 0;

  /* now consider the right border i = I/2, 0 < j <= J */
  for (k = 1; k <= si.J; k++)
    {
      translate_up (fij, d);
      /* now the f[k] are the coefficients of j -> F(I/2, k + j),
         with f[d] the constant coefficient */
      if (mpz_cmpabs (fij[d], max_norm) > 0)
        {
          mpz_abs (max_norm, fij[d]);
          i0 = si.I / 2;
          j0 = k;
        }
    }

  translate_vert_si (fij, d, - (long) si.J); /* back in (I/2, 0) */

  /* back in (0,0) */
  for (k = 1; k <= d; k++)
    mpz_div_2exp (fij[k], fij[k], (si.logI - 1) * k);

  /* translate in (0,J) */
  for (k = 0; k < d; k++)
    {
      /* multiply fij[k] by J^(d-k) */
      mpz_ui_pow_ui (qq, si.J, d - k);
      mpz_mul (fij[k], fij[k], qq);
    }

  /* translate in (I/2, J) */
  translate_hori_si (fij, d, si.I / 2);

  /* upper border: -I/2 <= i < I/2, j = J */
  for (k = 0; k < si.I; k++)
    {
      translate_left (fij, d);
      if (mpz_cmpabs (fij[0], max_norm) > 0)
        {
          mpz_abs (max_norm, fij[0]);
          i0 = (int) si.I / 2 - (k + 1);
          j0 = si.J;
        }
    }
  
  /* now in (-I/2, J), go back in (0,J) */
  translate_hori_si (fij, d, si.I / 2);

  /* now go back in (0,0) */
  for (k = 0; k < d; k++)
    {
      /* divide fij[k] by J^(d-k) */
      mpz_ui_pow_ui (qq, si.J, d - k);
      mpz_divexact (fij[k], fij[k], qq);
    }

  /* translate to i=-I/2: the coefficients of F(-I/2, j) are f[k]*(-I/2)^k */
  translate_left_2exp (fij, d, si.logI - 1);

  /* now consider the left border i = -I/2, 0 <= j <= J */
  for (k = 0; k < si.J; k++)
    {
      /* now the f[k] are the coefficients of j -> F(-I/2, k + j),
         with f[d] the constant coefficient */
      if (mpz_cmpabs (fij[d], max_norm) > 0)
        {
          mpz_abs (max_norm, fij[d]);
          i0 = - (int) si.I / 2;
          j0 = k;
        }
      translate_up (fij, d);
    }
  ret = log (mpz_get_d (max_norm));
  fprintf (stderr, "max norm F_q(%d,%u) has %1.2f bits [time=%lums]\n",
           i0, j0, ret * LOG_SCALE, cputime () - time);

  translate_vert_si (fij, d, - (int) si.J); /* back to (-I/2, 0) */
  for (k = 1; k <= d; k++) /* back to (0,0) */
    {
      mpz_div_2exp (fij[k], fij[k], (si.logI - 1) * k);
      if (k & 1)
        mpz_neg (fij[k], fij[k]);
    }

  mpz_clear (max_norm);
  mpz_clear (qq);
  return ret;
}

/* v <- f(i,j), where f is of degree d */
static void
eval_fij (mpz_t v, mpz_t *f, unsigned int d, long i, unsigned long j)
{
  unsigned int k;
  mpz_t jpow;

  mpz_init_set_ui (jpow, 1);
  mpz_set (v, f[d]);
  for (k = d; k-- > 0;)
    {
      mpz_mul_si (v, v, i);
      mpz_mul_ui (jpow, jpow, j);
      mpz_addmul (v, f[k], jpow);
    }
  mpz_clear (jpow);
}

/* initialize S[i+I/2+j*I] to round(log(|F_q(i,j)|)/log(rho))
   where log(maxnorm)/log(rho) = 255, i.e., rho = maxnorm^(1/255) */
static void
init_norms (unsigned char *S, cado_poly cpoly, sieve_info_t si)
{
  unsigned int i, j, k, l;
  unsigned int d = cpoly->degree;
  mpz_t *t;
  uint64_t time = cputime ();

  /* si.scale = LOG_MAX / log(max |F(a0*i+a1*j, b0*i+b1*j)|) */
  
  /* invariant: si.fij is the f(i,j) polynomial at (0, 0) */

  /* on each row (j fixed), the norms are a polynomial in i of degree d,
     which we evaluate using a table-of-difference method */

  t = malloc ((d + 1) * sizeof (mpz_t));
  for (k = 0; k <= d; k++)
    mpz_init (t[k]);

  for (i = j = 0; j < si.J; j++)
    {
      for (k = 0; k <= d; k++)
        eval_fij (t[k], si.fij, d, - ((long) si.I) / 2 + (long) k, j);
      
      /* initialize table of differences */
      for (k = 1; k <= d; k++)
        for (l = d; l >= k; l--)
          mpz_sub (t[l], t[l], t[l - 1]);

      /* now compute norms */
      for (k = 0; k < si.I; k++)
        {
          /* F_q(-I/2+k, j) is t[0] */
          S[i++] = (unsigned char) (log (fabs (mpz_get_d (t[0]))) * si.scale);

          /* update table of differences */
          for (l = 0; l < d; l++)
            mpz_add (t[l], t[l], t[l + 1]);
        }
    }

  for (k = 0; k <= d; k++)
    mpz_clear (t[k]);
  free (t);

  fprintf (stderr, "Computing norms took %lums\n", cputime () - time);
}

/* check that the double x fits into an int32_t */
#define fits_int32_t(x) \
  ((double) INT32_MIN <= (x)) && ((x) <= (double) INT32_MAX)

void
SkewGauss (sieve_info_t *si, double skewness)
{
  double a[2], b[2], q;

  a[0] = (double) si->q;
  ASSERT_ALWAYS(a[0] < 9007199254740992.0); /* si.q should be less than 2^53
                                               so that a[0] is exact */
  a[1] = 0.0;
  b[0] = (double) si->rho;
  skewness = rint (skewness);
  b[1] = skewness;
  while (1)
    {
      q = (a[0] * b[0] + a[1] * b[1]) / (b[0] * b[0] + b[1] * b[1]);
      q = rint (q);
      if (q == 0.0)
        {
          /* si->a0 and si->a1 should the larger values, in
             a[0] and b[0] */
          ASSERT_ALWAYS(fits_int32_t(a[0]));
          si->a0 = (int32_t) a[0];
          a[1] /= skewness;
          ASSERT_ALWAYS(fits_int32_t(a[1]));
          si->b0 = (int32_t) a[1];
          ASSERT_ALWAYS(fits_int32_t(b[0]));
          si->a1 = (int32_t) b[0];
          b[1] /= skewness;
          ASSERT_ALWAYS(fits_int32_t(b[1]));
          si->b1 = (int32_t) b[1];
          return;
        }
      a[0] -= q * b[0];
      a[1] -= q * b[1];
      q = (a[0] * b[0] + a[1] * b[1]) / (a[0] * a[0] + a[1] * a[1]);
      q = rint (q);
      if (q == 0.0)
        {
          /* exchange a and b */
          ASSERT_ALWAYS(fits_int32_t(a[0]));
          si->a1 = (int32_t) a[0];
          a[1] /= skewness;
          ASSERT_ALWAYS(fits_int32_t(a[1]));
          si->b1 = (int32_t) a[1];
          ASSERT_ALWAYS(fits_int32_t(b[0]));
          si->a0 = (int32_t) b[0];
          b[1] /= skewness;
          ASSERT_ALWAYS(fits_int32_t(b[1]));
          si->b0 = (int32_t) b[1];
          return;
        }
      b[0] -= q * a[0];
      b[1] -= q * a[1];
    }
}

/********************* factoring on rational side ****************************/

/* Those routines are a quick-and-dirty implementation which factors remaining
   norms on the rational side, after sieving on the algebraic side. It will
   probably become obsolete once sieving is done on the rational side.
*/

/* return ceil(log(k)/log(2)), assumes k >= 1 */
static unsigned int
ceil_log2 (unsigned long k)
{
  unsigned long l = 1;

  while (k > 1)
    {
      l ++;
      k = (k + 1) / 2;
    }
  return l;
}

/* Compute gcd(P,t) for all elements of the tree T.
   T[h][0] is the root of the tree. */
void
GcdTree (mpz_t **T, unsigned int *lT, unsigned long h, mpz_t P)
{
  mpz_t *G;
  unsigned long l, m;
#ifdef VERBOSE
  double tm = cputime ();
#endif
  
  G = (mpz_t*) malloc (lT[0] * sizeof (mpz_t));
  for (l = 0; l < lT[0]; l++)
    mpz_init (G[l]);

  mpz_gcd (G[0], T[h][0], P); /* gcd of the root */
#ifdef VERBOSE
  fprintf (stderr, "Tree root has %lu bits, prime product has %lu bits\n",
           mpz_sizeinbase (T[h][0], 2), mpz_sizeinbase (P, 2));
  fprintf (stderr, "Root gcd has size %lu bits\n", mpz_sizeinbase (G[0], 2));
  fprintf (stderr, "Root gcd took %1.0fms\n", cputime () - tm);
#endif
  for (l = h; l-- > 0;)
    { /* level l, we have lT[l] elements, the gcds of the upper level are
         in G[0...lT[l+1]-1] = G[0..ceil(lT[l]/2)-1] */
      for (m = lT[l]; m-- > 0;)
        mpz_gcd (G[m], T[l][m], G[m/2]);
    }
#ifdef VERBOSE
  fprintf (stderr, "Gcd tree took %1.0fms\n", cputime () - tm);
#endif

  /* now G[l] should divide T[0][l] for 0 <= l < lT[0] */
  for (l = 0; l < lT[0]; l++)
    {
      // ASSERT_ALWAYS (mpz_divisible_p (T[0][l], G[l]));
      do
        {
          mpz_divexact (T[0][l], T[0][l], G[l]);
          /* primes in G[l] might divide several times T[0][l] */
          mpz_gcd (G[l], T[0][l], G[l]);
        }
      while (mpz_cmp_ui (G[l], 1) != 0);
    }

  /* FIXME: now that the T[0][l] are smaller, should we recompute the
     product tree? Maybe from time to time. */

  for (l = 0; l < lT[0]; l++)
    mpz_clear (G[l]);
  free (G);
}

/* divide norm by all primes <= B, and print them in buf */
void
trial_divide (char *buf, mpz_t norm, const unsigned long B)
{
  unsigned long p, n = 0, P = B;
  double d;

  mpz_abs (norm, norm);
  d = sqrt (mpz_get_d (norm));
  if (d < (double) P)
    P = (unsigned long) d;
  for (p = 2; p <= P; p = getprime (p))
    {
      while (mpz_divisible_ui_p (norm, p))
        {
          mpz_divexact_ui (norm, norm, p);
          buf += sprintf (buf, "%lx,", p);
          n ++;
          d = sqrt (mpz_get_d (norm));
          if (d < (double) P)
            P = (unsigned long) d;
        }
      /* if norm <= p^2 then norm can only have one remaining prime,
         thus it suffices to check p <= sqrt(norm) */
    }
  if (mpz_cmp_ui (norm, B) <= 0)
    {
      buf += sprintf (buf, "%lx,", mpz_get_ui (norm));
      mpz_set_ui (norm, 1);
      n ++;
    }
  if (n == 0)
    buf[0] = '\0';
  else
    buf[-1] = '\0'; /* replace last ',' by end of string */
  getprime (0);
}

/* build a product tree with the rational norms,
   of height h = ceil(log2(reports)) */
static void
build_norms_rational (cado_poly cpoly, sieve_info_t *si, int *report_list,
                      unsigned long reports)
{
  mpz_t ci, cj, **T, norm;
  int i;
  unsigned int k, j, l;
  unsigned int *lT; /* lT[j] = length of T[j] */
  unsigned int h = ceil_log2 (reports);
  unsigned long p;         /* current prime */
  mpz_t *P = NULL;         /* prime accumulator */
  unsigned long aP = 0;    /* allocated cells in P[] */
  unsigned long lP;        /* number of primes in P[] */
#ifdef VERBOSE
  double tm = cputime ();
#endif
  unsigned long final_reports = 0;
  int64_t a;
  uint64_t b;
  char bufr[256], bufa[256];

  /* on the rational side, we still are in the q-lattice
     a = a0*i+a1*j, b = b0*i+b1*j, thus we have to consider
     the values of g[1]*(a0*i+a1*j) + g[0]*(b0*i+b1*j),
     which are (g[1]*a0+g[0]*b0)*i + (g[1]*a1+g[0]*b1)*j */

  mpz_init (ci); /* ci = g[1]*a0+g[0]*b0 */
  mpz_init (cj); /* cj = g[1]*a1+g[0]*b1 */
  mpz_init (norm);

  mpz_mul_si (ci, cpoly->g[1], si->a0);
  mpz_addmul_si (ci, cpoly->g[0], si->b0);
  mpz_mul_si (cj, cpoly->g[1], si->a1);
  mpz_addmul_si (cj, cpoly->g[0], si->b1);
  // gmp_fprintf (stderr, "ci=%Zd cj=%Zd\n", ci, cj);

  /* builds the product tree */
  T = (mpz_t**) malloc ((h + 1) * sizeof (mpz_t*));
  lT = (unsigned int*) malloc ((h + 1) * sizeof (unsigned int));
  for (l = 0; l <= h; l++)
    {
      /* T[0] is the level of the leaves, and T[h] is the root of the tree */
      lT[l] = 1 + ((reports - 1) >> l);
      T[l] = (mpz_t*) malloc (lT[l] * sizeof (mpz_t));
      for (k = 0; k < lT[l]; k++)
        {
          mpz_init (T[l][k]);
          if (l == 0) /* compute rational norms */
            {
              xToIJ (&i, &j, report_list[k], si);
              mpz_mul_si (T[l][k], ci, i);
              mpz_addmul_ui (T[l][k], cj, j);
            }
          else /* multiply two subtrees */
            {
              if (2 * k + 1 < lT[l-1])
                mpz_mul (T[l][k], T[l-1][2*k], T[l-1][2*k+1]);
              else
                mpz_set (T[l][k], T[l-1][2*k]);
            }
        }
    }
#ifdef VERBOSE
  size_t st;               /* size in bits of root of product tree */
  fprintf (stderr, "Building the product tree took %1.0fms\n",
           cputime () - tm);
  st = mpz_sizeinbase (T[h][0], 2);
  fprintf (stderr, "Root of product tree has size %lu bits\n", st);
#endif

  for (p = 2, lP = 0; p <= cpoly->rlim; p = getprime (p))
    {
      /* invariant: lP is the number of accumulated primes.
         At first step we have P=[2] then P=[1, 2*3],
         then [5, 2*3], then [1, 1, 2*3*5*7] */
      lP ++;
      if ((lP & (lP - 1)) == 0) /* lP = 2^m: need one more cell in P[] */
        {
          aP ++;
          P = (mpz_t*) realloc (P, aP * sizeof (mpz_t));
          mpz_init_set_ui (P[aP - 1], 1);
        }
      /* FIXME: accumulate t primes instead of one at a time in P[0] */
      if (lP & 1)
        mpz_set_ui (P[0], p);
      else /* lP even */
        {
          mpz_mul_ui (P[0], P[0], p);
          for (k = 0; ((lP >> k) & 1) == 0; k++)
            {
              mpz_mul (P[k+1], P[k+1], P[k]);
              mpz_set_ui (P[k], 1);
            }
          if (lP == (1UL << k))
            {
              if (mpz_sizeinbase (P[k], 2) > 65536)
                {
                  GcdTree (T, lT, h, P[k]);
                  mpz_set_ui (P[k], 1);
                  lP = 0;
                }
            }
        }
    }
  getprime (0);
  if (lP != 0) /* accumulate the remaining primes and call GcdTree() */
    {
      for (l = 1; l < aP; l++)
        mpz_mul (P[l], P[l], P[l-1]);
      GcdTree (T, lT, h, P[aP-1]);
    }

  /* check leftover norms that are smaller than the given bound */
  for (l = 0; l < reports; l++)
    {
      int large_size = mpz_sizeinbase (T[0][l], 2);
      /* if T[0][l] is a prime > 2^lpbr, we can ignore the relation */
      if (large_size <= cpoly->mfbr &&
          !(large_size > cpoly->lpbr && mpz_probab_prime_p (T[0][l], 1)))
        {
          /* report_list[l] is I/2+i+j*I */
          xToAB (&a, &b, report_list[l], si);
          eval_fij (norm, cpoly->g, 1, a, b);
          /* we know the product of large primes is T[0][l], thus we can
             divide them out */
          mpz_divexact (norm, norm, T[0][l]);
          trial_divide (bufr, norm, cpoly->rlim);
          /* since we divided out large primes, norm should be 1 here */
          ASSERT_ALWAYS(mpz_cmp_ui (norm, 1) == 0);
          eval_fij (norm, cpoly->f, cpoly->degree, a, b);
          /* we know that q divides the norm on the algebraic side */
          mpz_divexact_ui (norm, norm, si->q);
          /* FIXME: on the algebraic side, we could restrict to the
             factor base primes */
          trial_divide (bufa, norm, cpoly->alim);
          large_size = mpz_sizeinbase (norm, 2);
          if (large_size <= cpoly->mfba &&
              !(large_size > cpoly->lpba && mpz_probab_prime_p (norm, 1)))
            {
              final_reports ++;
              if (strlen (bufa) == 0)
                sprintf (bufa, "%lx", si->q);
              else
                sprintf (bufa + strlen (bufa), ",%lx", si->q);
              printf ("%ld,%lu:%s:%s\n", a, b, bufr, bufa);
            }
        }
    }
  fprintf (stderr, "Final reports: %lu\n", final_reports);

  for (l = 0; l < aP; l++)
    mpz_clear (P[l]);
  free (P);
  for (l = 0; l <= h; l++)
    {
      for (k = 0; k < lT[l]; k++)
        mpz_clear (T[l][k]);
      free (T[l]);
    }
  free (T);
  free (lT);
  mpz_clear (ci);
  mpz_clear (cj);
  mpz_clear (norm);
}

/* build a list report_list[] of values of k = i + I/2 + I*j where gcd(i,j)=1
   and norm is small on algebraic side, then pass it to build_norms_rational */
static void
factor_rational (cado_poly cpoly, sieve_info_t *si, unsigned char *S)
{
    /* sieve reports are those for which the remaining bit-size is less than
       lambda * lpb */
    unsigned char report_bound;
    unsigned char min=255;
    int kmin = 0;
    int i;
    unsigned int j, k;
    int *report_list = NULL; /* list of values of k = i + I/2 + I*j
                                          where norm small on algebraic side */
    unsigned long reports = 0;
    int64_t a;
    uint64_t b;

    report_bound = (unsigned char) (cpoly->alambda * cpoly->lpba
                                    * (si->scale / LOG_SCALE));
    // fprintf (stderr, "report_bound=%u\n", report_bound);

    for (k = 0; k < si->I * si->J; ++k)
      {
        if (S[k] < min)
          {
            kmin = k;
            min = S[k];
          }
        if (S[k] <= report_bound)
          {
            xToAB (&a, &b, k, si);
            if (gcd_ul ((unsigned long) labs(a), b) == 1UL)
              {
                reports ++;
                report_list = realloc (report_list, reports * sizeof (int));
                report_list[reports - 1] = k;
              }
          }
      }
    xToIJ(&i, &j, kmin, si);
    fprintf (stderr, "Min of %d is obtained for (i,j) = (%d,%d)\n", min, i, j);
    fprintf (stderr, "Number of reports on algebraic side: %lu\n", reports);

#if 0
    int stat[256];
    for (k = 0; k < 256; ++k)
        stat[k]=0;
    for (k = 0; k < si->I * si->J; ++k)
        stat[S[k]]++;
    for (k = 0; k < 256; ++k)
      fprintf (stderr, "[%d]%d ", k, stat[k]);
    fprintf (stderr, "\n");
#endif

    build_norms_rational (cpoly, si, report_list, reports);

    free (report_list);
}

/*************************** main program ************************************/

static void
usage (char *argv0)
{
  fprintf (stderr, "Usage: %s -poly xxx.poly -fb xxx.roots\n", argv0);
  exit (1);
}

// this main is just to play with rsa155
// and run with:
//   ./las rsa155.roots
int main(int argc, char ** argv) {
    char *argv0 = argv[0];
    sieve_info_t si;
    char *fbfilename = NULL, *polyfilename = NULL;
    cado_poly cpoly;
    double tm, t0 = seconds ();

    while (argc > 1 && argv[1][0] == '-')
      {
        if (argc > 2 && strcmp (argv[1], "-fb") == 0)
          {
            fbfilename = argv[2];
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-poly") == 0)
          {
            polyfilename = argv[2];
            argc -= 2;
            argv += 2;
          }
        else
          usage (argv0);
      }

    if (polyfilename == NULL || fbfilename == NULL)
      usage (argv0);

    if (!read_polynomial (cpoly, polyfilename))
      {
        fprintf (stderr, "Error reading polynomial file\n");
        exit (EXIT_FAILURE);
      }

    sieve_info_init (&si, cpoly->degree);

    /* those values do not depend on q */
    si.logI = 13;
    si.I = 1<<si.logI;
    si.J = 1<<(si.logI-1);
    
    /* special q (and root rho) */
    si.q = 16777291;
    si.rho = 4078255;
    SkewGauss (&si, cpoly->skew); /* computes a0, b0, a1, b1 from q, rho, and
                                   the skewness */
    si.scale = LOG_MAX / get_maxnorm (cpoly, si, si.fij);

    unsigned char * S;

    // One more is for the garbage during SSE
    S = (unsigned char *) malloc ((1+si.I*si.J)*sizeof(unsigned char));
    ASSERT_ALWAYS(S != NULL);

    init_norms (S, cpoly, si);

    factorbase_degn_t * fb;
    /* the logarithm scale is LOG_MAX / log(max norm) */
    tm = seconds ();
    fb = fb_read (fbfilename, si.scale, 0);
    ASSERT_ALWAYS(fb != NULL);
    fprintf (stderr, "Reading factor base took %f sec\n", seconds() - tm);

    tm = seconds ();
    //sieve_slow(S, fb, &si);
    sieve_random_access(S, fb, &si);
    fprintf (stderr, "Done sieving in %f sec\n", seconds()-tm);

    factor_rational (cpoly, &si, S);

    fprintf (stderr, "Total sieving time %f sec\n", seconds() - t0);

    sieve_info_clear (&si);
    clear_polynomial (cpoly);

    return 0;
}
