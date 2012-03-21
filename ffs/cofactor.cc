#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pXFactoring.h>

NTL_CLIENT

// Yurk, that's ugly.
// Anyway, this piece of C++ is not supposed to stay for long.
#define _Bool int

#include "fppol.h"
#include "types.h"
#include "macros.h"
#include "norm.h"
#include "string.h"

#ifdef USE_F3

void NTL2fppol(fppol_t r, zz_pX & P)
{
    fppol_set_zero(r);
    for (int i = deg(P); i >= 0; --i) {
        zz_p C = coeff(P, i);
        unsigned int c = rep(C);
        fp_t cc;
        cc[0] = c & 1u;
        cc[1] = (c >> 1) & 1u;
        fppol_set_coeff(r, cc, i);
    }
}

void fppol2NTL(zz_pX & R, fppol_t p)
{
    R = 0;
    for (int i = fppol_deg(p); i >= 0; --i) {
        fp_t cc;
        fppol_get_coeff(cc, p, i);
        unsigned int c = cc[0] + 2*cc[1];
        SetCoeff(R, i, c);
    }
}

int isBsmooth(zz_pX & p, int B) 
{
    vec_pair_zz_pX_long V;
    p /= LeadCoeff(p);
    SquareFreeDecomp(V, p);
    for (int i = 0; i < V.length(); ++i) {
        vec_pair_zz_pX_long T;
        zz_pX h = PowerXMod(3, V[i].a);
        NewDDF(T, V[i].a, h);
        for (int j = 0; j < T.length(); ++j) {
            if (T[j].b > B)
                return 0;
        }
    }
    return 1;
}
#define NTLPOL zz_pX
#define NTLPOLlong vec_pair_zz_pX_long

#else

void NTL2fppol(fppol_t r, GF2X & p)
{
    if (p == 0) {
        fppol_set_zero(r);
        return;
    }
    unsigned char *pbytes;
    int degp = deg(p);
    int n = 1 + degp/8;  // number of bytes.
    pbytes = (unsigned char *) malloc (n*sizeof(unsigned char));
    ASSERT_ALWAYS(pbytes != NULL);
    BytesFromGF2X(pbytes, p, n);
    unsigned int alloc = 1 + degp/64;
    if (r[0].alloc < alloc) {
        r[0].limbs = (fppol64_t *) realloc(r[0].limbs, alloc*sizeof(fppol64_t));
        r[0].alloc = alloc;
        ASSERT_ALWAYS(r[0].limbs != NULL);
    }
    memset(r[0].limbs, 0, alloc*sizeof(fppol64_t));
    for (int i = 0; i < n; ++i)
        r[0].limbs[i/8][0] |= ((uint64_t)pbytes[i]) << (8*(i%8));
    r[0].deg = degp;
    free(pbytes);
}

void fppol2NTL(GF2X & r, fppol_t p)
{
    unsigned char *pp;
    int d = p->deg;

    int num64 = 1 + d/64;  // exact number of relevant uint64_t.
    int numbytes = num64*8;
    pp = (unsigned char *)malloc(numbytes*sizeof(unsigned char));
    for (int i = 0; i < num64; ++i) {
        pp[8*i  ] = p->limbs[i][0];
        pp[8*i+1] = p->limbs[i][0] >> 8;
        pp[8*i+2] = p->limbs[i][0] >> 16;
        pp[8*i+3] = p->limbs[i][0] >> 24;
        pp[8*i+4] = p->limbs[i][0] >> 32;
        pp[8*i+5] = p->limbs[i][0] >> 40;
        pp[8*i+6] = p->limbs[i][0] >> 48;
        pp[8*i+7] = p->limbs[i][0] >> 56;
    }
    GF2XFromBytes(r, pp, numbytes);
    trunc(r, r, d+1); // just in case the top limb of p was not "clean".
    free(pp);
}

int isBsmooth(GF2X & p, int B) 
{
    vec_pair_GF2X_long V;
    SquareFreeDecomp(V, p);
    for (int i = 0; i < V.length(); ++i) {
        vec_pair_GF2X_long T;
        DDF(T, V[i].a);
        for (int j = 0; j < T.length(); ++j) {
            if (T[j].b > B)
                return 0;
        }
    }
    return 1;
}

#define NTLPOL GF2X
#define NTLPOLlong vec_pair_GF2X_long
#endif


// B[2]: smoothness bound on each side corresponding to F[2].
extern "C"
int factor_survivor(fppol_t a, fppol_t b, ffspol_t* F, int *B)
    // , qlat_t qlat)         TODO: will need to get (i,j) as inputs
{
    fppol_t Nab;

#ifdef USE_F3
    zz_p::init(3);
#endif

    NTLPOL NN;
    fppol_init(Nab);
    for (int twice = 0; twice < 2; twice++) {
        ffspol_norm(Nab, F[twice], a, b);
        fppol2NTL(NN, Nab);
        if (!isBsmooth(NN, B[twice])) {
            fppol_clear(Nab);
            return 0;
        }
    }

    // Yeah... we've got a relation!
    fppol_out(stdout, a); printf(",");
    fppol_out(stdout, b); printf(":");

    NTLPOLlong factors;
    for (int twice = 0; twice < 2; twice++) {
        ffspol_norm(Nab, F[twice], a, b);
        fppol2NTL(NN, Nab);
        NN /= LeadCoeff(NN);
        CanZass(factors, NN);
        // print factors.
        fppol_t toprint;
        fppol_init(toprint);
        for (int i = 0; i < factors.length(); ++i) {
            NTL2fppol(toprint, factors[i].a);
            for (int j = 0; j < factors[i].b; ++j) {
                fppol_out(stdout, toprint); 
                if ((i < factors.length()-1) || (j < factors[i].b-1))
                    printf(",");
            }
        }
        fppol_clear(toprint);
        if (twice == 0)
            printf(":");
    }
    printf("\n");
    fppol_clear(Nab);
    return 1;
}

