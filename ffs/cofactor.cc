#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>

NTL_CLIENT

#include "fppol.h"
#include "types.h"
#include "macros.h"
#include "norm.h"


void fppol2GF2x(GF2X & r, fppol_t p)
{
    unsigned char *pp;
    int d = p->deg;

    int num64 = 1 + d/64;  // exact number of relevant uint64_t.
    int numbytes = num64*8;
    pp = (unsigned char *)malloc(numbytes*sizeof(unsigned char));
    for (int i = 0; i < num64; ++i) {
        pp[8*i  ] = p->limb[i][0];
        pp[8*i+1] = p->limb[i][0] >> 8;
        pp[8*i+2] = p->limb[i][0] >> 16;
        pp[8*i+3] = p->limb[i][0] >> 24;
        pp[8*i+4] = p->limb[i][0] >> 32;
        pp[8*i+5] = p->limb[i][0] >> 40;
        pp[8*i+6] = p->limb[i][0] >> 48;
        pp[8*i+7] = p->limb[i][0] >> 56;
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

// B[2]: smoothness bound on each side corresponding to F[2].
extern "C"
int factor_survivor(fppol_t a, fppol_t b, ffspol_t* F, int *B)
    // , qlat_t qlat)         TODO: will need to get (i,j) as inputs
{
    fppol_t Nab;
    GF2X NN;
    fppol_init(Nab);
    for (int twice = 0; twice < 2; twice++) {
        ffspol_norm(Nab, &F[twice], a, b);
        fppol2GF2x(NN, Nab);
        if (!isBsmooth(NN, B[twice])) {
            fppol_clear(Nab);
            return 0;
        }
    }

    // Yeah... we've got a relation!
    printf("a b:"); 
    vec_pair_GF2X_long factors;
    for (int twice = 0; twice < 2; twice++) {
        ffspol_norm(Nab, &F[twice], a, b);
        fppol2GF2x(NN, Nab);
        CanZass(factors, NN);
        // print factors.
        GF2X::HexOutput = 1;
        cout << factors;
        if (twice == 0)
            printf(":");
    }
    printf("\n");
    fppol_clear(Nab);
    return 1;
}

