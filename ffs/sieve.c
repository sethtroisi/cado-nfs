#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "macros.h"

#include "types.h"
#include "fppol.h"
#include "ffspol.h"
#include "qlat.h"
#include "fb.h"
#include "latsieve.h"
#include "norm.h"
#include "cofactor.hh"

// usage: ./a.out q rho
int main(int argc, char **argv)
{
    ffspol_t ffspol[2]; 
    qlat_t qlat;
    factorbase_t FB[2];
    int I, J;  // strict bound on the degrees of the (i,j)
    unsigned char threshold[2] = { 40, 40};  // should not be fixed here.
    int lpb[2] = { 22, 22};  // should not be fixed here.
    I = 7; J = 7;
    int noerr;

    // Hardcoded GF(2^127) example.
    
    // F = x^5 + x*t + x + t^2  
    ffspol[0].deg = 5;
    ffspol[0].alloc = 6;
    ffspol[0].coeffs = (fppol_t *)malloc(6*sizeof(fppol_t));
    ASSERT_ALWAYS(ffspol[0].coeffs != NULL);
    for(int i = 0; i < 6; ++i)
        fppol_init(ffspol[0].coeffs[i]);
    fppol_set_ui(ffspol[0].coeffs[0], 4);
    fppol_set_ui(ffspol[0].coeffs[1], 3);
    fppol_set_ui(ffspol[0].coeffs[2], 0);
    fppol_set_ui(ffspol[0].coeffs[3], 0);
    fppol_set_ui(ffspol[0].coeffs[4], 0);
    fppol_set_ui(ffspol[0].coeffs[5], 1);

    // G = x*t^25 + x*t^24 + x*t^22 + x*t^20 + x*t^19 + x*t^15 + x*t^14 +
    // x*t^13 + x*t^12 + x*t^11 + x*t^10 + x*t^9 + x*t^7 + x*t^5 + x*t^2
    // + x*t + x + t^26 + t^21 + t^18 + t^16 + t^15 + t^13 + t^10 + t^9 +
    // t^8 + t^7 + t^6
    ffspol[1].deg = 1;
    ffspol[1].alloc = 2;
    ffspol[1].coeffs = (fppol_t *)malloc(2*sizeof(fppol_t));
    ASSERT_ALWAYS(ffspol[1].coeffs != NULL);
    for(int i = 0; i < 2; ++i)
        fppol_init(ffspol[1].coeffs[i]);
    fppol_set_ui(ffspol[1].coeffs[0], 69576640);
    fppol_set_ui(ffspol[1].coeffs[1], 56164007);

    // Read q, rho on command line:
    ASSERT_ALWAYS(argc == 3);
    noerr = sq_set_str(qlat->q, argv[1]);
    ASSERT_ALWAYS(noerr);
    noerr = sq_set_str(qlat->rho, argv[2]);
    ASSERT_ALWAYS(noerr);

    // Reduce the q-lattice
    noerr = skewGauss(qlat, 0);
    assert (noerr);
    print_qlat_info(qlat);

    // Read the factor bases
    noerr = fbread(&FB[0], "Aroots");
    assert (noerr);
    noerr = fbread(&FB[1], "Rroots");
    assert (noerr);

    // Allocate and init the sieve space
    unsigned char *S;
    S = (unsigned char *) malloc ((1<<I)*(1<<J)*sizeof(unsigned char));
    ASSERT_ALWAYS(S != NULL);
    memset(S, 0, (1<<I)*(1<<J));
    
    for (int twice = 0; twice < 2; twice++) {
        int side = twice; // put 1-twice here, to do the rational side first.

        // Norm initialization.
        // convention: if a position contains 255, it must stay like
        // this. It means that the other side is hopeless.
        init_norms(S, ffspol[side], I, J, qlat);

        // sieve
        sieveFB(S, FB[side], I, J, qlat);

        // mark survivors
        unsigned char *Sptr = S;
        for (unsigned int i = 0; i < (1U<<I); ++i)
            for (unsigned int j = 0; j < (1U<<J); ++j, ++Sptr)
            {
                if (*Sptr > threshold[side])
                    *Sptr = 255; 
            }
    }

    // survivors cofactorization
    {
        unsigned char *Sptr = S;
        fppol_t a, b;
        ij_t ii, jj;
        fppol_init(a);
        fppol_init(b);
        for (unsigned int i = 0; i < (1U<<I); ++i)
            for (unsigned int j = 0; j < (1U<<J); ++j, ++Sptr)
            {
                if (*Sptr != 255) {
                    ii[0] = i;
                    jj[0] = j;
                    ij2ab(a, b, ii, jj, qlat);
                    factor_survivor(a, b, ffspol, lpb);
                }
            }
        fppol_clear(a);
        fppol_clear(b);
    }

    return EXIT_SUCCESS;
}
