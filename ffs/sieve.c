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
#include "timing.h"

// #define EXPENSIVE_CHECK

// usage: ./a.out q rho
int main(int argc, char **argv)
{
    ffspol_t ffspol[2]; 
    qlat_t qlat;
    factorbase_t FB[2];
    int I, J;  // strict bound on the degrees of the (i,j)
    unsigned int II, JJ; // corresponding maximum integer
    unsigned char threshold[2] = { 50, 50};  // should not be fixed here.
    int lpb[2] = { 25, 25};  // should not be fixed here.
    I = 9; J = 9;
    int noerr;
    
    II = 1; JJ = 1;
    for (int i = 0; i < I; ++i)
        II *= FP_SIZE;
    for (int i = 0; i < J; ++i)
        JJ *= FP_SIZE;

    ffspol_init(ffspol[0]);
    ffspol_init(ffspol[1]);

#ifdef USE_F2
    // Hardcoded GF(2^127) example.
    /*
       F = x^5 + x^4 + x^3 + x^2 + x + t^2
       G =  x*t^25 + x*t^24 + x*t^23 + x*t^21 + x*t^17 + x*t^16 + x*t^15
       + x*t^13 + x*t^10 + x*t^9 + x*t^7 + x*t^5 + x*t^3 + x + t^26 +
       t^25 + t^24 + t^23 + t^21 + t^20 + t^17 + t^14 + t^13 + t^11 + t^9
       + t^6 + t^4 + t^3 + t^2
       Typical special q: 3f29c8d 1e3dc3
    */
    ffspol_set_str(ffspol[0], "4,1,1,1,1,1");
    ffspol_set_str(ffspol[1], "7b26a5c,3a3a6a9");
#else
#ifdef USE_F3
    // Hardcoded GF(3^97) example.
    /*
       F = x^5 + x*t + x + t^2;
       G = x*t^25 + x*t^23 + 2*x*t^22 + 2*x*t^21 + 2*x*t^20 + x*t^18 +
       2*x*t^17 + x*t^16 + 2*x*t^14 + 2*x*t^13 + 2*x*t^11 + x*t^10 +
       2*x*t^9 + 2*x*t^8 + x*t^7 + 2*x*t^6 + x*t^5 + x*t^4 + 2*x*t^2 +
       x*t + 2*x + t^26 + 2*t^25 + t^23 + t^22 + 2*t^20 + t^19 + t^18 +
       2*t^17 + t^16 + t^12 + t^11 + 2*t^10 + t^8 + t^7 + 2*t^6 + t^4 +
       2*t^3 + t + 2;
       Typical special q: 41409965 2aaa4696
     */
    ffspol_set_str(ffspol[0], "10,5,0,0,0,1");
    ffspol_set_str(ffspol[1], "18525901616186,46a19289a6526");
#else
// USE_F3 or USE_F2 must be defined
    ASSERT_ALWAYS(0);
#endif
#endif

    // Read q, rho on command line:
    ASSERT_ALWAYS(argc == 3);
    noerr = sq_set_str(qlat->q, argv[1]);
    ASSERT_ALWAYS(noerr);
    noerr = sq_set_str(qlat->rho, argv[2]);
    ASSERT_ALWAYS(noerr);

    qlat->side = 0; // assume that the special-q is on the algebraic side.

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
    S = (unsigned char *) malloc((II*JJ)*sizeof(unsigned char));
    ASSERT_ALWAYS(S != NULL);
    memset(S, 0, (II*JJ));
    S[0] = 255;
    double t_norms = 0;
    double t_sieve = 0;
    double t_cofact = 0;
    int nrels = 0;
    
    for (int twice = 0; twice < 2; twice++) {
        int side = twice; // put 1-twice here, to do the rational side first.

        // Norm initialization.
        // convention: if a position contains 255, it must stay like
        // this. It means that the other side is hopeless.
        t_norms -= seconds();
        init_norms(S, ffspol[side], I, J, qlat, qlat->side == side);
        t_norms += seconds();

#ifdef EXPENSIVE_CHECK
        // Check special-q divides the norm
        // if (side == qlat->side) {
        if (1) {
            ij_t i, j;
            fppol_t a, b;
            fppol_init(a);
            fppol_init(b);
            fppol_t norm, bigq, rem;
            fppol_init(norm);
            fppol_init(bigq);
            fppol_init(rem);
            fppol_set_sq(bigq, qlat->q);
            for (unsigned int ii = 0; ii < II; ii++) {
                if (!ij_set_ui(i, ii))
                    continue;
                for (unsigned int jj = 0; jj < JJ; jj++) {
                    int position = ii + II*jj;
                    if (position == 0)
                        continue;
                    if (!ij_set_ui(j, jj))
                        continue;
                    if (S[position] != 255) {
                        ij2ab(a, b, i, j, qlat);
                        ffspol_norm(norm, &ffspol[side], a, b);
                        //ASSERT_ALWAYS(fppol_deg(norm) == S[position]);
                        if (side == qlat->side) {
                            fppol_rem(rem, norm, bigq);
                            ASSERT_ALWAYS(fppol_is_zero(rem));
                        }
                    }
                }
            }
            fppol_clear(a);
            fppol_clear(b);
            fppol_clear(norm);
            fppol_clear(bigq);
            fppol_clear(rem);
        }
#endif

        // sieve
        t_sieve -= seconds();
        sieveFB(S, FB[side], I, J, qlat);
        t_sieve += seconds();

        // mark survivors
        // no need to check if this is a valid position
        unsigned char *Sptr = S;
        for (unsigned int j = 0; j < JJ; ++j)
            for (unsigned int i = 0; i < II; ++i, ++Sptr)
            {
                if (*Sptr > threshold[side])
                    *Sptr = 255; 
            }
    }


    t_cofact -= seconds();
    // survivors cofactorization
    {
        unsigned char *Sptr = S;
        fppol_t a, b;
        ij_t ii, jj, gg;
        fppol_init(a);
        fppol_init(b);
        for (unsigned int j = 0; j < JJ; ++j) {
            if (!ij_set_ui(jj, j))
                continue;
            for (unsigned int i = 0; i < II; ++i, ++Sptr)
            {
                if (!ij_set_ui(ii, i))
                    continue;
                if (*Sptr != 255) {
    //                printf("i,j = %u %u\n", i, j);
                    ij_gcd(gg, ii, jj);
                    if (ij_deg(gg) != 0 && i != 1 && j != 1)
                        continue;
                    ij2ab(a, b, ii, jj, qlat);
                    nrels += factor_survivor(a, b, ffspol, lpb);
                }
            }
        fppol_clear(a);
        fppol_clear(b);
    }
    t_cofact += seconds();

    fprintf(stdout, "# Total: %d relations found\n", nrels);
    fprintf(stdout, "# Time spent: %1.1f s (norms); %1.1f s (sieve); %1.1f s (cofact)\n", t_norms, t_sieve, t_cofact);
    fprintf(stdout, "# Rate: %1.4f s/rel\n", (t_norms+t_sieve+t_cofact)/nrels);

    ffspol_clear(ffspol[0]);
    ffspol_clear(ffspol[1]);
    free(S);

    return EXIT_SUCCESS;
}
