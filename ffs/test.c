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

#ifdef TEST_MAIN

// usage: ./a.out q rho
int main(int argc, char **argv)
{
    ffspol_t ffspol[2]; 
    qlat_t qlat;
    factorbase_t FB[2];
    int I, J;  // strict bound on the degrees of the (i,j)
    unsigned char threshold[2] = { 40, 40};  // should not be fixed here.
    I = 13; J = 12;

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

    // Read a and b on command line
    fppol64_t aa, bb;
    errno = 0;
    aa[0] = strtoul(argv[1], NULL, 10);
    ASSERT_ALWAYS(errno == 0);
    bb[0] = strtoul(argv[2], NULL, 10);
    ASSERT_ALWAYS(errno == 0);

    fppol_t a, b;
    fppol_init(a);
    fppol_init(b);
    fppol_set_64(a, aa);
    fppol_set_64(b, bb);

    int smoothness[2] = {30, 30};
    factor_survivor(a, b, ffspol, smoothness); 

    fppol_clear(a);
    fppol_clear(b);



    return EXIT_SUCCESS;
}

#endif
