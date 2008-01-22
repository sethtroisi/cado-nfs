#include "gf2x.h"
#include <stdio.h>
#include <stdlib.h>

static void mul_basecase(ulong * c, const ulong * a,
			 long na, const ulong * b, long nb)
{
    if (na == nb) {
	switch (na) {
	case 0:		// This can occur in call from KarMul
	    return;
	case 1:
	    mul1(c, a[0], b[0]);
	    return;
#ifndef TUNE_MUL1
	case 2:
	    mul2(c, a, b);
	    return;
	case 3:
	    mul3(c, a, b);
	    return;
	case 4:
	    mul4(c, a, b);
	    return;
	case 5:
	    mul5(c, a, b);
	    return;
	case 6:
	    mul6(c, a, b);
	    return;
	case 7:
	    mul7(c, a, b);
	    return;
	case 8:
	    mul8(c, a, b);
	    return;
	case 9:
	    mul9(c, a, b);
	    return;
#endif
	default:
	    fprintf(stderr, "basecase.c: ran off end of switch\n"
		    "na=nb=%ld ; decrease MUL_KARA_THRESHOLD\n", na);
	    exit(1);
	}
    } else {
	if (na < nb) {
	    long i;
	    c[nb] = mul_1_n(c, b, nb, a[0]);
	    for (i = 1; i < na; i++) {
		c[nb + i] = addmul_1_n(c + i, c + i, b, nb, a[i]);
	    }
	} else {
	    mul_basecase(c, b, nb, a, na);
	}
    }
}
