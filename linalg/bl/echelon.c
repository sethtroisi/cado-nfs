#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"
#include "mat_ops.h"
#include "echelon.h"
#include <limits.h>

#include <math.h>
#include <string.h>

#define	WBITS	(CHAR_BIT * sizeof(unsigned long))
#define	iceildiv(x,y)	(((x)+(y)-1)/(y))
/////////////////////////////////////////////////////
//
//
//     Forme Echelon
//
//
////////////////////////////////////////////////////




void ConstructPermBit(unsigned long n, unsigned long *US, unsigned long *p)
{
    unsigned long i, j, k, SizeRowA;

    SizeRowA = iceildiv(n, WBITS);


//Alloc
    unsigned long *IN;
    IN = Allocmn(n, n);

    unsigned long *Col;
    Col = Allocmn(1, n);

    unsigned long *p1;
    p1 = Allocmn(n, n);

    memset(IN, 0, n * iceildiv(n, WBITS) * sizeof(unsigned long));
    for (i = 0; i < n; ++i) {
	IN[i * SizeRowA + i / WBITS] = (1UL << i % WBITS);
    }

    unsigned long posin = n - US[0];
    unsigned long posout = 0;
    unsigned long test = 0;

    for (i = 0; i < n; ++i) {
	for (j = 1; j <= US[0]; ++j) {
	    if (US[j] == i) {
		test = 1;
		break;
	    }
	}



	SelectColumn(n, n, IN, i, Col);


	if (test == 1) {
	    for (k = 0; k < SizeRowA; ++k) {
		p1[posin * SizeRowA + k] = Col[k];
	    }
	    posin += 1;
	} else {
	    for (k = 0; k < SizeRowA; ++k) {
		p1[posout * SizeRowA + k] = Col[k];
	    }
	    posout += 1;
	}
	test = 0;
    }


    TransposeBit(n, n, p1, p);


    free(IN);
    free(Col);
    free(p1);


}



void Construction_S_WinvBit(unsigned long n, unsigned long *m,
			    unsigned long *USet, unsigned long *S,
			    unsigned long *Set0, unsigned long *Winv,
			    unsigned long *SGPerm)
{
    unsigned long i, j, posS, k, h, SizeRowA;


    unsigned long *P;

    P = Allocmn(n, n);


    ConstructPermBit(n, USet, P);

    unsigned long *Q;

    Q = Allocmn(n, n);

    unsigned long *T;

    T = Allocmn(n, n);

    unsigned long *Resnn;

    Resnn = Allocmn(n, n);

    TransposeBit(n, n, P, Q);

    DMultBit(n, n, n, Q, m, Resnn);


    DMultBit(n, n, n, Resnn, P, T);


    unsigned long *M;
    M = Allocmn(n, 2 * n);
    memset(M, 0, n * iceildiv(2 * n, WBITS) * sizeof(unsigned long));

    SizeRowA = iceildiv(n, WBITS);

    unsigned long *IN;


    IN = Allocmn(n, n);

    memset(IN, 0, n * iceildiv(n, WBITS) * sizeof(unsigned long));

    for (i = 0; i < n; ++i) {
	IN[i * SizeRowA + i / WBITS] = (1UL << i % WBITS);
    }

    unsigned long SizeRow2A;

    SizeRow2A = iceildiv(n + n, WBITS);

    for (i = 0; i < n; ++i) {
	for (j = 0; j < SizeRow2A; ++j) {

	    if (j < SizeRowA - 1) {

		M[i * SizeRow2A + j] = T[i * SizeRowA + j];
	    }
	    if (j == SizeRowA - 1) {

		if ((n % WBITS > i) && (j - SizeRowA + 1 == i / WBITS))
		    M[i * SizeRow2A + SizeRowA - 1] = (1UL << (i));
		M[i * SizeRow2A + SizeRowA - 1] =
		    (M[i * SizeRow2A + SizeRowA - 1] << n % WBITS);
		M[i * SizeRow2A + SizeRowA - 1] |=
		    T[i * SizeRowA + SizeRowA - 1];
	    }
	    if (j > SizeRowA - 1) {

		if ((i < WBITS - n % WBITS) && (n % WBITS != 0))
		    M[i * SizeRow2A + j] =
			(1UL << (n % WBITS - i % WBITS - 1));
		else {
		    if ((n % WBITS) != 0) {
			M[i * SizeRow2A + j] =
			    (1UL << (i % WBITS - WBITS + n % WBITS));
		    } else if (j - SizeRowA == i / WBITS) {
			M[i * SizeRow2A + j] = (1UL << (i));
		    }


		}

	    }

	}
    }



    posS = 1;

    unsigned long *Set;
    Set = malloc((n + 1) * sizeof(unsigned long));


    for (j = 0; j < n; ++j) {
	{
	    if ((M[j * SizeRow2A + j / WBITS] >> j & 1) == 0) {
		for (k = j; k < n; ++k) {
		    if ((M[k * SizeRow2A + j / WBITS] >> j & 1) != 0) {
			SwapRowsBit(n, n + n, M, j, k);
			break;
		    }
		}
	    }
	}
	if ((M[j * SizeRow2A + j / WBITS] >> j & 1) != 0) {
	    Set[posS] = j;
	    ++posS;

	    MultiplyRowBit(n, n + n, M,
			   (M[j * SizeRow2A + j / WBITS] >> j & 1), j);

	    for (h = 0; h < n; ++h) {
		if (h != j) {
		    AddRowBit(n, n + n, M,
			      (M[h * SizeRow2A + j / WBITS] >> j & 1), j, h);
		}
	    }
	} else {
	    for (k = j; k < n; ++k) {
		if ((M[j * SizeRow2A + (j + n) / WBITS] >> (j + n) & 1) == 0) {
		    if ((M[k * SizeRow2A + (j + n) / WBITS] >> (j + n) & 1) !=
			0) {
			SwapRowsBit(n, n + n, M, j, k);
			break;
		    }
		}
	    }
	    if ((M[j * SizeRow2A + (j + n) / WBITS] >> (j + n) & 1) != 0) {
		for (h = 0; h < n; ++h) {
		    if (h != j) {
			AddRowBit(n, n + n, M,
				  (M[h * SizeRow2A + (j + n) / WBITS] >>
				   (j + n) & 1), j, h);
		    }
		}
		MultiplyRowBit(n, n + n, M, 0, j);
	    }
	}

    }


    unsigned long *Si;
    Si = Allocmn(n, n);
    memset(Si, 0, n * iceildiv(n, WBITS) * sizeof(unsigned long));
    for (i = 0; i < n; ++i) {
	Si[i * SizeRowA + i / WBITS] = (1UL << i % WBITS);
    }

    Set[0] = posS - 1;
    Set0[0] = posS - 1;

//construction of Winv corrected by permutation


    memset(IN, 0, n * iceildiv(n, WBITS) * sizeof(unsigned long));


    unsigned long *C;

    C = Allocmn(n, n);

    for (i = 0; i < n; ++i) {
	SelectColumn(n, n + n, M, n + i, C);
	for (j = 0; j < SizeRowA; ++j)
	    IN[i * SizeRowA + j] = C[j];
    }


    TransposeBit(n, n, IN, T);
    DMultBit(n, n, n, P, T, Resnn);
    TransposeBit(n, n, P, Q);
    DMultBit(n, n, n, Resnn, Q, Winv);


//Desalloc



    free(T);
    free(Q);
    free(P);
    free(Resnn);
    free(M);
    free(IN);

//construction of the permutation given by the initial data


    unsigned long posin = n - USet[0];
    unsigned long posout = 0;

    unsigned long *GPerm;
    GPerm = malloc(n * sizeof(unsigned long));

    unsigned long test = 0;

    for (i = 1; i <= n; ++i) {
	for (j = 1; j <= USet[0]; ++j) {
	    if (USet[j] == i - 1) {
		test = 1;
		break;
	    }
	}
	if (test == 1) {
	    GPerm[posin] = i;
	    posin += 1;
	} else {
	    GPerm[posout] = i;
	    posout += 1;
	}
	test = 0;
    }



//Construction of the image by the previous permutation of the set of columns to use in next step

    for (i = 1; i <= Set[0]; ++i) {
	SGPerm[i] = GPerm[Set[i]] - 1;
    }
    SGPerm[0] = Set[0];


//Construction of Si


    for (i = 0; i < SGPerm[0]; ++i)
	for (j = 0; j < SizeRowA; ++j) {
	    S[i * SizeRowA + j] = Si[SGPerm[i + 1] * SizeRowA + j];
	}

    free(Set);
    free(GPerm);
    free(C);
    free(Si);

}




////////////////////////////////////////////////////////
//
//
//  Functions to older versions of Lanczos
//
///////////////////////////////////////////////////////



/*

Construct a permutation matrix p such that when we multiply p by a the columns of a indexed in the set UseSet (1,,k) appear in the right. 

*/

void ConstructPerm(unsigned long n, unsigned long US[n + 1],
		   unsigned long **p)
{
    unsigned long i, j, k;
    unsigned long IN[n][n];
    for (i = 0; i < n; ++i)
	for (j = 0; j < n; ++j)
	    if (i == j) {
		IN[i][j] = 1;
	    } else {
		IN[i][j] = 0;
	    }
    unsigned long posin = n - US[0];
    unsigned long posout = 0;
    unsigned long test = 0;
    for (i = 0; i < n; ++i) {
	for (j = 1; j <= US[0]; ++j) {
	    if (US[j] == i) {
		test = 1;
		break;
	    }
	}
	if (test == 1) {
	    for (k = 0; k < n; ++k) {
		p[k][posin] = IN[k][i];
	    }
	    posin += 1;
	} else {
	    for (k = 0; k < n; ++k) {
		p[k][posout] = IN[k][i];
	    }
	    posout += 1;
	}
	test = 0;
    }

}












/*

Compute S and Winv

*/

void Construction_S_Winv(unsigned long n, unsigned long **m,
			 unsigned long USet[n + 1], unsigned long **S,
			 unsigned long Set0[1], unsigned long **Winv,
			 unsigned long **SGPerm)
{
    unsigned long i, j, posS, k, h;


    unsigned long **P;

    P = Alloc2(n, n);

    ConstructPerm(n, USet, P);


    unsigned long **Q;

    Q = Alloc2(n, n);

    unsigned long **T;

    T = Alloc2(n, n);

    unsigned long **Resnn;

    Resnn = Alloc2(n, n);

    Transpose(n, n, P, Q);

    DMult(n, n, n, Q, m, Resnn);

    DMult(n, n, n, Resnn, P, T);
    unsigned long M[n][n + n];

    for (i = 0; i < n; ++i) {
	for (j = 0; j < 2 * n; ++j) {
	    if (j < n) {
		M[i][j] = T[i][j];
	    } else {
		if (i == j - n) {
		    M[i][j] = 1;
		} else {
		    M[i][j] = 0;
		}
	    }
	}
    }
    posS = 1;

    unsigned long Set[n + 1];


    for (j = 0; j < n; ++j) {
	{
	    if (M[j][j] == 0) {
		for (k = j; k < n; ++k) {
		    if (M[k][j] != 0) {
			SwapRows(n, n + n, M, j, k);
			break;
		    }
		}
	    }
	}
	if (M[j][j] != 0) {
	    Set[posS] = j;
	    ++posS;
	    MultiplyRow(n, n + n, M, M[j][j], j);

	    for (h = 0; h < n; ++h) {
		if (h != j) {
		    AddRow(n, n + n, M, M[h][j], j, h);
		}
	    }
	} else {
	    for (k = j; k < n; ++k) {
		if (M[j][j + n] == 0) {
		    if (M[k][j + n] != 0) {
			SwapRows(n, n + n, M, j, k);
			break;
		    }
		}
	    }
	    if (M[j][j + n] != 0) {
		for (h = 0; h < n; ++h) {
		    if (h != j) {
			AddRow(n, n + n, M, -M[h][j + n], j, h);
		    }
		}
		MultiplyRow(n, n + n, M, 0, j);
	    }
	}
    }


    unsigned long Si[n][n];
    for (i = 0; i < n; ++i)
	for (j = 0; j < n; ++j) {
	    if (i == j) {
		Si[i][j] = 1;
	    } else {
		Si[i][j] = 0;
	    }
	}
    Set[0] = posS - 1;
    Set0[0] = posS - 1;





//construction of Winv corrected by permutation

    for (i = 0; i < n; ++i)
	for (j = 0; j < n; ++j) {
	    T[i][j] = M[i][j + n];
	}

    DMult(n, n, n, P, T, Resnn);
    Transpose(n, n, P, Q);
    DMult(n, n, n, Resnn, Q, Winv);




//Desalloc


    FreeAlloc2(T, n);
    FreeAlloc2(Q, n);
    FreeAlloc2(P, n);
    FreeAlloc2(Resnn, n);




//construction of the permutation given by the initial data


    unsigned long posin = n - USet[0];
    unsigned long posout = 0;
    unsigned long GPerm[n];
    unsigned long test = 0;
    for (i = 1; i <= n; ++i) {
	for (j = 1; j <= USet[0]; ++j) {
	    if (USet[j] == i - 1) {
		test = 1;
		break;
	    }
	}
	if (test == 1) {
	    GPerm[posin] = i;
	    posin += 1;
	} else {
	    GPerm[posout] = i;
	    posout += 1;
	}
	test = 0;
    }


//Construction of the image by the previous permutation of the set of columns to use in next step

    for (i = 1; i <= Set[0]; ++i) {
	SGPerm[i][0] = GPerm[Set[i]] - 1;
    }
    SGPerm[0][0] = Set[0];




//Construction of Si

    for (i = 0; i < SGPerm[0][0]; ++i)
	for (j = 0; j < n; ++j) {
	    S[i][j] = Si[SGPerm[i + 1][0]][j];
	}


}
