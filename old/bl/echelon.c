#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"
#include "mat_ops.h"
#include "echelon.h"
#include "ReadWrite.h"
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







void ConstructPermBit_new(unsigned long n, unsigned long *US, DenseMatrix p)
{
    unsigned long i, j, k, SizeRowA;
    SizeRowA = iceildiv(n, WBITS);
    DenseMatrix IN,Col,p1;
    IN->Data = Allocmn(n, n);IN->Nrows=n;IN->Ncols=n;
    Col->Data = Allocmn(1, n);Col->Nrows=1;Col->Ncols=n;
    p1->Data = Allocmn(n, n);p1->Nrows=n;p1->Ncols=n;
    memset(IN->Data, 0, n * iceildiv(n, WBITS) * sizeof(unsigned long));
    for (i = 0; i < n; ++i) {
	IN->Data[i * SizeRowA + i / WBITS] = (1UL << i % WBITS);
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
	SelectColumn_new(IN, i, Col);
	if (test == 1) {
	    for (k = 0; k < SizeRowA; ++k) {
		p1->Data[posin * SizeRowA + k] = Col->Data[k];
	    }
	    posin += 1;
	} else {
	    for (k = 0; k < SizeRowA; ++k) {
		p1->Data[posout * SizeRowA + k] = Col->Data[k];
	    }
	    posout += 1;
	}
	test = 0;
    }
    TransposeBit_new(p1, p);
    free(IN->Data);
    free(Col->Data);
    free(p1->Data);
}








void Construction_S_WinvBit_new(DenseMatrix m,
			    unsigned long *USet, DenseMatrix S,
			    unsigned long *Set0, DenseMatrix Winv,
			    unsigned long *SGPerm)
{
    unsigned long i, j, posS, k, h, SizeRowA;
    unsigned long n=m->Nrows;

    DenseMatrix P,Q,T,Res;
    Res->Data=Allocmn(n,n);
    P->Data = Allocmn(n,n);P->Nrows=n;P->Ncols=n;
    Q->Data = Allocmn(n,n);Q->Nrows=n;Q->Ncols=n;
    T->Data = Allocmn(n,n);T->Nrows=n;T->Ncols=n;
    Res->Nrows=n;
    Res->Ncols=n;
    ConstructPermBit_new(n,USet, P);
    TransposeBit_new(P, Q);
    DMultBit_new(Q, m,Res);
    DMultBit_new(Res, P, T);

    DenseMatrix M;
    M->Data = Allocmn(n, 2 * n);M->Nrows=n;M->Ncols=2*n;
    memset(M->Data, 0, n * iceildiv(2 * n, WBITS) * sizeof(unsigned long));
    SizeRowA = iceildiv(n, WBITS);
    DenseMatrix IN;
    IN->Data = Allocmn(n, n);IN->Nrows=n;IN->Ncols=n;
    memset(IN->Data, 0, n * iceildiv(n, WBITS) * sizeof(unsigned long));
    for (i = 0; i < n; ++i) {
	IN->Data[i * SizeRowA + i / WBITS] = (1UL << i % WBITS);
    }
    unsigned long SizeRow2A;
    SizeRow2A = iceildiv(n + n, WBITS);
    for (i = 0; i < n; ++i) {
	for (j = 0; j < SizeRow2A; ++j) {

	    if (j < SizeRowA - 1) {

		M->Data[i * SizeRow2A + j] = T->Data[i * SizeRowA + j];
	    }
	    if (j == SizeRowA - 1) {

		if ((n % WBITS > i) && (j - SizeRowA + 1 == i / WBITS))
		    M->Data[i * SizeRow2A + SizeRowA - 1] = (1UL << (i));
		M->Data[i * SizeRow2A + SizeRowA - 1] =
		    (M->Data[i * SizeRow2A + SizeRowA - 1] << n % WBITS);
		M->Data[i * SizeRow2A + SizeRowA - 1] |=
		    T->Data[i * SizeRowA + SizeRowA - 1];
	    }
	    if (j > SizeRowA - 1) {

		if ((i < WBITS - n % WBITS) && (n % WBITS != 0))
		    M->Data[i * SizeRow2A + j] =
			(1UL << (n % WBITS - i % WBITS - 1));
		else {
		    if ((n % WBITS) != 0) {
			M->Data[i * SizeRow2A + j] =
			    (1UL << (i % WBITS - WBITS + n % WBITS));
		    } else if (j - SizeRowA == i / WBITS) {
			M->Data[i * SizeRow2A + j] = (1UL << (i));
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
	    if ((M->Data[j * SizeRow2A + j / WBITS] >> j & 1) == 0) {
		for (k = j; k < n; ++k) {
		    if ((M->Data[k * SizeRow2A + j / WBITS] >> j & 1) != 0) {
			SwapRowsBit_new(M, j, k);
			break;
		    }
		}
	    }
	}
	if ((M->Data[j * SizeRow2A + j / WBITS] >> j & 1) != 0) {
	    Set[posS] = j;
	    ++posS;

	    MultiplyRowBit_new( M,
			   (M->Data[j * SizeRow2A + j / WBITS] >> j & 1), j);

	    for (h = 0; h < n; ++h) {
		if (h != j) {
		    AddRowBit_new( M,
			      (M->Data[h * SizeRow2A + j / WBITS] >> j & 1), j, h);
		}
	    }
	} else {
	    for (k = j; k < n; ++k) {
		if ((M->Data[j * SizeRow2A + (j + n) / WBITS] >> (j + n) & 1) == 0) {
		    if ((M->Data[k * SizeRow2A + (j + n) / WBITS] >> (j + n) & 1) !=
			0) {
			SwapRowsBit_new( M, j, k);
			break;
		    }
		}
	    }
	    if ((M->Data[j * SizeRow2A + (j + n) / WBITS] >> (j + n) & 1) != 0) {
		for (h = 0; h < n; ++h) {
		    if (h != j) {
			AddRowBit_new( M,
				  (M->Data[h * SizeRow2A + (j + n) / WBITS] >>
				   (j + n) & 1), j, h);
		    }
		}
		MultiplyRowBit_new(M, 0, j);
	    }
	}

    }
    DenseMatrix Si;
    Si->Data = Allocmn(n, n);Si->Nrows=n;Si->Ncols=n;
    memset(Si->Data, 0, n * iceildiv(n, WBITS) * sizeof(unsigned long));
    for (i = 0; i < n; ++i) {
	Si->Data[i * SizeRowA + i / WBITS] = (1UL << i % WBITS);
    }
    Set[0] = posS - 1;
    Set0[0] = posS - 1;

//construction of Winv corrected by permutation

    memset(IN->Data, 0, n * iceildiv(n, WBITS) * sizeof(unsigned long));
    DenseMatrix C;
    C->Data = Allocmn(n, n);C->Nrows=n;C->Ncols=n;
    for (i = 0; i < n; ++i) {
	SelectColumn_new(M, n + i, C);
	for (j = 0; j < SizeRowA; ++j)
	    IN->Data[i * SizeRowA + j] = C->Data[j];
    }

    TransposeBit_new(IN, T);
    DMultBit_new(P,T, Res);
    TransposeBit_new(P, Q);
    DMultBit_new(Res, Q, Winv);

//Desalloc

    free(T->Data);
    free(Q->Data);
    free(P->Data);
    free(Res->Data);
    free(M->Data);
    free(IN->Data);


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
	    S->Data[i * SizeRowA + j] = Si->Data[SGPerm[i + 1] * SizeRowA + j];
	}

    free(Set);
    free(GPerm);
    free(C->Data);
    free(Si->Data);

}








/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
























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
