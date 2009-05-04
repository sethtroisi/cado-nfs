#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "random.h"
#include "alloc.h"
#include <limits.h>

#define	WBITS	(CHAR_BIT * sizeof(unsigned long))
#define	iceildiv(x,y)	(((x)+(y)-1)/(y))


/////////////////////////////////////
//
//
// Matrix generation functions
//
//
/////////////////////////////////////


/*

Random element in F_2

*/

FILE *random_file = NULL;

void init_random()
{
    if (random_file != NULL)
	return;
    random_file = fopen("/dev/urandom", "r");
}

void close_random()
{
    if (random_file == NULL)
	return;
    fclose(random_file);
}

unsigned long random_coeff()
{
    if (random_file == NULL)
	init_random();
#if 0
    int b = random() % 2;
#endif
    int b;
    fread(&b, sizeof(unsigned long), 1, random_file);
    return b & 1;
}



unsigned long random_coeff2()
{
    if (random_file == NULL)
	init_random();
#if 0
    int b = random() % 2;
#endif
    unsigned long b;
    fread(&b, sizeof(unsigned long *), 1, random_file);

    return b;
}



void RandomDMatrixBit(unsigned long N,
		      unsigned long m, unsigned long n, unsigned long *a)
{
    int i, j;

    if (n % N == 0)
	i = n / N * m;
    else
	i = n / N * m + m;

    for (j = 0; j < i; ++j) {
	a[j] = random_coeff2();
    }


}

unsigned long *RandomDMatrixBitTest(unsigned long m,
				    unsigned long n, unsigned long *a)
{
    int j;
    for (j = 0; j < iceildiv(n, WBITS) * m; ++j) {
	a[j] = random_coeff2();
    }
    return a;
}


void RandomDMatrixBitTest_new(unsigned long m,
				    unsigned long n,DenseMatrix a)
{
    int j;
    for (j = 0; j < iceildiv(n, WBITS) * m; ++j) {
	a->Data[j] = random_coeff2();
    }
    a->Nrows=m;
    a->Ncols=n;
//    return a;
}



#if 0
void RandomDMatrixS(DenseMatrix a, unsigned long m,unsigned long n)
{
    int j;
    for (j = 0; j < iceildiv(n, WBITS) * m; ++j) {
	a->Data[j] = random_coeff2();
    }
    a->Nrows=m;
    a->Ncols=n;
}

#endif



void CreateRandomMatrixFile(char *f, unsigned long N, unsigned long m,
			    unsigned long n, unsigned long k)
{
    unsigned long *L, *L1;
    unsigned long i, j, pos, pos1, col, t, t1;
    FILE *File;
    L = malloc(N * (k + 20));
    L1 = malloc(N * (k + 20));
    File = fopen(f, "w+");
    for (j = 0; j < m; ++j) {
	pos = 1;
	L[0] = 0;
	pos1 = 0;
	for (i = 0; i <= n; ++i) {
	    col = random_coeff2() % n;
	    if (col > L[L[0]]) {
		L[L[0] + 1] = col;
		pos++;
		L[0] = L[0] + 1;
	    } else {
		if (col < L[L[0]]) {
		    for (t = 1; t < L[0]; ++t) {

			if (L[t] >= col)
			    break;
		    }
		    if (L[t] != col) {
			for (t1 = 0; t1 < t; ++t1) {
			    L1[t1] = L[t1];
			}
			L1[t] = col;
			L[0] = L[0] + 1;
			for (t1 = t + 1; t1 <= L[0]; ++t1) {
			    L1[t1] = L[t1 - 1];
			}

			for (t1 = 1; t1 <= L[0]; ++t1) {
			    L[t1] = L1[t1];
			}
		    }
		}

	    }
	    if ((L[0] > k + (random_coeff2() % 10)) || (L[L[0]] == n))
		break;
	}
	for (i = 0; i <= L[0]; ++i) {
	    fprintf(File, "%lu ", L[i]);
	};
	fprintf(File, "\n");
    }

    free(L);
    free(L1);
    fclose(File);
}


void CreateRandomMatrix(unsigned long m, unsigned long n, unsigned long k,
			unsigned long *a, unsigned long dev)
{
    unsigned long *L, *L1;
    unsigned long i, j, pos, pos1, col, t, t1, LengA;
    L = malloc(sizeof(unsigned long) * (k + dev + 1));
    L1 = malloc(sizeof(unsigned long) * (k + dev + 1));
    LengA = 0;
    for (j = 0; j < m; ++j) {
	pos = 1;
	L[0] = 0;
	pos1 = 0;
	for (i = 0; i <= n; ++i) {
	    col = random_coeff2() % n;
	    if (col > L[L[0]]) {
		L[L[0] + 1] = col;
		pos++;
		L[0] = L[0] + 1;
	    } else {
		if (col < L[L[0]]) {
		    for (t = 1; t < L[0]; ++t) {

			if (L[t] >= col)
			    break;
		    }
		    if (L[t] != col) {
			for (t1 = 0; t1 < t; ++t1) {
			    L1[t1] = L[t1];
			}
			L1[t] = col;
			L[0] = L[0] + 1;
			for (t1 = t + 1; t1 <= L[0]; ++t1) {
			    L1[t1] = L[t1 - 1];
			}

			for (t1 = 1; t1 <= L[0]; ++t1) {
			    L[t1] = L1[t1];
			}
		    }
		}

	    }
	    if ((L[0] > k + (random_coeff2() % dev)) || (L[L[0]] == n))
		break;
	}


	for (i = 0; i <= L[0]; ++i)
	    a[LengA + i] = L[i];
	LengA += L[0] + 1;
    }
    free(L);
    free(L1);
}





void CreateRandomSymmetricMatrix(unsigned long n, unsigned long k,
				 unsigned long *A, unsigned long dev)
{
    unsigned long *L;
    unsigned long *L1;

//unsigned long *A;

    unsigned long i, j, col, t, t1, LengA, SizeRowA;

    SizeRowA = k + dev + 1;

    L = malloc(sizeof(unsigned long) * (SizeRowA + 1));
    L1 = malloc(sizeof(unsigned long) * (SizeRowA + 1));

//A=malloc(sizeof(unsigned long)*SizeRowA*n);

    memset(A, 0, sizeof(unsigned long) * SizeRowA * n);

    LengA = 0;

    for (j = 0; j < n; ++j) {


	for (t = 0; t <= A[j * SizeRowA]; ++t)
	    L[t] = A[j * SizeRowA + t];


	for (i = A[j * SizeRowA]; i <= n; ++i) {
	    col = (random_coeff2() % n);


	    if ((A[col * SizeRowA] <= k + (random_coeff2() % dev))
		&& (col >= j)) {

		if (col > L[L[0]]) {
		    L[L[0] + 1] = col;
		    L[0] = L[0] + 1;

		    A[col * SizeRowA + A[col * SizeRowA] + 1] = j;
		    A[col * SizeRowA] = A[col * SizeRowA] + 1;
		} else {
		    if (col < L[L[0]]) {
			for (t = 1; t < L[0]; ++t) {
			    if (L[t] >= col)
				break;
			}
			if (L[t] != col) {
			    for (t1 = 0; t1 < t; ++t1) {
				L1[t1] = L[t1];
			    }
			    L1[t] = col;
			    L[0] = L[0] + 1;
			    for (t1 = t + 1; t1 <= L[0]; ++t1) {
				L1[t1] = L[t1 - 1];
			    }
			    for (t1 = 1; t1 <= L[0]; ++t1) {
				L[t1] = L1[t1];
			    }
			    A[col * SizeRowA + A[col * SizeRowA] + 1] = j;
			    A[col * SizeRowA] = A[col * SizeRowA] + 1;
			}
		    }
		}

	    }
	    if ((L[0] > k + (random_coeff2() % dev)) || (L[L[0]] == n))
		break;
	}

	for (i = 0; i <= L[0]; ++i)
	    A[LengA + i] = L[i];
	LengA += L[0] + 1;
    }

    free(L);

    free(L1);



}

unsigned long *CreateRandomSymmetricMatrixTest(unsigned long m,
					       unsigned long n,
					       unsigned long k)
{
    unsigned long *a;
    unsigned long *L;
    unsigned long *L1;
    unsigned long *A;
    unsigned long i, j, col, t, t1, LengA, SizeRowA, VAR;



    a = malloc(m * (k + 10) * sizeof(unsigned long));

    VAR = 11;

    SizeRowA = k + VAR;

    L = malloc(sizeof(unsigned long) * (k + VAR));
    L1 = malloc(sizeof(unsigned long) * (k + VAR));
    A = malloc(sizeof(unsigned long) * (k + VAR) * m);


    for (i = 0; i < (k + VAR) * m; ++i)
	A[i] = 0;

    LengA = 0;


    for (j = 0; j < m; ++j) {


	for (t = 0; t <= A[j * SizeRowA]; ++t)
	    L[t] = A[j * SizeRowA + t];


	for (i = A[j * SizeRowA]; i <= n; ++i) {
	    col = (random_coeff2() % n);


	    if ((A[col * SizeRowA] <= k + (random_coeff2() % 10))
		&& (col >= j)) {

		if (col > L[L[0]]) {
		    L[L[0] + 1] = col;
		    L[0] = L[0] + 1;

		    A[col * SizeRowA + A[col * SizeRowA] + 1] = j;
		    A[col * SizeRowA] = A[col * SizeRowA] + 1;
		} else {
		    if (col < L[L[0]]) {
			for (t = 1; t < L[0]; ++t) {
			    if (L[t] >= col)
				break;
			}
			if (L[t] != col) {
			    for (t1 = 0; t1 < t; ++t1) {
				L1[t1] = L[t1];
			    }
			    L1[t] = col;
			    L[0] = L[0] + 1;
			    for (t1 = t + 1; t1 <= L[0]; ++t1) {
				L1[t1] = L[t1 - 1];
			    }
			    for (t1 = 1; t1 <= L[0]; ++t1) {
				L[t1] = L1[t1];
			    }
			    A[col * SizeRowA + A[col * SizeRowA] + 1] = j;
			    A[col * SizeRowA] = A[col * SizeRowA] + 1;
			}
		    }
		}

	    }
	    if ((L[0] > k + (random_coeff2() % 10)) || (L[L[0]] == n))
		break;
	}

	for (i = 0; i <= L[0]; ++i)
	    a[LengA + i] = L[i];
	LengA += L[0] + 1;
    }

    free(L);
    free(L1);
    free(A);
    return a;
}



////////////////////////////////////////////////////
//
//  To older versions of Lanczos
//
/////////////////////////////////////////////////

/*
Generates a dense m x n matrix
*/


void RandomDMatrix(unsigned long m, unsigned long n, unsigned long **a)
{
    unsigned long i, j;
    for (i = 0; i < m; ++i)
	for (j = 0; j < n; ++j) {
	    a[i][j] = random_coeff();
	}
//return 0;
}


/*
Generates a sparce m x n matrix
*/


void RandomSMatrix(unsigned long m, unsigned long n, unsigned long **a)
{
    unsigned long i, j;
    for (i = 0; i < m; ++i) {

	a[i][0] = 0;
	for (j = 1; j < n + 1; j++) {
	    if (random_coeff() == 1) {
		a[i][0] = a[i][0] + 1;
		a[i][a[i][0]] = j - 1;

	    }
	}
    }
    a[m][0] = -1;

//return 0;
}



/*
Generates a symmetric sparce n x n matrix
*/
/*

void RandomSymmetricSMatrix(unsigned long n, unsigned long **a)
{
	unsigned long i, j, l;
	for (i = 0; i < n; ++i) {
		a[i][0] = 0;
	}

	for (i = 0; i < n; ++i) {
		l = a[i][a[i][0]] + 1;
		for (j = l + 1; j < n; j++) {
			if (random_coeff() == 1) {
				if (i != j) {
					a[i][0] = a[i][0] + 1;
					a[i][a[i][0]] = j;
					a[j][0] = a[j][0] + 1;
					a[j][a[j][0]] = i;
				}


			}
		}
	}
	a[n][0] = -1;


}

*/


/*
Generates a symmetric dense n x n matrix
*/


void RandomSymetricDMatrix(unsigned long n, unsigned long **a)
{
    unsigned long i, j;
    for (i = 0; i < n; ++i)
	for (j = i; j < n; ++j) {
	    a[i][j] = random_coeff();
	    a[j][i] = a[i][j];
	}


}
