#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "alloc.h"


#define	iceildiv(x,y)	(((x)+(y)-1)/(y))
#define	WBITS	(CHAR_BIT * sizeof(unsigned long))


//////////////////////////////////////
//
//    Memory  Allocation
//   
//
//////////////////////////////////////


unsigned long *Allocmn(unsigned long m, unsigned long n)
{
    unsigned long limbs_per_matrix;
    unsigned long *a;

    limbs_per_matrix = iceildiv(n, WBITS) * m;
    limbs_per_matrix += limbs_per_matrix&1;

    a = malloc(limbs_per_matrix * sizeof(unsigned long));

    if (a == NULL) {
	fprintf(stderr, "out of memory\n");
	return 0;
    }
    return a;
}

unsigned long **Allocmn3(unsigned long m, unsigned long n)
{
    unsigned long **a;

    a = malloc(3 * sizeof(unsigned long *));
    a[0] = Allocmn(m, n);
    a[1] = Allocmn(m, n);
    a[2] = Allocmn(m, n);
    return a;
}



void FreeAllocmn3(unsigned long **a)
{
    unsigned long i;

    for (i = 0; i < 3; ++i)
	free(a[i]);
    free(a);
}

void FreeAllocmn3_new(DenseMatrix a[3])
{
    unsigned long i;

    for (i = 0; i < 3; ++i){
	free(a[i]->Data);}
    
}


///////////////////////////////////////////
//
// To older versions of Lanczos
//
//////////////////////////////////////////


unsigned long ***Alloc3Test(unsigned long ***a, unsigned long nrows,
			    unsigned long ncolumns)
{
    unsigned long k, i;
    a = malloc(3 * sizeof(int *));
    if (a == NULL) {
	fprintf(stderr, "out of memory\n");
	return 0;
    }

    for (k = 0; k < 3; ++k) {
	a[k] = malloc(nrows * sizeof(int));
	if (a[k] == NULL) {
	    fprintf(stderr, "out of memory\n");
	    return 0;
	}
    }
    for (k = 0; k < 3; ++k) {
	for (i = 0; i < nrows; i++) {
	    a[k][i] = malloc(ncolumns / 8);
	    if (a[k][i] == NULL) {
		fprintf(stderr, "out of memory\n");
		return 0;
	    }
	}
    }
    return a;
}


unsigned long **Alloc2Test(unsigned long **a, unsigned long nrows,
			   unsigned long ncolumns)
{
    unsigned long i;
    a = malloc(nrows * sizeof(int *));
    if (a == NULL) {
	fprintf(stderr, "out of memory\n");
	return 0;
    }
    for (i = 0; i < nrows; i++) {
	a[i] = malloc(ncolumns / 8);
	if (a[i] == NULL) {
	    fprintf(stderr, "out of memory\n");
	    return 0;
	}
    }

    return a;
}

void FreeAlloc2(unsigned long **a, unsigned long nrows)
{
    unsigned long i;

    for (i = 0; i < nrows; ++i)
	free(a[i]);
    free(a);
}



void FreeAlloc3(unsigned long ***a, unsigned long nrows)
{
    unsigned long k, i;
    for (k = 0; k < 3; ++k)
	for (i = 0; i < nrows; ++i)
	    free(a[k][i]);
    for (k = 0; k < 3; ++k)
	free(a[k]);
    free(a);
}




unsigned long ***Alloc3(unsigned long nrows, unsigned long ncolumns)
{
    unsigned long ***a;
    unsigned long k, i;
    a = malloc(3 * sizeof(unsigned long **));
    if (a == NULL) {
	fprintf(stderr, "out of memory\n");
	return 0;
    }

    for (k = 0; k < 3; ++k) {
	a[k] = malloc(nrows * sizeof(unsigned long *));
	if (a[k] == NULL) {
	    fprintf(stderr, "out of memory\n");
	    return 0;
	}
    }
    for (k = 0; k < 3; ++k) {
	for (i = 0; i < nrows; i++) {
	    a[k][i] = malloc(ncolumns * sizeof(unsigned long));
	    if (a[k][i] == NULL) {
		fprintf(stderr, "out of memory\n");
		return 0;
	    }
	}
    }
    return a;
}


unsigned long **Alloc2(unsigned long nrows, unsigned long ncolumns)
{
    unsigned long **a;
    unsigned long i;
    a = malloc(nrows * sizeof(unsigned long *));
    if (a == NULL) {
	fprintf(stderr, "out of memory\n");
	return 0;
    }
    for (i = 0; i < nrows; i++) {
	a[i] = malloc(ncolumns * sizeof(unsigned long));
	if (a[i] == NULL) {
	    fprintf(stderr, "out of memory\n");
	    return 0;
	}
    }

    return a;
}




unsigned long *AllocmnOld(unsigned long N, unsigned long m, unsigned long n,
			  unsigned long *a)
{
    unsigned long i;

    if (n % N == 0)
	i = n / N * m;
    else
	i = (n / N + 1) * m;

// printf("in allocmn %lu  %lu %lu %lu \n",i,m,n,(n/N+1)*m);

    a = malloc(WBITS * i);

// printf("in 2 allocmn %lu\n",i);
    if (a == NULL) {
	fprintf(stderr, "out of memory\n");
	return 0;
    }
    return a;
}



unsigned long **Allocmn3Old(unsigned long N, unsigned long m, unsigned long n,
			    unsigned long **a)
{
    unsigned long i, j;

    if (n % N == 0)
	i = n / N * m;
    else
	i = (n / N + 1) * m;

    a = malloc(WBITS * 3);
    if (a == NULL) {
	fprintf(stderr, "out of memory\n");
	return 0;
    }

    for (j = 0; j < 3; ++j) {

	a[j] = malloc(WBITS * i);
	if (a[j] == NULL) {
	    fprintf(stderr, "out of memory\n");
	    return 0;
	}

    }
    return a;
}
