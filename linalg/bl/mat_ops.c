#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include "mat_ops.h"
#include "alloc.h"
#include "ReadWrite.h"
#include "alloc.h"

#define	WBITS	(CHAR_BIT * sizeof(unsigned long))
#define	iceildiv(x,y)	(((x)+(y)-1)/(y))




//////////////////////////////////////
//
//
// Operations with matrices
//
//
//////////////////////////////////////





///////////////////////////////////////////////////////////////////////
//
//  Matrix operations on lines in Bit notation
//    
//
//////////////////////////////////////////////////////////////////



unsigned long GaussElimBit(unsigned long m, unsigned long n, unsigned long *a,
			   unsigned long *b, unsigned long *c,
			   unsigned long *d, unsigned long *Lend)
{
    unsigned long teste, Count, SizeRowA, i, j, k, t, Posji, Poskt, Pivotit;
    SizeRowA = iceildiv(n, WBITS);
    for (i = 0; i < m * SizeRowA; i++)
	return 0;
}








/*
unsigned long GaussElimBit(unsigned long m, unsigned long n, unsigned long *a, unsigned long *b, unsigned long *c, unsigned long *d,unsigned long *Lend)
{
	unsigned long teste,Count,SizeRowA,i,j,k,t,Posji,Poskt,Pivotit;
	SizeRowA=iceildiv(n,WBITS);
	unsigned long *p;
	p=Allocmn(m,m);
	unsigned long *Resnm;
	Resnm=Allocmn(n,m);
	teste=0;
	memset(c,0,m*iceildiv(m,WBITS)*sizeof(unsigned long));
        //memset(c1,0,m*iceildiv(m,WBITS)*sizeof(unsigned long));
        
        for (i = 0; i <m ; ++i)
	{
         c[i*iceildiv(m,WBITS)+i/WBITS]=(1UL<<i%WBITS);
         //c1[i*iceildiv(m,WBITS)+i/WBITS]=(1UL<<i%WBITS);
	}
        for (i=0; i<m*SizeRowA;i++) {b[i]=a[i];}
        
        Count=0;
        for (i=0; i<m; ++i)
        {t=i;
        do {
            Pivotit=((b[i*SizeRowA+ t / WBITS]>> (t % WBITS)) & 1UL);
            k=0;
            if (Pivotit==0)
        	{
        	for (k=i+1; k<m; ++k)
       		{
          	Poskt=((b[k*SizeRowA+ t / WBITS]>> (t % WBITS)) & 1UL);
          	if (Poskt==1) 
            	{
            	SwapRowsBit(m,n,b,i,k);
            	SwapRowsBit(m,m,c,i,k);
            	break;
            	}  
        	}
        	}        	
            if (k==m) 
            {
            t++;
            }
            //printf("toto %lu %lu\n",t,k);
           } while (k==m); 	
        if (t>n && teste==0) {Lend[0]=i;teste=1;}      	                
        for (j=i+1; j<m; ++j)
        {
        Posji=((b[j*SizeRowA+ i / WBITS]>> (i % WBITS)) & 1UL);
        if (Posji==1)
        {
        AddRowBit(m,n,b,1,i,j); 
        AddRowBit(m,m,c,1,i,j);
        }
        
        }
        }
        k=0;
        for (i=0; i<m-Lend[0];i++) 
        {
        
        d[i]=i+Lend[0];
        }
        free(p);
        free(Resnm);
	return 0;
}
*/

unsigned long InverseList(unsigned long N, unsigned long Lengm,
			  const unsigned long *m, unsigned long *b)
{
    unsigned long j, i;
    unsigned long k = 0;
    for (i = 0; i < N; i++) {
	for (j = 0; j < Lengm; j++) {
	    if (m[j] == i) {
		break;
	    }
	}
	if (j >= Lengm) {
	    b[k] = i;
	    k++;
	}
    }
    return 0;
}

unsigned long SelectLinesBit(unsigned long n, const unsigned long *a,
			     unsigned long m, unsigned long *b)
{
    unsigned long SizeRowA, i, j;

    SizeRowA = iceildiv(n, WBITS);

    for (i = 0; i < m; ++i)
	for (j = 0; j < SizeRowA; ++j)
	    b[i * SizeRowA + j] = a[i * SizeRowA + j];
    return 0;
}

unsigned long SelectLinesListBit(unsigned long n, const unsigned long *a,
				 unsigned long *m, unsigned long Lengm,
				 unsigned long *b)
{
    unsigned long SizeRowA, i, j;

    SizeRowA = iceildiv(n, WBITS);

    for (i = 0; i < Lengm; ++i) {
	for (j = 0; j < SizeRowA; ++j)
	    b[i * SizeRowA + j] = a[m[i] * SizeRowA + j];
    }
    return 0;
}





unsigned long AddRowBit(unsigned long m, unsigned long n, unsigned long *a,
			unsigned long u, unsigned long i, unsigned long j)
{
    unsigned long p, SizeRowA;
    SizeRowA = iceildiv(n, WBITS);

    if (u != 0) {


	for (p = 0; p < SizeRowA; ++p) {
	    a[j * SizeRowA + p] ^= a[i * SizeRowA + p];
	}
    }
    return 0;
}



unsigned long SwapRowsBit(unsigned long m, unsigned long n, unsigned long *a,
			  unsigned long i, unsigned long j)
{
    unsigned long k, p, SizeRowA;

    SizeRowA = iceildiv(n, WBITS);

    for (p = 0; p < SizeRowA; ++p) {
	k = a[i * SizeRowA + p];
	a[i * SizeRowA + p] = a[j * SizeRowA + p];
	a[j * SizeRowA + p] = k;
    }
    return 0;
}





unsigned long MultiplyRowBit(unsigned long m, unsigned long n,
			     unsigned long *a, unsigned long u,
			     unsigned long j)
{
    unsigned long p, SizeRowA;

    SizeRowA = iceildiv(n, WBITS);

    if (u == 0) {
	for (p = 0; p < SizeRowA; ++p) {
	    a[j * SizeRowA + p] = 0;
	}
    }
    return 0;
}




/*

Dense Matrix + Dense Matrix = Dense Matrix

*/

void DMatrixSumBit(unsigned long m,
		   unsigned long n,
		   const unsigned long *A,
		   const unsigned long *B, unsigned long *C)
{
    unsigned long i, SizeA;

    SizeA = iceildiv(n, WBITS) * m;

    for (i = 0; i < SizeA; i++) {
	C[i] = A[i] ^ B[i];
    }
}


/*

Dense Matrix x Dense Matrix = Dense Matrix

*/



void DMultBit(unsigned long m,
	      unsigned long n,
	      unsigned long p,
	      const unsigned long *A,
	      const unsigned long *B, unsigned long *C)
{
    unsigned long i, j, t;
    unsigned long LrowA, LrowB;

    LrowA = iceildiv(n, WBITS);
    LrowB = iceildiv(p, WBITS);

    const unsigned long *A1 = A;

    for (i = 0; i < m; i++) {

	unsigned long P;
	unsigned int k, count;

	for (t = 0; t < LrowB; t++) {

	    const unsigned long *sptr = B;
	    count = 0;
	    C[i * LrowB + t] = 0;
	    sptr += t;
	    for (j = 0; j < n / WBITS; j++) {
		P = A1[j];
		for (k = 0; k < WBITS; k++) {

		    if (P & 1UL)
			C[i * LrowB + t] ^= *sptr;



		    sptr += LrowB;
		    P >>= 1UL;
		}
	    }

	    if (n % WBITS) {
		P = A1[j];
		for (k = 0; k < n % WBITS; k++) {
		    if (P & 1UL)
			C[i * LrowB + t] ^= *sptr;



		    sptr += LrowB;
		    P >>= 1UL;
		}
	    }

	}

	A1 += LrowA;
    }
}




void TVUBit(unsigned long m,
	    unsigned long n,
	    const unsigned long *A, const unsigned long *B, unsigned long *C)
{
    unsigned long i, P, k;

    memset(C, 0, WBITS * sizeof(unsigned long));

    for (i = 0; i < m; i++) {

	P = *A++;
	for (k = 0; k < WBITS; k++) {

	    //if (P & 1UL) C[k] ^= B[i];

	    C[k] ^= (B[i] & -(P & 1UL));

	    P >>= 1UL;
	}

    }


}



void VUBit(unsigned long m,
	   unsigned long n,
	   const unsigned long *A, const unsigned long *B, unsigned long *C)
{
    unsigned long i, P, k;

    memset(C, 0, m * sizeof(unsigned long));

    for (i = 0; i < m; i++) {

	P = *A++;
	for (k = 0; k < n; k++) {

	    //if (P & 1UL) C[k] ^= B[i];

	    C[i] ^= (B[k] & -(P & 1UL));

	    P >>= 1UL;
	}

    }


}





inline void addmul(uint64_t a, uint64_t b, uint64_t * w)
{
    unsigned int i;
#if 0
    /* Dans un sens */
    for (i = 0; i < 64; i++) {
	*w++ ^= b & -(a & 1);
	a >>= 1;
    }
#endif
#if 0
    /* Dans l'autre -- va plus vite. */
    for (i = 0; i < 64; i++) {
	*w++ ^= b & (((int64_t) a) >> 63);
	a <<= 1;
    }
#endif
#if 0
    /* Ã€ peu prÃ¨s comme la mÃ©thode 1, mais pas mieux */
    typedef uint64_t mvec_t[2];
    mvec_t mb[4] = {
	{0, 0}, {b, 0}, {0, b}, {b, b},
    };
    for (i = 0; i < 64; i += 2) {
	const uint64_t *y = mb[a & 3];
	*w++ ^= y[0];
	*w++ ^= y[1];
	a >>= 2;
    }
#endif
#if 1
    /* Avec des sse-2 */
    typedef uint64_t sse2_t __attribute__ ((vector_size(16)));
    sse2_t mb[4] = {
	(sse2_t) {0, 0},
	(sse2_t) {b, 0},
	(sse2_t) {0, b},
	(sse2_t) {b, b},
    };
    sse2_t *sw = (sse2_t *) w;
    for (i = 0; i < 64; i += 2) {
	*sw++ ^= mb[a & 3];
	a >>= 2;
    }
#endif
}






void TVUBit_v2(unsigned long m,
	       unsigned long n, uint64_t * A, uint64_t * B, uint64_t * C)
{
    unsigned long i, P, j;

    memset(C, 0, WBITS * sizeof(unsigned long));

    for (j = 0; j < 16; j++) {
	unsigned int i;
	memset(C, 0, 64 * sizeof(uint64_t));
	for (i = 0; i < m; i++) {
	    addmul(A[i], B[i], C);
	}
    }

}


















/*

Sparse Matrix x Dense Matrix = Dense Matrix

*/





/*
void SMultDmatrixBitNew(unsigned long NrLines,unsigned long NrCols_b, unsigned long *a, unsigned long *b, unsigned long *c)
{
	unsigned long k, j, i, s,count,BgLine,SizeRowsb,t;
         
       
        SizeRowsb=iceildiv(NrCols_b,WBITS);
        memset(c,0,NrLines*SizeRowsb*sizeof(unsigned long));
        
        BgLine=0;count=0;
        i=0;
        for (i=0;i<NrLines; ++i)
	{k=a[BgLine];
            s = 0;
	      for (j = 1; j <= k; ++j)
		{for (t = 0; t < SizeRowsb; ++t)
			c[i*SizeRowsb+t]^=b[a[BgLine+j]*SizeRowsb+t];count++;
		}
              BgLine+=k+1;
             
          } 
}

*/




/*

Select column k from dense matrix V

*/



void SelectColumn(unsigned long m,
		  unsigned long n,
		  const unsigned long *V, unsigned long k, unsigned long *C)
{
    unsigned long i, step, L;


    L = iceildiv(m, WBITS);

    step = 0;
//for (i = 0; i <L ; ++i) C[i]=0;
    memset(C, 0, L * sizeof(unsigned long));
    for (i = 0; i < m; ++i) {
	if (n % WBITS == 0) {
	    if (((V[(n / WBITS) * i + k / WBITS] >> (k % WBITS)) & 1) != 0) {
		C[i / WBITS] |= (1UL << step);

	    }
	} else {
	    if (((V[(n / WBITS + 1) * i + k / WBITS] >> (k % WBITS)) & 1) !=
		0) {
		C[i / WBITS] |= (1UL << step);
	    }
	}
	if (((i + 1) / WBITS) > (i / WBITS))

	    step = 0;

	else
	    step++;
    }
}




/*

Transpose dense matrix

*/

void TransposeBit(unsigned long m, unsigned long n, const unsigned long *a,
		  unsigned long *b)
{
    unsigned long i, j, SizeColumnB;


    SizeColumnB = iceildiv(m, WBITS);

    unsigned long *C;
    C = malloc(SizeColumnB * sizeof(unsigned long));

    for (i = 0; i < n; ++i) {

	SelectColumn(m, n, a, i, C);
	for (j = 0; j < SizeColumnB; ++j) {
	    b[i * SizeColumnB + j] = C[j];
	}
    }
    free(C);
}

/*

dense plus dense matrix

*/

void DSumBit(unsigned long m, unsigned long n, unsigned long *a,
	     unsigned long *b, unsigned long *c)
{
    unsigned long k, j, SizeRowA;
    SizeRowA = iceildiv(n, WBITS);

    for (k = 0; k < m; ++k)
	for (j = 0; j < SizeRowA; ++j) {
	    c[k * SizeRowA + j] = (a[k * SizeRowA + j] ^ b[k * SizeRowA + j]);
	}


}

/*

Transpose sparse times sparse times dense block

*/


void TSSMultDmatrix(unsigned long m, unsigned long n, unsigned long p,
		    unsigned long *a, unsigned long *b, unsigned long *c)
{
    unsigned long k, j, i, s, BgLine, SizeRowsb, t;
    unsigned long *Resmp;
    Resmp = Allocmn(m, p);

    SizeRowsb = iceildiv(p, WBITS);

    memset(c, 0, n * SizeRowsb * sizeof(unsigned long));
    memset(Resmp, 0, m * SizeRowsb * sizeof(unsigned long));

    BgLine = 0;
    i = 0;
    for (i = 0; i < m; ++i) {
	k = a[BgLine];
	s = 0;
	for (j = 1; j <= k; ++j) {
	    for (t = 0; t < SizeRowsb; ++t)
		Resmp[i * SizeRowsb + t] ^= b[a[BgLine + j] * SizeRowsb + t];
	}
	BgLine += k + 1;
    }

    BgLine = 0;
    i = 0;
    for (i = 0; i < m; ++i) {
	k = a[BgLine];
	s = 0;
	for (j = 1; j <= k; ++j) {
	    for (t = 0; t < SizeRowsb; ++t)
		c[a[BgLine + j] * SizeRowsb + t] ^= Resmp[i * SizeRowsb + t];
	}
	BgLine += k + 1;
    }



    free(Resmp);
}



void SMultDmatrixBit(unsigned long m, unsigned long n, unsigned long p,
		     unsigned long *a, unsigned long *b, unsigned long *c,
		     unsigned long start, unsigned long size)
{
/*	unsigned long k, j, i, s,BgLine,SizeRowsb,t;
	       
        SizeRowsb=iceildiv(p,WBITS);
        
        memset(c,0,m*SizeRowsb*sizeof(unsigned long));
        
        BgLine=0;
        i=0;
        for (i=start;i<start+size; ++i)
	{k=a[BgLine];
            s = 0;
	      for (j = 1; j <= k; ++j)
		{for (t = 0; t < SizeRowsb; ++t)
			c[i*SizeRowsb+t]^=b[a[BgLine+j]*SizeRowsb+t];
		}
              BgLine+=k+1;  
         }
       
 */
    unsigned long k, i;

    //assert(SizeRowsb == 1);

    memset(c, 0, m * sizeof(unsigned long));

    i = 0;
    for (i = start; i < start + size; ++i) {
	for (k = *a++; k--;) {
	    c[i] ^= b[*a++];
	}
    }



}


void SMultDmatrixBitNew(unsigned long NrLines, unsigned long NrCols_b,
			unsigned long *a, unsigned long *b, unsigned long *c)
{
#if 0
    unsigned long k, j, i, s, count, BgLine, SizeRowsb, t;
    SizeRowsb = iceildiv(NrCols_b, WBITS);
    memset(c, 0, NrLines * SizeRowsb * sizeof(unsigned long));

    BgLine = 0;
    count = 0;
    i = 0;
    for (i = 0; i < NrLines; ++i) {
	k = a[BgLine];
	s = 0;
	for (j = 1; j <= k; ++j) {
	    for (t = 0; t < SizeRowsb; ++t)
		c[i * SizeRowsb + t] ^= b[a[BgLine + j] * SizeRowsb + t];
	    count++;
	}
	BgLine += k + 1;

    }
#endif

    unsigned long k, i;

    //assert(SizeRowsb == 1);
    memset(c, 0, NrLines * sizeof(unsigned long));

    i = 0;
    for (i = 0; i < NrLines; ++i) {
	for (k = *a++; k--;) {
	    c[i] ^= b[*a++];
	}
    }
}







void TSMultDmatrixBit(unsigned long m, unsigned long n, unsigned long p,
		      unsigned long *a, unsigned long *b, unsigned long *c,
		      unsigned long start, unsigned long size)
{
    /*unsigned long k, j, i, s,BgLine,SizeRowsb,t;

       SizeRowsb=iceildiv(p,WBITS);

       memset(c,0,n*SizeRowsb*sizeof(unsigned long));


       BgLine=0;
       i=0;
       for (i=start;i<start+size; ++i)
       {k=a[BgLine];
       s = 0;
       for (j = 1; j <= k; ++j)
       {for (t = 0; t < SizeRowsb; ++t)
       c[a[BgLine+j]*SizeRowsb+t]^=b[i*SizeRowsb+t];
       }
       BgLine+=k+1;  
       }
     */

    unsigned long k, i;
    memset(c, 0, n * sizeof(unsigned long));

    i = 0;
    for (i = start; i < start + size; ++i) {
	for (k = *a++; k--;) {
	    c[*a++] ^= b[i];
	}
    }




}











/*

void    TSMultDmatrixBit(unsigned long m,unsigned long n,unsigned long p, unsigned long *a, unsigned long *b, unsigned long *c)
{
	unsigned long k, j, i, s,BgLine,SizeRowsb,t;
	       
        SizeRowsb=iceildiv(p,WBITS);
        
        memset(c,0,n*SizeRowsb*sizeof(unsigned long));
        
        
        BgLine=0;
        i=0;
        for (i=0;i<m; ++i)
	{k=a[BgLine];
            s = 0;
	      for (j = 1; j <= k; ++j)
		{for (t = 0; t < SizeRowsb; ++t)
			c[a[BgLine+j]*SizeRowsb+t]^=b[i*SizeRowsb+t];
		}
              BgLine+=k+1;  
         }
}

*/






/*

Sparse to Dense marix convertion

*/

void StoDBit(unsigned long m, unsigned long n, unsigned long *a,
	     unsigned long *b)
{
    unsigned long i, j, BgLine, SizeRowb;
    BgLine = 0;
    SizeRowb = iceildiv(n, WBITS);

    for (i = 0; i < m * SizeRowb; ++i) {
	b[i] = 0;
    }
    // memset(b,0,m*SizeRowb*sizeof(unsigned long));

    for (i = 0; i < m; ++i) {
	for (j = 1; j <= a[BgLine]; ++j) {

	    b[i * SizeRowb + (a[BgLine + j] / WBITS)] |=
		(1UL << (a[BgLine + j] % WBITS));

	}

	BgLine += a[BgLine] + 1;
    }
}


//////////////////////////////////////////////////
//
//  Functions to older versios of Lanczos
//  
//
///////////////////////////////////////////////////



/////////////////////////////////////////////////////
// 
//
//   Operations on Matrices (non dynamicaly created)
//
/////////////////////////////////////////////////////



/*
Add u times Row i to row j of a  m x n matrix 
*/

unsigned long AddRow(unsigned long m, unsigned long n, unsigned long a[m][n],
		     unsigned long u, unsigned long i, unsigned long j)
{
    unsigned long p;
    for (p = 0; p < n; ++p) {
	a[j][p] = (a[j][p] + a[i][p] * u) % 2;
    }
    return 0;
}


/*
Swap Rows i and j of a  m x n matrix 
*/

unsigned long SwapRows(unsigned long m, unsigned long n,
		       unsigned long a[m][n], unsigned long i,
		       unsigned long j)
{
    unsigned long k[n];
    unsigned long p;
    for (p = 0; p < n; ++p) {
	k[p] = a[i][p];
	a[i][p] = a[j][p];
	a[j][p] = k[p];
    }
    return 0;
}

/*

Multiply the i row of a  m x n matrix  by u  

*/

unsigned long MultiplyRow(unsigned long m, unsigned long n,
			  unsigned long a[m][n], unsigned long u,
			  unsigned long j)
{
    unsigned long p;
    for (p = 0; p < n; ++p) {
	a[j][p] = (a[j][p] * u) % 2;
    }
    return 0;
}




/*
void SMultDmatrixBit(unsigned long NrLines, unsigned long *a, unsigned long *b, unsigned long *c)
{
	unsigned long k, j, i, s,count,BgLine;
        i=0;
        BgLine=0;count=0;
	do {k=a[BgLine];
            s = 0;
	      for (j = 1; j <= k; ++j)
		{
			c[i]^=b[a[BgLine+j]];count++;
		}
              BgLine+=k+1;
             i++;
          } while(i<=NrLines);

}
*/









/*
Number of rows of a matrix. -1 in position 0 is the stop marker
*/

unsigned long NumbRows(unsigned long **a)
{
    unsigned long k;
    k = 0;
    while (a[k][0] != -1) {
	++k;
    }
    return k;
}


/*

Sparced matrix, dense matrix multiplication with result a dense matrix

*/

void SMultDmatrix(unsigned long m, unsigned long **a, unsigned long **b,
		  unsigned long **c)
{
    unsigned long k, j, i, s, NR;

    NR = NumbRows(a);

    //prunsigned longf("====NR=====%d\n",NR);
    for (i = 0; i < NR; ++i) {
	for (k = 0; k < m; ++k)
	{
	    s = 0;
	    for (j = 1; j <= a[i][0]; ++j) {
		if (b[a[i][j]][k] == 1) {
		    if (s == 1) {
			s = 0;
		    } else {
			s = 1;
		    }
		}
		c[i][k] = s;
	    }
	}
    }

}



/*

Convertion from Sparce to Dense matrix

*/


void StoD(unsigned long m, unsigned long n, unsigned long **a,
	  unsigned long **b)
{
    unsigned long i, j;
    for (i = 0; i < m; ++i) {
	for (j = 0; j < n; ++j) {
	    b[i][j] = 0;
	}
    }
    for (i = 0; i < m; ++i)
	for (j = 1; j <= a[i][0]; ++j)
	    b[i][a[i][j]] = 1;

}



/*

Convertion from Dense to Sparce matrix

*/


void DtoS(unsigned long m, unsigned long n, unsigned long **b,
	  unsigned long **a)
{
    unsigned long i, j;
    for (i = 0; i < m; ++i) {
	a[i][0] = 0;
	for (j = 0; j < n; ++j) {
	    if (b[i][j] == 1) {
		a[i][0] += 1;
		a[i][a[i][0]] = j;
	    }
	}
    }
    a[m][0] = -1;
}






/*

Dense Matrix sum 

*/

void DSum(unsigned long m, unsigned long n, unsigned long **a,
	  unsigned long **b, unsigned long **c)
{
    unsigned long k, j;
    for (k = 0; k < m; ++k)
	for (j = 0; j < n; ++j) {
	    if (a[k][j] == b[k][j]) {
		c[k][j] = 0;
	    } else {
		c[k][j] = 1;
	    }
	}


}


/*

Dense Matrix multiplication 

*/

void DMult(unsigned long m, unsigned long n, unsigned long p,
	   unsigned long **a, unsigned long **b, unsigned long **c)
{
    unsigned long k, j, i, s;
    for (k = 0; k < m; ++k)
	for (j = 0; j < p; ++j) {
	    s = 0;
	    for (i = 0; i < n; ++i) {
		s += (a[k][i] * b[i][j]) % 2;
	    }
	    c[k][j] = s % 2;
	}

}




/*

Dense Matrix transposition 

*/

void Transpose(unsigned long m, unsigned long n, unsigned long **a,
	       unsigned long **b)
{
    unsigned long i, j;
    for (i = 0; i < m; ++i)
	for (j = 0; j < n; ++j) {
	    b[j][i] = a[i][j];
	}
//return 0;
}


/*

Select the first m lines of the matrix a 

*/


unsigned long SelectLines(unsigned long n, unsigned long **a, unsigned long m,
			  unsigned long **b)
{
    unsigned long i, j;
    for (i = 0; i < m; ++i)
	for (j = 0; j < n; ++j)
	    b[i][j] = a[i][j];
    return 0;
}


/*

Count the exact number of coeffs in sparce matrix 

*/

unsigned long NumbCoeffSMatrix(unsigned long *a, unsigned long m)
{
    unsigned long i, Sum, Bgline;
    Bgline = 0;
    Sum = 0;
    for (i = 0; i < m; i++) {
	Sum += a[Bgline];
	Bgline += a[Bgline] + 1;
    }
    return Sum;
}
