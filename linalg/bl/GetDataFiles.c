
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "basic.h"
#include "random.h"
#include "alloc.h"
#include "mat_ops.h"
#include "echelon.h"
#include "ReadWrite.h"



#define	iceildiv(x,y)	(((x)+(y)-1)/(y))
#define	WBITS	(CHAR_BIT * sizeof(unsigned long))

unsigned long *SMultDmatrixBitTest(unsigned long m, unsigned long n,
				   unsigned long p, unsigned long *a,
				   unsigned long *b)
{
    unsigned long k, j, i, s, BgLine, SizeRowsb, t;
    unsigned long *c;
    c = Allocmn(m, p);

    SizeRowsb = iceildiv(p, WBITS);

    memset(c, 0, m * SizeRowsb * sizeof(unsigned long));


    BgLine = 0;
    i = 0;
    for (i = 0; i < m; ++i) {
	k = a[BgLine];
	s = 0;
	for (j = 1; j <= k; ++j) {
	    for (t = 0; t < SizeRowsb; ++t)
		c[i * SizeRowsb + t] ^= b[a[BgLine + j] * SizeRowsb + t];
	}
	BgLine += k + 1;
    }

    return c;
}



int main(int argc, char *argv[])
{

    unsigned long m;
    unsigned long n;
    unsigned long Block;
    unsigned long k, dev;

    if (argc == 1) {
	m = 999;
	n = 1000;
	Block = 64;
	k = 100;
	dev = 10;
    } else {
	if (argc == 3) {
	    m = atoi(argv[1]);
	    n = atoi(argv[2]);
	    Block = WBITS;
	    k = n / 1000;
	    dev = 10;
	} else if (argc == 4) {
	    m = atoi(argv[1]);
	    n = atoi(argv[2]);
	    Block = atoi(argv[3]);
	    k = n / 1000;
	    dev = 10;
	} else if (argc == 5) {
	    m = atoi(argv[1]);
	    n = atoi(argv[2]);
	    Block = atoi(argv[3]);
	    k = atoi(argv[4]);
	    dev = 10;
	} else if (argc == 6) {
	    m = atoi(argv[1]);
	    n = atoi(argv[2]);
	    Block = atoi(argv[3]);
	    k = atoi(argv[4]);
	    dev = atoi(argv[5]);
	}

    }




//      m=100000;
//      n =100000;                
    //      Block=64;
    //      k=100;
    //      dev=10;


    printf("m= %lu   n= %lu  Block= %lu   k= %lu   dev= %lu\n", m, n, Block,
	   k, dev);
    unsigned long *A;
    unsigned long *C;
    char *f, *f1, *f2;
    f = "./Block.txt";
    f1 = "./Sparse.txt";
    f2 = "./Amatrix.txt";
    C = Allocmn(n, Block);
    //A1=Allocmn(n,Block);  

    A = malloc(n * (k + dev + 1) * sizeof(unsigned long));

    printf("Size of %lu times %lu Matrix C= %lu Bytes\n", n, Block,
	   n * iceildiv(Block, WBITS) * sizeof(unsigned long));
    //   RandomDMatrixBitTest(n,Block,C);
    ReadSMatrixFileBlock(f1, A, 0, n);
    WriteSMatrix(A, n, n, 'A', f2);
    //DisplaySMatrix(n,A);   
    //  displayMatrix(C,n,Block,'c'); 
    ReadDMatrixFile(f, Block, C);
    displayMatrix(C, n, Block, 'C');
    printf("Size of %lu times %lu Matrix A=%lu Bytes\n", n, n,
	   n * (k + dev + 1) * sizeof(unsigned long));

    //   CreateRandomSymmetricMatrix(n,k,A,dev);

    //   WriteSMatrixSparse(A,n,f1);

    //   WriteMatrixDense(C,n,Block,f);


    free(C);
    free(A);
    close_random();
    printf("Done\n");

    return 0;
}

	       // printf("toto\n");
	     //   DisplaySMatrix(n,A);

		 //printf("toto\n");
	     //   DisplaySMatrix(A);
	     //   printf("toto\n");
	     //   A1=malloc((10)*(k+dev+1)*sizeof(unsigned long));
	     //   ReadSMatrixFileBlock(f1,A1,5,10);
	     //   DisplaySMatrix(10,A1);
	      // ReadDMatrixFile(f,Block,C);
	     //    displayMatrix(C,n,Block,'C');
	     //    A1=SMultDmatrixBitTest(n,n,Block,A,C);         
	    //     displayMatrix(A1,n,Block,'A');
	     //    A2=malloc((n/3)*(k+dev+1)*sizeof(unsigned long));        
	     //    ReadSMatrixFileBlock(f1,A2,0,n/3);
	     //    A3=SMultDmatrixBitTest(n/3,n,Block,A2,C); 
	     //    displayMatrix(A3,n/3,Block,'a');              
	     //    ReadSMatrixFileBlock(f1,A2,n/3,n/3);
	     //    A3=SMultDmatrixBitTest(n/3,n,Block,A2,C);
	     //    displayMatrix(A3,n/3,Block,'b');
	     //    ReadSMatrixFileBlock(f1,A2,2*n/3,n/3+n%3);
	     //    A3=SMultDmatrixBitTest(n/3+n%3,n,Block,A2,C);
	     //                     displayMatrix(A3,n/3+n%3,Block,'c');    
	      //Lanczos(m,n,Block,A1,C,X1);
