#include <stdint.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <time.h>
//#include "struct.h"
#include "basic.h"
#include "random.h"
#include "alloc.h"
#include "mat_ops.h"
#include "echelon.h"
#include "ReadWrite.h"
#include "timing.h"
#include "ProcLanczos.h"



#define	iceildiv(x,y)	(((x)+(y)-1)/(y))
#define	WBITS	(CHAR_BIT * sizeof(unsigned long))


int main(int argc, char *argv[])
{

    unsigned long m = atoi(argv[1]);;

    unsigned long T3,T1, T2,i;
    float DiffTime2,DiffTime1;



// Definition of Block size 

    unsigned long Block = WBITS;


    unsigned long *d1, *c1,*c2;
    // unsigned long *e;
    uint64_t *Y, *c, *d;
    Y = Allocmn(Block,m + 1);
    c = Allocmn(m+1,Block);
    c2 = Allocmn(m+1,Block);
    c1 = Allocmn( Block,m + 1);
   // e = Allocmn(Block, m + 1);
    d = Allocmn(Block,Block);
    d1 = Allocmn(Block, Block);

    RandomDMatrixBitTest(Block, m, Y);
    RandomDMatrixBitTest(Block, Block, d);

DiffTime1=0;
DiffTime2=0;

unsigned long **V;
    V = Allocmn3(m, Block);

for (i = 0; i < m * iceildiv(Block, WBITS); ++i) {
	V[0][i] = Y[i];
	V[1][i] = Y[i];
	V[2][i] = Y[i];
    }


for (i=0; i<1; i++) {

T1=microseconds();

    VUBit(m, Block, Y, d, c);

    //DMultBit(m,Block, 1, Y, d, c2);

  // TransposeBit(Block,Block,d,d1);

T2=microseconds();

   d1=Allocmn(Block,Block);
   TransposeBit(Block,Block,d,d1);

   //   Refresh_ArrayBit(m, Block, V, Y);

 //   TVUBit(m, Block, Y, Y, c);
  //   VUBit_v2(m, Block, Y, d, c2);

T3=microseconds();

DiffTime1+=T2-T1;
DiffTime2+=T3-T2;


/*
displayMatrix(Y,m,Block,'Y');
displayMatrix(d,Block,Block,'d');

displayMatrix(c,m,Block,'c');
*/
//displayMatrix(c2,m,Block,'C');
/*
unsigned long *r;r=Allocmn(3,Block);r[0]=10;r[1]=99;r[2]=0;


for (i = 0; i <WBITS; i++) {
            r[2]^=((__builtin_parity(r[0]^d[i]) & 1UL) << i);
	    displayMatrixScreen(&r[2],1,Block);

	}

displayMatrix(&r[0],1,Block,'v');
displayMatrix(d,Block,Block,'B');
displayMatrix(&r[2],1,Block,'p');


//displayMatrixScreen(r,3,Block);

printf("Paridade de %lu ^ %lu=%lu",r[0],r[1],__builtin_parity(r[0]^r[1]));
*/
//printf("-----------------------------------------------------------------------------------------------------\n");

//displayMatrixScreen(c2,m,Block);

}

printf("Time for alg1 VxBolck is = %f s Time for alg2 TVxV is = %f s\n",DiffTime1/1000000,DiffTime2/1000000 );

return 0;
}



//     displayMatrixScreen(c,Block,m);


//    for (t = 0; t < m; ++t) {


 
//printf("Position of Pivot at column %lu  is  %lu s\n",t,Pivot(Block,m,Y,2,t));


//if (p==0) {Tm1=wallclock_microseconds();}

//Lanczos(M,Block,Ker);


//Testes a multiplicacao


//memset(c,0,(M.Ncols+1)*sizeof(unsigned long));



//displaySMatrixNew(M.Data,M.Nrows,M.Ncols,'a');

//displayMatrix(Y,m,Block,'Y');

//	Tm1 = clock();

//TransposeBit(m, Block, Y, e);
//	DMultBit(m, Block, Block, Y, d, c1);

//	VUBit(m, Block, Y, d, c);



//	Tm2 = clock();
//	DiffTime = Tm2 - Tm1;
//	printf("Total Time for TVU = %f s\n",
//	       DiffTime / 1.0 / CLOCKS_PER_SEC);


//    }

//    displayMatrix(c, m, Block, 'c');
//    displayMatrix(c1, m, Block, 'e');

//displayMatrix(d,Block,Block,'d');

