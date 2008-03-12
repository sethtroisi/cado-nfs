#include <stdint.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include "struct.h"
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

    unsigned long Tm1, Tm2;
    float DiffTime;



// Definition of Block size 

    unsigned long Block = WBITS;


    unsigned long *e, t, *d1, *c1;
    uint64_t *Y, *c, *d;
    Y = Allocmn(m + 1, Block);
    c = Allocmn(m + 1, Block);
    c1 = Allocmn(m + 1, Block);
    e = Allocmn(Block, m + 1);
    d = Allocmn(Block, Block);
    d1 = Allocmn(Block, Block);

    RandomDMatrixBitTest(m, Block, Y);
    RandomDMatrixBitTest(Block, Block, d);

//displayMatrix(Y,m,Block,'Y');
//displayMatrix(c,m,Block,'c');


    for (t = 0; t < 1; ++t) {

//if (p==0) {Tm1=wallclock_microseconds();}

//Lanczos(M,Block,Ker);


//Testes a multiplicacao


//memset(c,0,(M.Ncols+1)*sizeof(unsigned long));



//displaySMatrixNew(M.Data,M.Nrows,M.Ncols,'a');

//displayMatrix(Y,m,Block,'Y');

	Tm1 = clock();

//TransposeBit(m, Block, Y, e);
	DMultBit(m, Block, Block, Y, d, c1);

	VUBit(m, Block, Y, d, c);



	Tm2 = clock();
	DiffTime = Tm2 - Tm1;
	printf("Total Time for TVU = %f s\n",
	       DiffTime / 1.0 / CLOCKS_PER_SEC);


    }

    displayMatrix(c, m, Block, 'c');
    displayMatrix(c1, m, Block, 'e');

//displayMatrix(d,Block,Block,'d');


    return 0;
}
