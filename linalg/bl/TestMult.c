#include "mpi_select.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
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
    SparseMatrix M;
    DenseMatrix Ker;
    char *Fl;

    Fl = malloc(3 * sizeof(unsigned long));
    Fl = argv[1];

    unsigned long Tm1, Tm2;
    float DiffTime;



// Definition of Block size 

    unsigned long Block = WBITS;

// Test if MPI is initialized if not initializes it 

    int size, p, flag;
    flag = 0;
    MPI_Initialized(&flag);
    if (!flag) {
	MPI_Init(&argc, &argv);	//MPI initialize 
    }
// end Test if MPI is initialized if not initializes it 



    MPI_Comm_size(MPI_COMM_WORLD, &size);	//get the number of processes (size=1 if not using MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &p);	//get the process ranks (p=0 if not using MPI)



// Get the sparce matrix from file 


    unsigned long *Num;
    Num = malloc(3 * sizeof(unsigned long));
    ReadSMatrixFileData(Fl, Num);
    M->Nrows = Num[0];
    M->Ncols = Num[1];
    M->Weight = Num[2];


// For MPI  Number of Lines of the matrix to read in each process

    unsigned long BlockSize, t;
    BlockSize = SizeBlock(size, p, M->Nrows);

// end for MPI   Number of Lines of the matrix to read in each process   


    M->Data =
	malloc((M->Nrows / size + M->Nrows % size) * (M->Weight +
						    1) *
	       sizeof(unsigned long));
    ReadSMatrixFileBlockNew(Fl, M->Data, p * (M->Nrows / size), BlockSize);

    printf("Number of Coeffs = %lu \n", NumbCoeffSMatrix(M->Data, M->Nrows));

// end Get the sparce matrix from file


    unsigned long *Y, *c;
    Y = Allocmn(M->Ncols + 1, Block);
    c = Allocmn(M->Ncols + 1, Block);



    Ker->Data = Allocmn(M->Ncols, Block);
    Ker->Ncols = Block;
    Ker->Nrows = M->Ncols;

    for (t = 0; t < 10; ++t) {

//if (p==0) {Tm1=wallclock_microseconds();}

//Lanczos(M,Block,Ker);


//Testes a multiplicacao


//memset(c,0,(M->Ncols+1)*sizeof(unsigned long));

	RandomDMatrixBitTest(M->Ncols, Block, Y);

//displaySMatrixNew(M->Data,M->Nrows,M->Ncols,'a');

//displayMatrix(Y,M->Ncols,Block,'Y');

	if (p == 0) {
	    Tm1 = wallclock_microseconds();
	}

	SMultDmatrixBit(M->Nrows, M->Ncols, Block, M->Data, Y, c, 0, M->Nrows);

//displayMatrix(c,M->Ncols,Block,'c');


	if (p == 0) {
	    Tm2 = wallclock_microseconds();
	    DiffTime = Tm2 - Tm1;
	    printf("Total Time for Lanczos = %f s\n", DiffTime / 1000000);
	}

    }

    MPI_Finalize();
    return 0;
}
