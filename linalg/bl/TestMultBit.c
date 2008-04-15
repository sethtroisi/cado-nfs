#include <stdint.h>
#include "mpi_select.h"
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
    SparseMatrix M;
    DenseMatrix Ker;
    
    char *Fl;

    Fl = argv[1];


    // unsigned long Tm1, Tm2;
    // float DiffTime;
    unsigned long Block = WBITS;

// Test if MPI is initialized if not initializes it 

    int size, p, flag;
    flag = 0;
    MPI_Initialized(&flag);
    if (!flag) {
	MPI_Init(&argc, &argv);	//MPI initialize 
    }
// end Test if MPI is initialized if not initializes it 


    init_random();

    MPI_Comm_size(MPI_COMM_WORLD, &size);	//get the number of processes (size=1 if not using MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &p);	//get the process ranks (p=0 if not using MPI)
    // MPI_Status status;


// Get the sparse matrix from file 


    unsigned long *Num;
    Num = malloc(2 * sizeof(unsigned long));
    ReadSMatrixFileData(Fl, Num);
    M->Nrows = Num[0];
    M->Ncols = Num[1];
    
   free(Num);

    // For MPI  Number of Lines of the matrix to read in each process
    unsigned long BlockSize,i;
    BlockSize = SizeBlock(size, p, M->Nrows);

    // end for MPI   Number of Lines of the matrix to read in each process   
    
    unsigned long *NumberCoeffBlocks=malloc(size*sizeof(unsigned long));

    CoeffperBlock(NumberCoeffBlocks,Fl);


    if (p==0) {
        for (i=0; i<size; i++) {
            printf("Job %d Block  %lu  has size %lu \n",
                    p, i, NumberCoeffBlocks[i]);
        }
    }

   M->Data =malloc(NumberCoeffBlocks[p]*sizeof(unsigned long));


    ReadSMatrixFileBlockNew(Fl, M->Data, p * (M->Nrows / size),BlockSize);


// end Get the sparce matrix from file

    DenseMatrix Result;
    
    Result->Data = Allocmn(M->Nrows, Block);
    Ker->Data = Allocmn(M->Ncols, Block);


//for tests

RandomDMatrixBitTest(M->Ncols,Block, Ker->Data);


Ker->Nrows=M->Ncols;
Ker->Ncols=Block;

unsigned long *KerPar;
KerPar=malloc(2*sizeof(unsigned long));
KerPar[0]=Ker->Nrows;
KerPar[1]=Ker->Ncols;

// unsigned long m=M->Nrows;
// unsigned long n=M->Ncols;

//printf("p = %lu   m = %lu   n = %lu!!\n",p,m,n);


//printf("p=%lu  m  %lu   n %lu \n",p,m,n);

 

//SMatrix_Vector(Result,M,Ker);


//Test_SMatrix_Vector(Result->Data, M->Data, Ker->Data,m,n);

printf("Multiplication done in %lu!!\n",p);


//free(NumberCoeffBlocks);
//free(KerPar);
//free(Result->Data);
//free(Ker->Data);
//free(M->Data);


/*
  unsigned long *ResmN_Dist1;
    ResmN_Dist1 = Allocmn(SizeBlock(size, p, m), Block);


printf("p=  %lu Toto\n",p);

    unsigned long j, *d;
    d = malloc(m * sizeof(unsigned long));


    MPI_Bcast(Ker->Data, n, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);


    //MPI_Bcast(NC, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);


    SMultDmatrixBit(m, n, Block, M->Data, Ker->Data, ResmN_Dist1, p * (m / size),
		    SizeBlock(size, p, m));


    if (p != 0) {
	MPI_Send(ResmN_Dist1, SizeBlock(size, p, m), MPI_UNSIGNED_LONG, 0, 1,
		 MPI_COMM_WORLD);
    };

    if (p == 0) {

	for (i = 0; i < SizeBlock(size, 0, m); i++) {
	    d[i] = ResmN_Dist1[i];
	}

	unsigned long Bg = SizeBlock(size, p, m);

	for (i = 1; i < size; i++) {


	    MPI_Recv(ResmN_Dist1, SizeBlock(size, i, m), MPI_UNSIGNED_LONG, i,
		     1, MPI_COMM_WORLD, &status);

	    for (j = 0; j < SizeBlock(size, i, m); j++) {
		d[j + Bg] = ResmN_Dist1[j];
	    }

	    Bg += SizeBlock(size, i, m);

	}
    }

*/
    

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

//      Tm1 = clock();

//TransposeBit(m, Block, Y, e);
//      DMultBit(m, Block, Block, Y, d, c1);

//      VUBit(m, Block, Y, d, c);



//      Tm2 = clock();
//      DiffTime = Tm2 - Tm1;
//      printf("Total Time for TVU = %f s\n",
//             DiffTime / 1.0 / CLOCKS_PER_SEC);


//    }

//    displayMatrix(c, m, Block, 'c');
//    displayMatrix(c1, m, Block, 'e');

//displayMatrix(d,Block,Block,'d');
