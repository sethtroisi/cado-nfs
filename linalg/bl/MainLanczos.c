#include <stdint.h>
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

    Fl = argv[1];


    unsigned long Tm1, Tm2;
    float DiffTime;
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
    //M->Weight = Num[2];

    free(Num);

    //printf("%s \n",Fl2);

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

 //   M->Data =
//	malloc((M->Nrows / size + M->Nrows % size) * (M->Weight + 1) *
//	       sizeof(unsigned long));

	M->Data =
	malloc(NumberCoeffBlocks[p]*sizeof(unsigned long));


    ReadSMatrixFileBlockNew(Fl, M->Data, p * (M->Nrows / size),BlockSize);

// end Get the sparce matrix from file

    DenseMatrix Result;

    Result->Data = Allocmn(M->Nrows, Block);
    Ker->Data = Allocmn(M->Ncols, Block);


    if (p == 0) {
	Tm1 = microseconds();
    }

    Lanczos(Ker, M, Block);

    if (p == 0) {
	Tm2 = microseconds();
	DiffTime = Tm2 - Tm1;
	printf
	    ("Total Time for Lanczos = %f s   Size of block in kernel is %lu\n",
	     DiffTime / 1000000, Ker->Ncols);
    }


    SMatrix_Vector(Result,M,Ker);

//    Test_SMatrix_Vector(Result->Data, M->Data, Ker->Data,M->Nrows,M->Ncols);

    free(M->Data);

    if (p == 0) {
	if (TestZero(M->Nrows, Block, Result->Data)) {
	    printf
		("There were NO errors During the process. The vectors are in the kernel of the matrix  \n");
	} else {
	    printf("There was an error in the process \n");
	}

    }

    if ((p == 0) && (argc == 3)) {
	char *Fl2 = argv[2];
	WriteBlockMatrix(Ker, Fl2);
	printf("Kernel written to file %s\n", Fl2);
    }

    free(Result->Data);
    free(Ker->Data);
    free(NumberCoeffBlocks);

    close_random();

    return 0;
}






//    printf("P= %lu   Ker %lu\n",p,Ker->Nrows);
    

//if (p==0) {displayMatrix(Ker->Data,Ker->Nrows,Block,'K');}

//if (p==0) {displaySMatrixNew(M->Data,M->Nrows,M->Ncols,'M');}

