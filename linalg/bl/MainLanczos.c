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

void usage()
{
    fprintf(stderr, "Usage: MainLanczos [options] <matrix> [<kernel>]\n"
            "Accepted options:\n"
            "--column-dependency\tCompute a column dependency (default)\n"
            "--row-dependency\tCompute a row dependency\n");
    exit(1);
}

int main(int argc, char *argv[])
{
    SparseMatrix M;
    DenseMatrix Ker;
    char *matrix_file_name = NULL;
    char *kernel_file_name = NULL;

    argc--,argv++;
    for( ; argc ; argc--,argv++) {
        if (strcmp(argv[0], "--column-dependency") == 0) {
            /* do nothing, column dep. is the default */
        } else if (strcmp(argv[0], "--row-dependency") == 0) {
            // transpose_mat = 1;
        } else if (matrix_file_name == NULL) {
            matrix_file_name = argv[0];
        } else if (kernel_file_name == NULL) {
            kernel_file_name = argv[0];
        } else {
            fprintf(stderr, "Unexpected command line argument: %s\n", argv[0]);
            usage();
        }
    }
    if (matrix_file_name == NULL) {
        usage();
    }







    unsigned long Tm1, Tm2;
    float DiffTime;
    unsigned long Block = WBITS;

// Test if MPI is initialized if not initializes it 

    int nb_processes, p, flag;
    flag = 0;
    MPI_Initialized(&flag);
    if (!flag) {
	MPI_Init(&argc, &argv);	//MPI initialize 
    }
// end Test if MPI is initialized if not initializes it 


    init_random();

    MPI_Comm_size(MPI_COMM_WORLD, &nb_processes);	//get the number of processes (nb_processes=1 if not using MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &p);	//get the process ranks (p=0 if not using MPI)
    // MPI_Status status;


// Get the sparse matrix from file 


     {
         unsigned long Num[2];
         ReadSMatrixDimensions(matrix_file_name, Num);
         M->Nrows = Num[0];
         M->Ncols = Num[1];
     }

    unsigned long i;
    
    PrepareMatrixSlices(M, matrix_file_name);

    if (p == 0) {
        for (i = 0; i < nb_processes; i++) {
            printf("Job %d Block %lu has size %lu coeffs, %lu rows\n", p, i,
                    M->slices[i]->nbcoeffs, M->slices[i]->i1 - M->slices[i]->i0);
        }
    }

    ReadSMatrixSlice(M, matrix_file_name);

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
	    ("Total Time for Lanczos = %f s   nb_processes of block in kernel is %lu\n",
	     DiffTime / 1000000, Ker->Ncols);
    }


    SMatrix_Vector(Result,M,Ker);

    if (p == 0) {
	if (TestZero(M->Nrows, Block, Result->Data)) {
	    printf
		("There were NO errors During the process. The vectors are in the kernel of the matrix  \n");
	} else {
	    printf("There was an error in the process \n");
	}

    }

    if ((p == 0) && kernel_file_name) {
	WriteBlockMatrix(Ker, kernel_file_name);
	printf("Kernel written to file %s\n", kernel_file_name);
    }

    free(Result->Data);
    free(Ker->Data);
    free(M->Data);
    free(M->slices);


    close_random();

    return 0;
}






//    printf("P= %lu   Ker %lu\n",p,Ker->Nrows);
    

//if (p==0) {displayMatrix(Ker->Data,Ker->Nrows,Block,'K');}

//if (p==0) {displaySMatrixNew(M->Data,M->Nrows,M->Ncols,'M');}

