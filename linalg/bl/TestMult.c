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
    fprintf(stderr, "Usage: TestMult <matrix>\n");
    exit(1);
}


int main(int argc, char *argv[])
{
    SparseMatrix M;
    // DenseMatrix Ker;
    char *matrix_file_name = NULL;

    argc--,argv++;
    for( ; argc ; argc--,argv++) {
        if (matrix_file_name == NULL) {
            matrix_file_name = argv[0];
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

    int nb_processes, p, flag;
    flag = 0;
    MPI_Initialized(&flag);
    if (!flag) {
	MPI_Init(&argc, &argv);
    }

    init_random();

    MPI_Comm_size(MPI_COMM_WORLD, &nb_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);


// Get the sparse matrix from file 

     {
         unsigned long Num[2];
         ReadSMatrixDimensions(matrix_file_name, Num);
         M->Nrows = Num[0];
         M->Ncols = Num[1];
     }

    unsigned long i;
    int t;
    
    PrepareMatrixSlices(M, matrix_file_name);

    if (p == 0) {
        for (i = 0; i < nb_processes; i++) {
            printf("Job %d Block %lu has size %lu coeffs,"
                    " %lu rows ; i0=%lu, i1=%lu\n",
                    p, i, M->slices[i]->nbcoeffs,
                    M->slices[i]->i1 - M->slices[i]->i0,
                    M->slices[i]->i0,
                    M->slices[i]->i1);
        }
    }


    ReadSMatrixSlice(M, matrix_file_name);


    unsigned long *Y, *c;
    Y = Allocmn(M->Ncols + 1, Block);
    c = Allocmn(M->Ncols + 1, Block);



    // Ker->Data = Allocmn(M->Ncols, Block);
    // Ker->Ncols = Block;
    // Ker->Nrows = M->Ncols;

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

        matrix_slice_ptr me = M->slices[p];
        SMultDmatrixBit(M->Nrows, M->Ncols, Block, M->Data, Y, c,
            me->i0, me->i1 - me->i0);

//displayMatrix(c,M->Ncols,Block,'c');


	if (p == 0) {
	    Tm2 = wallclock_microseconds();
	    DiffTime = Tm2 - Tm1;
	    printf("matrix multiplication = %.2f s\n", DiffTime / 1000000);
	}

    }

    MPI_Finalize();
    return 0;
}
