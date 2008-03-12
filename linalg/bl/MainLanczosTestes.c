#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "struct.h"
#include "basic.h"
#include "random.h"
#include "alloc.h"
#include "mat_ops.h"
#include "echelon.h"
#include "ReadWrite.h"
#include "timing.h"
#include "ProcLanczosTestes.h"



#define	iceildiv(x,y)	(((x)+(y)-1)/(y))
#define	WBITS	(CHAR_BIT * sizeof(unsigned long))



int main(int argc, char *argv[])
{
    struct SparseMatrix M;
    struct DenseMatrix Ker;
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
    M.Nrows = Num[0];
    M.Ncols = Num[1];
    M.Weight = Num[2];


// For MPI  Number of Lines of the matrix to read in each process

    unsigned long BlockSize;
    BlockSize = SizeBlock(size, p, M.Nrows);

// end for MPI   Number of Lines of the matrix to read in each process   


    M.Data =
	malloc((M.Nrows / size + M.Nrows % size) * (M.Weight +
						    1) *
	       sizeof(unsigned long));
    ReadSMatrixFileBlockNew(Fl, M.Data, p * (M.Nrows / size), BlockSize);

// end Get the sparce matrix from file


    Ker.Data = Allocmn(M.Ncols, Block);
    Ker.Ncols = Block;
    Ker.Nrows = M.Ncols;

    if (p == 0) {
	Tm1 = wallclock_microseconds();
    }



    Lanczos(M, Block, Ker);



    if (p == 0) {
	Tm2 = wallclock_microseconds();
	DiffTime = Tm2 - Tm1;
	printf("Total Time for Lanczos = %f s\n", DiffTime / 1000000);

	char *f;
	f = malloc(10 * sizeof(unsigned long));
	sprintf(f, "ControlFile_%lu_%lu_%lu.txt", M.Nrows, M.Ncols, Block);
	FILE *File;
	File = fopen(f, "a");
	fprintf(File, "Total Time for Lanczos = %f s\n", DiffTime / 1000000);
	fclose(File);

    }

    MPI_Finalize();
    return 0;
}
