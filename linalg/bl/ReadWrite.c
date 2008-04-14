#include "mpi_select.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "struct.h"
#include "ReadWrite.h"


#define	iceildiv(x,y)	(((x)+(y)-1)/(y))
#define	WBITS	(CHAR_BIT * sizeof(unsigned long))



unsigned long SizeBlock(unsigned long Nproc, unsigned long proc,
			unsigned long m)
{
    unsigned long SizeABlock;
    if (proc == Nproc - 1) {
	SizeABlock = m / Nproc + (m % Nproc);
    } else {
	SizeABlock = m / Nproc;
    }
    return SizeABlock;
}

void displayMatrixV(unsigned long **M,
		    unsigned long m, unsigned long n, char c)
{
    unsigned long i, k;
    printf("%c:=Matrix(GF(2),%lu,%lu,[", c, m, n);
    for (k = 0; k < m; ++k) {

	for (i = 0; i < n; ++i) {
	    if ((i != n - 1) || (k != m - 1))
		printf("%lu, ", M[k][i]);
	    else
		printf("%lu ", M[k][i]);
	}
    }
    printf("]);\n");

}





void displayMatrix(const unsigned long *matrix,
		   unsigned long m, unsigned long n, char c)
{
    unsigned long i;
    unsigned long width = iceildiv(n, WBITS);

    printf("%c:=Matrix(GF(2),%lu,%lu,[", c, m, n);
    for (i = 0; i < m; i++) {
	unsigned int j;
	unsigned int bit;
	for (j = 0; j < n; j++) {
	    bit = matrix[i * width + j / WBITS] >> (j % WBITS);
	    bit &= 1;
	    if (i != 0 || j != 0)
		printf(", ");
	    printf("%d", bit);
	}
    }
    printf("]);\n");

}



void displayMatrixScreen(const unsigned long *matrix,
			 unsigned long m, unsigned long n)
{
    unsigned long i;
    unsigned long width = iceildiv(n, WBITS);


    for (i = 0; i < m; i++) {
	unsigned int j;
	unsigned int bit;
	for (j = 0; j < n; j++) {
	    bit = matrix[i * width + j / WBITS] >> (j % WBITS);
	    bit &= 1;
	    printf("%d ", bit);
	}
	printf("\n");
    }
    printf("\n");

}




void WriteFileMatrix(const unsigned long *matrix,
		     unsigned long m, unsigned long n, char c)
{
    unsigned long i;
    unsigned long width = iceildiv(n, WBITS);
    FILE *File;
    File = fopen("TestesR.txt", "w");
    fprintf(File, "%c:=Matrix(GF(2),%lu,%lu,[", c, m, n);
    for (i = 0; i < m; i++) {
	unsigned int j;
	unsigned int bit;
	for (j = 0; j < n; j++) {
	    bit = matrix[i * width + j / WBITS] >> (j % WBITS);
	    bit &= 1;
	    if (i != 0 || j != 0)
		fprintf(File, ", ");
	    fprintf(File, "%d", bit);
	}
    }
    fprintf(File, "]);\n");
    fclose(File);
}




void displaySMatrix(const unsigned long *matrix,
		    unsigned long m, unsigned long n, char c)
{
    unsigned long i;
    unsigned long width = iceildiv(n, WBITS);

    printf("%c:=Matrix(GF(2),%ld,%ld,[", c, m, n);
    for (i = 0; i < m; i++) {
	unsigned int j;
	unsigned int bit;
	for (j = 0; j < n; j++) {
	    bit = matrix[i * width + j / WBITS] >> (j % WBITS);
	    bit &= 1;
	    if (i != 0 || j != 0)
		printf(", ");
	    printf("%d", bit);
	}
    }
    printf("]);\n");

}

void displaySMatrixNew(const unsigned long *matrix,
		       unsigned long m, unsigned long n, char c)
{
    unsigned long i, k;
    //unsigned long width = iceildiv(n, WBITS);

    printf("%c:=Matrix(GF(2),%ld,%ld,[", c, m, n);
    unsigned long BgLine, Line;
    BgLine = 0;
    for (Line = 0; Line < m; Line++) {

	for (i = 0; i < n; i++) {
	    for (k = BgLine + 1; k <= BgLine + matrix[BgLine]; k++) {
		if (matrix[k] > i) {
		    printf("%d", 0);
		    if (i != n - 1 || Line != m - 1) {
			printf(",");
		    };
		    break;
		}
		if (matrix[k] == i) {
		    printf("%d", 1);
		    if (i != n - 1 || Line != m - 1) {
			printf(",");
		    };
		    break;
		}
	    }
	    if (i > matrix[BgLine + matrix[BgLine]]) {
		printf("%d", 0);
		if (i != n - 1 || Line != m - 1) {
		    printf(",");
		};
	    }
	}
	//printf("BgLine= %d\n",BgLine);
	BgLine += matrix[BgLine] + 1;
	printf("\n");
    }
    printf("]);\n");
}


void WriteFileSMatrixNew(const unsigned long *matrix,
			 unsigned long m, unsigned long n, char c)
{
    unsigned long i, k;
    //unsigned long width = iceildiv(n, WBITS);
    FILE *File;
    File = fopen("TestesS.txt", "w");
    fprintf(File, "%c:=Matrix(GF(2),%ld,%ld,[", c, m, n);
    unsigned long BgLine, Line;
    BgLine = 0;
    for (Line = 0; Line < m; Line++) {

	for (i = 0; i < n; i++) {
	    for (k = BgLine + 1; k <= BgLine + matrix[BgLine]; k++) {
		if (matrix[k] > i) {
		    fprintf(File, "%d", 0);
		    if (i != n - 1 || Line != m - 1) {
			fprintf(File, ",");
		    };
		    break;
		}
		if (matrix[k] == i) {
		    fprintf(File, "%d", 1);
		    if (i != n - 1 || Line != m - 1) {
			fprintf(File, ",");
		    };
		    break;
		}
	    }
	    if (i > matrix[BgLine + matrix[BgLine]]) {
		fprintf(File, "%d", 0);
		if (i != n - 1 || Line != m - 1) {
		    fprintf(File, ",");
		};
	    }
	}
	//printf("BgLine= %d\n",BgLine);
	BgLine += matrix[BgLine] + 1;
	fprintf(File, "\n");
    }
    fprintf(File, "]);\n");
    fclose(File);
}



/*

void DisplaySMatrix(
		unsigned long * A
                )
{
unsigned long i,Test,BgLine;
BgLine=0;Test=A[BgLine];
do
{for (i=BgLine;i<=BgLine+Test; ++i)
      {
       printf("%lu ",A[i]);
       }
       printf("\n"); 
       BgLine=i;Test=A[BgLine];
} while (Test!=0);
}

*/


void DisplaySMatrix(unsigned long m, unsigned long *A)
{
    unsigned long j, i, BgLine;
    BgLine = 0;
    for (j = 0; j < m; ++j) {
	for (i = 0; i <= A[BgLine]; ++i) {
	    printf("%lu ", A[BgLine + i]);
	}
	printf("\n");
	BgLine += i;
    }
}






unsigned long SizeSMatrixFile(char *f)
{
    unsigned long i;
    char c;
    FILE *File;
    static char str[256];
    int test, Lim, count;

    File = fopen(f, "r");
    memset(str, 0, 4);
    i = 0;
    count = 0;
    test = 0;
    do {
	c = fgetc(File);
	if (c == '\n') {
	    count = 0;
	}
	if (c != ' ') {
	    strncat(str, &c, 1);
	} else {
	    if (count == 0) {
		Lim = atoi(str);

	    }
	    if (count == Lim) {
		i += Lim + 1;

		count = 1;
	    }
	    memset(str, 0, 4);
	    test = 1;
	    count += 1;
	}
    } while (c != EOF);
//printf("Sum=%d \n",i);
    fclose(File);
//free(File);
    return i;
}





void ReadSMatrixFile(char *f, unsigned long *M)
{
    int i;
    char c;
    FILE *File;
    static char str[128];

    File = fopen(f, "r");

    memset(str, 0, 4);
    i = 0;
    do {
	c = fgetc(File);

	if (c != ' ') {
	    strncat(str, &c, 1);
	} else {
	    M[i] = atoi(str);
	    i++;

	    memset(str, 0, 4);
	}
    } while (c != EOF);
//printf("Sum=%d \n",i);
    fclose(File);
//free(str);

}



void ReadSMatrixFileData(char *f, unsigned long *data)
{
    int i;
    char c;
    FILE *File;
    static char str[128];

    File = fopen(f, "r");

    memset(str, 0, 4);
    i = 0;
    do {
	c = fgetc(File);

	if (c != ' ') {
	    strncat(str, &c, 1);
	} else {
	    data[i] = atoi(str);
	    i++;
	    memset(str, 0, 4);
	}
    } while (i < 2);
//printf("Sum=%d \n",i);
    fclose(File);
//free(str);

}



void ReadSMatrixFileNew(char *f, unsigned long *M)
{
    int i;
    char c;
    FILE *File;
    static char str[128];

    File = fopen(f, "r");

    memset(str, 0, 4);
    i = 0;
    do {
	c = fgetc(File);

	if (c != ' ') {
	    strncat(str, &c, 1);
	} else {
	    if (i > 2) {
		M[i - 3] = atoi(str);
	    }
	    i++;
	    memset(str, 0, 4);
	}
    } while (c != EOF);
//printf("Sum=%d \n",i);
    fclose(File);
//free(str);

}



void ReadSMatrixFileBlock(char *f, unsigned long *M, unsigned long k,
			  unsigned long size)
{
    unsigned long l, q, i, j, i1;
    FILE *File;


    File = fopen(f, "r");


    for (j = 0; j < k; ++j) {
	fscanf(File, "%lu", &l);
	for (i = 0; i < l; ++i) {
	    fscanf(File, "%lu", &q);
	}
    }
    i = 0;
    j = 0;

    for (j = 0; j < size; ++j) {
	fscanf(File, "%lu", &l);
	M[i] = l;
	i++;
	for (i1 = 0; i1 < l; ++i1) {
	    fscanf(File, "%lu", &q);
	    M[i] = q;
	    i++;
	}
    }

    fclose(File);
//free(File);

}




void ReadSMatrixFileBlockNew(char *f, unsigned long *M, unsigned long k,
			     unsigned long size)
{
    unsigned long l, q, i, j, i1;
    FILE *File;


    File = fopen(f, "r");

    for (i = 0; i < 3; ++i) {
	fscanf(File, "%lu", &q);
    }

    for (j = 0; j < k; ++j) {
	fscanf(File, "%lu", &l);
	for (i = 0; i < l; ++i) {
	    fscanf(File, "%lu", &q);
	}
    }
    i = 0;
    j = 0;

    for (j = 0; j < size; ++j) {
	fscanf(File, "%lu", &l);
	M[i] = l;
	i++;
	for (i1 = 0; i1 < l; ++i1) {
	    fscanf(File, "%lu", &q);
	    M[i] = q;
	    i++;
	}
    }

    fclose(File);
//free(File);

}






/*

Count the number of coeff to be read for each block and broadcast the output.

*/
void CoeffperBlock( unsigned long *NumberCoeffBlocks, char *f)
{

    int p,size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);	//get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &p);	//get the process ranks

if (p==0){

    unsigned long Nrows,l, q, i, j, i1,*LengthBlocks;
    
    memset(NumberCoeffBlocks, 0, size * sizeof(unsigned long));

    LengthBlocks=malloc(size*sizeof(unsigned long));

    memset(LengthBlocks, 0, size * sizeof(unsigned long));
    FILE *File;
    File = fopen(f, "r");

    for (i = 0; i < 3; ++i) {
	fscanf(File, "%lu", &q);
        if (i==0) {Nrows=q;}
    }

    for (j = 0; j < size; ++j){if (j==0) {LengthBlocks[j]+=SizeBlock(size,j,Nrows);} else {LengthBlocks[j]+=(LengthBlocks[j-1]+SizeBlock(size,j,Nrows));} printf("Size blocks is %lu \n",LengthBlocks[j]);}  

    for (j = 0; j < size; ++j){NumberCoeffBlocks[j]+=SizeBlock(size,j,Nrows);}
   
    unsigned long Block=0;

    for (j = 0; j < Nrows; ++j) {
	fscanf(File, "%lu", &l);
        if (j<LengthBlocks[Block]){
            NumberCoeffBlocks[Block]+=l;
        } else {
        Block++;
        NumberCoeffBlocks[Block]+=(l);
	}	
	for (i1 = 0; i1 < l; ++i1) {
	    fscanf(File, "%lu", &q);
	}
    }
    fclose(File);
   }
   MPI_Bcast(NumberCoeffBlocks, size, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
}











void WriteSMatrix(unsigned long *A,
		  unsigned long m, unsigned long n, char c, char *f)
{
    unsigned long i, BgLine;
    FILE *File;
    File = fopen(f, "a");
    BgLine = 0;
    //unsigned long width = iceildiv(n, WBITS);

    fprintf(File, "%c:=Matrix(GF(2),%ld,%ld,[", c, m, n);
    for (i = 0; i < m; i++) {
	unsigned int j;
	unsigned int k;
	for (j = 0; j < n; j++) {
	    for (k = 1; k <= A[BgLine]; k++) {
		if (A[BgLine + k] == j) {
		    if (i != 0 || j != 0)
			fprintf(File, ", ");
		    fprintf(File, "%d", 1);
		    break;
		}
	    }
	    if (k > A[BgLine]) {
		if (i != 0 || j != 0)
		    fprintf(File, ", ");
		fprintf(File, "%d", 0);
	    }
	}

	BgLine += A[BgLine] + 1;

	// printf("BL=%lu\n",BgLine);
    }
    fprintf(File, "]);\n");
    fclose(File);
}






void WriteMatrix(const unsigned long *matrix,
		 unsigned long m, unsigned long n, char c, char *f)
{
    FILE *File;
    File = fopen(f, "a");
    unsigned long i;
    unsigned long width = iceildiv(n, WBITS);

    fprintf(File, "%c:=Matrix(GF(2),%ld,%ld,[", c, m, n);
    for (i = 0; i < m; i++) {
	unsigned int j;
	unsigned int bit;
	for (j = 0; j < n; j++) {
	    bit = matrix[i * width + j / WBITS] >> (j % WBITS);
	    bit &= 1;
	    if (i != 0 || j != 0)
		fprintf(File, ", ");
	    fprintf(File, "%d", bit);
	}
    }
    fprintf(File, "]);\n");
    fclose(File);

}

void WriteSMatrixSparse(unsigned long *A, unsigned long m, char *f)
{
    unsigned long i, BgLine;
    FILE *File;
    File = fopen(f, "w");
    BgLine = 0;
    //unsigned long width = iceildiv(n, WBITS);

    //fprintf(File,"%c:=Matrix(GF(2),%ld,%ld,[",c,m,n);
    for (i = 0; i < m; i++) {
	unsigned int j;

	for (j = 0; j <= A[BgLine]; j++) {
	    fprintf(File, "%lu ", A[BgLine + j]);
	}

	BgLine += A[BgLine] + 1;
	fprintf(File, "\n");
	// printf("BL=%lu\n",BgLine);
    }
    //fprintf(File,"]);\n");
    fclose(File);
}




void WriteSMatrixSparseNew(unsigned long *A,
			   unsigned long m, unsigned long n, unsigned long w,
			   char *f)
{
    unsigned long i, BgLine;
    FILE *File;
    File = fopen(f, "w");
    BgLine = 0;

//    fprintf(File, "%lu %lu %lu \n", m, n, w);

     fprintf(File, "%lu %lu \n", m, n);

    for (i = 0; i < m; i++) {
	unsigned int j;

	for (j = 0; j <= A[BgLine]; j++) {
	    fprintf(File, "%lu ", A[BgLine + j]);
	}

	BgLine += A[BgLine] + 1;
	fprintf(File, "\n");
	// printf("BL=%lu\n",BgLine);
    }
    //fprintf(File,"]);\n");
    fclose(File);
}



void WriteMatrixDense(const unsigned long *matrix,
		      unsigned long m, unsigned long n, char *f)
{
    FILE *File;
    File = fopen(f, "w");
    unsigned long i;
    unsigned long width = iceildiv(n, WBITS);

    for (i = 0; i < m; i++) {
	unsigned int j;
	unsigned int bit;
	for (j = 0; j < n; j++) {
	    bit = matrix[i * width + j / WBITS] >> (j % WBITS);
	    bit &= 1;
	    //if (i != 0 || j != 0) 
	    fprintf(File, "%d ", bit);
	}
    }

    fclose(File);

}


void ReadDMatrixFile(char *f, unsigned long n, unsigned long *M)
{
    int i, j, L;
    char c;
    FILE *File;
    static char str[128];
    L = iceildiv(n, WBITS);
    File = fopen(f, "r");
    M[0] = 0;
    memset(str, 0, 4);
    j = 0;
    i = 0;
    do {
	c = fgetc(File);
	if (c != ' ') {
	    strncat(str, &c, 1);
	} else {
	    if (atoi(str) == 1) {
		M[i] |= 1UL << (j % WBITS);
	    }
	    j++;
	    if (j % WBITS == 0) {
		i++;
		M[i] = 0;
	    }

	    memset(str, 0, 4);
	}
    } while (c != EOF);
//printf("Sum=%d \n",i);
    fclose(File);
}





void WriteBlockMatrix(DenseMatrix K, char *f)
{

int size,p;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    

if (p==0){

   

    FILE *File;
    File = fopen(f, "w");
    unsigned long i;
    for (i = 0; i < K->Nrows; i++) {
	
	fprintf(File, "%lu \n", K->Data[i]);

    }
    fclose(File);
}
    
}





/*
GetSparseMatrix(char *f,unsigned long *M)
{
unsigned long *Num;
Num=malloc(3*sizeof(unsigned long));


ReadSMatrixFileData(f,Num);



M=malloc((Num[0])*(Num[2]+1)*sizeof(unsigned long));
ReadSMatrixFileNew(f,M.Data);
return 0;
}

*/
