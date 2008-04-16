#include "mpi_select.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include "struct.h"
#include "ReadWrite.h"
#include "mat_ops.h"


#define	iceildiv(x,y)	(((x)+(y)-1)/(y))
#define	WBITS	(CHAR_BIT * sizeof(unsigned long))


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



void ReadSMatrixDimensions(char *filename, unsigned long *data)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Cannot open %s: %s\n", filename, strerror(errno));
        exit(1);
    }

    if (fscanf(f, "%lu %lu", &(data[0]), &(data[1])) != 2) {
        fprintf(stderr, "Parse error while reading matrix header\n");
        exit(1);
    }

    fclose(f);
}



#if 0
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
#endif

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



#if 0
void ReadSMatrixFileBlockNew(char *f, unsigned long *M, unsigned long k,
			     unsigned long size)
{
    unsigned long l, q, i, j, i1;
    FILE *File;


    File = fopen(f, "r");

    for (i = 0; i < 2; ++i) {
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
#endif

void ReadSMatrixSlice(SparseMatrix M, char *filename)
{
    FILE * f = fopen(filename, "r");
    unsigned long i;
    unsigned long * ptr;
    int j;
    MPI_Comm_rank(MPI_COMM_WORLD, &j);
    matrix_slice_ptr s = M->slices[j];


    if (f == NULL) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        exit(1);
    }
    if (fseeko(f, s->start, SEEK_SET) < 0) {
        fprintf(stderr, "fseek(%s,%ld): %s\n",
                filename, s->start, strerror(errno));
    }

    M->Data = malloc((s->nbcoeffs + s->i1 - s->i0) * sizeof(unsigned long));

    ptr = M->Data;

    for(i = s->i0 ; i < s->i1 ; i++) {
        unsigned long l;
        fscanf(f, "%lu", &l);
        *ptr++ = l;
        for( ; l-- ; ) {
            fscanf(f, "%lu", ptr++);
        }
    }

    assert(ptr - M->Data == (s->nbcoeffs + s->i1 - s->i0));

    fclose(f);
}


/* imported from bw/cxx/balance.c ; this reads a table of row weights. */
struct row {
    int         w;
    off_t       o;
};

static struct row * read_matrix(unsigned long * p_nr, const char * filename)
{
    char buf[1024];

    FILE * f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "%s: %s\n", filename, strerror(errno));
        exit(1);
    }

    char * ptr;
    int siz;
    off_t o;

#define READ_LINE() do {						\
        o = ftello(f);                                                  \
        ptr = fgets(buf, sizeof(buf), f);				\
        if (ptr == NULL) {						\
            fprintf(stderr, "Unexpected %s at position %ld in %s\n",    \
                    feof(f) ? "EOF" : "error", o, filename);            \
            exit(1);							\
        }								\
        siz = strlen(ptr);						\
    } while (0)

#define GOBBLE_LONG_LINE() do {                                         \
    while (!(siz < (sizeof(buf)-1) || buf[sizeof(buf)-2] == '\n')) {	\
        ptr = fgets(buf, sizeof(buf), f);				\
        if (ptr == NULL) {						\
            fprintf(stderr, "Unexpected %s at position %ld in %s\n",    \
                    feof(f) ? "EOF" : "error", ftello(f), filename);    \
            exit(1);							\
        }								\
        siz = strlen(ptr);						\
    }                                                                   \
    } while (0)								\

    
    uint32_t nr, nc;
    uint32_t i;

    /* Read matrix header */
    READ_LINE();
    if (sscanf(buf, "%" SCNu32 " %" SCNu32, &nr, &nc) != 2) {
        fprintf(stderr, "Parse error while reading header: %s\n", buf);
        exit(1);
    }
    GOBBLE_LONG_LINE();

    struct row * row_table = malloc((nr+1) * sizeof(struct row));
    if (row_table == NULL) abort();

    fprintf(stderr, "Reading row table...");
    fflush(stderr);

    for(i = 0 ; i < nr ; i++) {
        int w;
        READ_LINE();
        if (sscanf(buf, "%d", &w) != 1) {
            fprintf(stderr, "Parse error while reading line %"PRIu32" of %s\n",
                    i+1, filename);
            exit(1);
        }
        row_table[i].w = w;
        row_table[i].o = o;
        GOBBLE_LONG_LINE();
    }
    row_table[nr].o = ftello(f);

    fprintf(stderr, "done\n");
    fflush(stderr);

    fclose(f);

    *p_nr = nr;

    return row_table;
}





/* Count the number of coeff to be read for each block and broadcast the
 * output.
 */
void PrepareMatrixSlices(SparseMatrix M, char *filename)
{
    int p, nb_processes;
    int i;
    unsigned long * offsets;
    unsigned long * weights;
    unsigned long * indices;

    MPI_Comm_size(MPI_COMM_WORLD, &nb_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);

    offsets = malloc((nb_processes+1) * sizeof(unsigned long));
    weights = malloc( nb_processes * sizeof(unsigned long));
    indices = malloc((nb_processes+1) * sizeof(unsigned long));


    /* Only one thread reads the matrix */
    if (p == 0) {
        unsigned long nrows;
        struct row * row_table = read_matrix(&nrows, filename);
        int j;
	for (j = 0; j < nb_processes + 1; ++j) {
            unsigned long i0 = (j * nrows) / nb_processes;
            indices[j] = i0;
            offsets[j] = row_table[i0].o;
            /* Count all the coeffs of the previous slice */
            if (j) {
                unsigned long i;
                weights[j-1] = 0;
                for(i = indices[j-1] ; i < i0 ; i++) {
                    weights[j-1] += row_table[i].w;
                }
            }
        }
    }
    MPI_Bcast(offsets, nb_processes + 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(weights, nb_processes, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(indices, nb_processes, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    /* reconstruct the nice data */

    M->slices = malloc(nb_processes*sizeof(matrix_slice));
    for(i = 0 ; i < nb_processes ; i++) {
        M->slices[i]->i0 = indices[i];
        M->slices[i]->i1 = indices[i+1];
        M->slices[i]->start = offsets[i];
        M->slices[i]->end = offsets[i+1];
        M->slices[i]->nbcoeffs = weights[i];
    }

    /* Dispose the scrap paper */
    free(offsets);
    free(weights);
    free(indices);
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

    int size, p;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);


    if (p == 0) {
	FILE *File;
	File = fopen(f, "w");
	unsigned long i;
	for (i = 0; i < K->Nrows; i++) {
	    fprintf(File, "%016lx \n", K->Data[i]);
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
