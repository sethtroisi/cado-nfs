/*
 * Program: merge
 * Author : F. Morain
 * Purpose: merging relations
 * 
 * Algorithm: digest and interpolation from Cavallar.
 *
 */

// TODO: use unified compact lists...!
// TODO: reintroduce mat->weight

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "utils/utils.h"
#include "files.h"
#include "gzip.h"
#include "sparse.h"
#include "merge.h"
#include "prune.h"

#define DEBUG 0
#define TEX 0

#define TRACE_COL -1 // 253224 // 231 // put to -1 if not...!
#define TRACE_ROW -1 // put to -1 if not...!

#define USE_CONNECT 0

#define MERGE_LEVEL_MAX 256 // maximum level for a merge; such a large value
                            // is only useful when not using BW

#ifdef USE_MPI
unsigned int mpi_index = 0;
#endif

#define M_STRATEGY 3 // 0: finish that mergelevel
                     // 1: change if min weight < mergelevel
                     // 2: jump to minimal possible mergelevel
                     // 3: perform one merge, then check for next min weight

// 0 means dummy filtering using slow methods
// 1 means:
//  * fast with optimized-though-memory-consuming data structure;
//  * heavy columns are killed.
// 2 means:
//  * fast as above;
//  * heavy columns are present, but not in S[]; this yields more accurate
//    weights.
// 3 means:
//  * reintroduce columns whose weight is <= mergelevelmax (not suggested
//    for large numbers).
#define USE_MERGE_FAST 2

typedef struct {
  int len;    /* number of non-zero entries */
  INT *val;   /* array of size len containing column indices of non-zero
                 entries */
} rel_t;

static int 
cmp(const void *p, const void *q) {
    int x = *((int *)p);
    int y = *((int *)q);
    return (x <= y ? -1 : 1);
}

dclist
dclistCreate(INT j)
{
    dclist dcl = (dclist)malloc(sizeof(struct dclist));

    dcl->j = j;
    dcl->prev = NULL;
    dcl->next = NULL;
    return dcl;
}

void
dclistPrint(FILE *file, dclist dcl)
{
    while(dcl != NULL){
      fprintf(file, " -> %ld", (long int) dcl->j);
	dcl = dcl->next;
    }
}

int
dclistLength(dclist dcl)
{
    int l = 0;

    while(dcl != NULL){
	l++;
	dcl = dcl->next;
    }
    return l;
}

void
dclistTex(FILE *file, dclist dcl)
{
    fprintf(stderr, "\\begin{array}{r}\n");
    while(dcl != NULL){
      fprintf(file, "%ld\\\\\n", (long int) dcl->j);
	dcl = dcl->next;
    }
    fprintf(stderr, "\\end{array}\n");
}

/* insert j in doubly-chained list dcl (between cell of dcl and dcl->next),
   and return pointer at cell containing j */
static dclist
dclistInsert(dclist dcl, INT j)
{
    dclist newdcl = dclistCreate(j);

    newdcl->next = dcl->next;
    newdcl->prev = dcl;
    if(newdcl->next != NULL)
	newdcl->next->prev = newdcl;
    dcl->next = newdcl;
    return newdcl;
}

// connect dcl to S
void
dclistConnect(dclist S, dclist dcl)
{
    dcl->next = S->next;
    dcl->prev = S;
    S->next = dcl;
}

// TODO: ncols could be a new renumbered ncols...
void
initMat(sparse_mat_t *mat, int cwmax, int rwmax, int mergelevelmax)
{
#if USE_MERGE_FAST
    // cwmax+2 to prevent bangs
    dclist *S = (dclist *)malloc((cwmax+2) * sizeof(dclist));
    dclist *A = (dclist *)malloc(mat->ncols * sizeof(dclist));
    INT **R;
    int j;
#endif
#if USE_TAB == 0
    mat->data = (rel_t*) malloc (mat->nrows * sizeof (rel_t));
    ASSERT_ALWAYS(mat->data != NULL);
#else
    mat->rows = (INT **)malloc(mat->nrows * sizeof (INT *));
    ASSERT_ALWAYS(mat->rows != NULL);
#endif
    mat->wt = (int*) malloc (mat->ncols * sizeof (int));
    memset(mat->wt, 0, mat->ncols * sizeof (int));
    mat->rwmax = rwmax;
    mat->cwmax = cwmax;
    mat->mergelevelmax = mergelevelmax;
#if USE_MERGE_FAST
    // S[0] has a meaning at least for temporary reasons
    for(j = 0; j <= cwmax+1; j++)
	S[j] = dclistCreate(-1);
    // TODO_MPI: that R could well be reduced to R[jmin..jmax[
    R = (INT **)malloc(mat->ncols * sizeof(INT *));
    ASSERT_ALWAYS(R != NULL);
    mat->S = S;
    mat->A = A;
    mat->R = R;
#endif
}

/* Initializes the SWAR data structure.
   Inputs:
      mat->ncols - number of columns of matrix mat
      mat->wt[j] - weight of column j, for 0 <= j < mat->ncols
      mat->cwmax - weight bound
      mat->S[w]  - lists containing only one element (-1), 0 <= w <= mat->cwmax
   Outputs:
      mat->S[w]  - doubly-chained list with columns j of weight w
 */
void
fillSWAR(sparse_mat_t *mat, INT jmin, INT jmax)
{
    INT j, *Rj;

    for(j = jmin; j < jmax; j++){
#  if DEBUG >= 1
	fprintf(stderr, "Treating column %d\n", j);
#  endif
	if(mat->wt[j] <= mat->cwmax){
	    mat->A[j] = dclistInsert (mat->S[mat->wt[j]], j);
#  if DEBUG >= 1
	    fprintf(stderr, "Inserting %d in S[%d]:", j, mat->wt[j]);
	    dclistPrint(stderr, mat->S[mat->wt[j]]->next);
	    fprintf(stderr, "\n");
#  endif
	    Rj = (INT *)malloc((mat->wt[j]+1) * sizeof(INT));
	    Rj[0] = 0; // last index used
	    mat->R[j] = Rj;
	}
	else{
#if USE_MERGE_FAST <= 1
	    mat->wt[j] = -1;
#else
	    mat->wt[j] = -mat->wt[j]; // trick!!!
#endif
	    mat->A[j] = NULL; // TODO: renumber j's?????
	    mat->R[j] = NULL;
	}
    }
}

// TODO
void
closeSWAR(/*sparse_mat_t *mat*/)
{
}

void
printStats(sparse_mat_t *mat)
{
    int w;

    for(w = 0; w <= mat->cwmax; w++)
	fprintf(stderr, "I found %d primes of weight %d\n",
		dclistLength(mat->S[w]->next), w);
}

void
checkData(sparse_mat_t *mat)
{
    int j, nbj = 0;

    for(j = 0; j < mat->ncols; j++)
	if(mat->wt[j])
	    nbj++;
    if(mat->rem_ncols != nbj){
	fprintf(stderr, "rem_ncols=%d nbj=%d\n", mat->rem_ncols, nbj);
	exit(1);
    }
}

// dump for debugging reasons
void
printSWAR(sparse_mat_t *mat, int ncols)
{
    int j, w;
    INT k;

    fprintf(stderr, "===== S is\n");
    for(w = 0; w <= mat->cwmax; w++){
	fprintf(stderr, "  S[%d] ", w);
	dclistPrint(stderr, mat->S[w]->next);
	fprintf(stderr, "\n");
    }
    fprintf(stderr, "===== Wj is\n");
    for(j = 0; j < ncols; j++)
	fprintf(stderr, "  Wj[%d]=%d\n", j, mat->wt[j]);
    fprintf(stderr, "===== R is\n");
    for(j = 0; j < ncols; j++){
	if(mat->R[j] == NULL)
	    continue;
	fprintf(stderr, "  R[%d]=", j);
	for(k = 1; k <= mat->R[j][0]; k++)
          fprintf(stderr, " %ld", (long int) mat->R[j][k]);
	fprintf(stderr, "\n");
    }
}

void
matrix2tex(sparse_mat_t *mat)
{
    int i, j, k, *tab;

    tab = (int *)malloc(mat->ncols * sizeof(int));
    fprintf(stderr, "$$\\begin{array}{l|");
    for(j = 0; j < mat->ncols; j++)
	fprintf(stderr, "c");
    fprintf(stderr, "}\n");
    for(i = 0; i < mat->nrows; i++){
	memset(tab, 0, mat->ncols * sizeof(int));
	for(k = 0; k < lengthRow(mat, i); k++)
	    tab[cell(mat, i, k)] = 1;
	fprintf(stderr, "R_{%d}", i);
	for(j = 0; j < mat->ncols; j++)
	    fprintf(stderr, "& %d", tab[j]);
	fprintf(stderr, "\\\\\n");
    }
    fprintf(stderr, "\\hline\n");
    fprintf(stderr, "\\text{weight}");
    for(j = 0; j < mat->ncols; j++)
	fprintf(stderr, "& %d", mat->wt[j]);
    fprintf(stderr, "\\\\\n");
    fprintf(stderr, "\\end{array}$$\n");
    free(tab);
}

void
texSWAR(sparse_mat_t *mat)
{
    int j, w;
    INT k;

    fprintf(stderr, "\\end{verbatim}\n");
    matrix2tex(mat);
    fprintf(stderr, "$$\\begin{array}{|");
    for(w = 0; w <= mat->cwmax; w++)
	fprintf(stderr, "r|");
    fprintf(stderr, "|r|}\\hline\n");
    for(w = 0; w <= mat->cwmax+1; w++){
	fprintf(stderr, "S[%d]", w);
	if(w < mat->cwmax+1)
	    fprintf(stderr, "&");
    }
    fprintf(stderr, "\\\\\\hline\n");
    for(w = 0; w <= mat->cwmax+1; w++){
	dclistTex(stderr, mat->S[w]->next);
	if(w < mat->cwmax+1)
            fprintf(stderr, "&");
    }
    fprintf(stderr, "\\\\\\hline\n");
    fprintf(stderr, "\\end{array}$$\n");
    fprintf(stderr, "$$\\begin{array}{|c|");
    for(j = 0; j < mat->ncols; j++)
	fprintf(stderr, "r|");
    fprintf(stderr, "}\\hline\n");
    fprintf(stderr, "j");
    for(j = 0; j < mat->ncols; j++)
	fprintf(stderr, "& %d", j);
    fprintf(stderr, "\\\\\\hline\n");
    fprintf(stderr, "W[j]");
    for(j = 0; j < mat->ncols; j++)
	fprintf(stderr, "& %d", mat->wt[j]);
    fprintf(stderr, "\\\\\\hline\n\\end{array}$$\n");
    fprintf(stderr, "\\begin{verbatim}\n");
    for(j = 0; j < mat->ncols; j++){
	if(mat->R[j] == NULL)
	    continue;
	fprintf(stderr, "  R[%d]=", j);
	for(k = 1; k <= mat->R[j][0]; k++)
          fprintf(stderr, " %ld", (long int) mat->R[j][k]);
	fprintf(stderr, "\n");
    }
}

// terrific hack: everybody on the same line
// the first is to be destroyed in replay!!!
void
reportn(FILE *outfile, INT *ind, int n)
{
    int i;

#if DEBUG >= 1
    fprintf(stderr, "Reporting for n=%d\n", n);
#endif
#ifdef USE_MPI
    fprintf(outfile, "REPORT %d ", mpi_index);
#endif
    for(i = 0; i < n; i++){
	fprintf(outfile, "%ld", (long int) ind[i]);
	if(i < n-1)
	    fprintf(outfile, " ");
    }
    fprintf(outfile, "\n");
}

void
report1(FILE *outfile, INT i)
{
    reportn(outfile, &i, 1);
}

void
report2(FILE *outfile, INT i1, INT i2)
{
    INT tmp[2];
    tmp[0] = i1;
    tmp[1] = i2;
    reportn(outfile, tmp, 2);
}

/* compute the weight mat->wt[j] of each column j, for the set of relations
   in file purgedfile.
*/
void
initWeightFromFile(sparse_mat_t *mat, FILE *purgedfile, int jmin, int jmax)
{
    int i, j, ret, x, nc;

    memset(mat->wt, 0, mat->ncols * sizeof(int));
    for(i = 0; i < mat->nrows; i++){
	ret = fscanf(purgedfile, "%d", &j); // unused index to rels file
	ASSERT_ALWAYS (ret == 1);
	ret = fscanf (purgedfile, "%d", &nc);
	ASSERT_ALWAYS (ret == 1);
	for(j = 0; j < nc; j++){
	    ret = fscanf(purgedfile, PURGE_INT_FORMAT, &x);
	    ASSERT_ALWAYS (ret == 1);
	    if((x >= jmin) && (x < jmax))
		mat->wt[x]++;
	}
    }
}

#define BUF_LEN 100

/* Reads a matrix file, and puts in mat->wt[j], 0 <= j < ncols, the
   weight of column j (adapted from matsort.c).

   We skip columns that are too heavy.

   mat->wt is already filled and put to - values when needed. So don't touch.
*/
void
readmat(sparse_mat_t *mat, FILE *file, int jmin, int jmax)
{
    int ret;
    int i, j, k;
    int nc, buf[BUF_LEN];
    INT ibuf;

    ret = fscanf (file, "%d %d", &(mat->nrows), &(mat->ncols));
    ASSERT_ALWAYS (ret == 2);
    
    fprintf(stderr, "Reading matrix of %d rows and %d columns: excess is %d\n",
	    mat->nrows, mat->ncols, mat->nrows - mat->ncols);
    mat->rem_nrows = mat->nrows;
    mat->weight = 0;

    for (i = 0; i < mat->nrows; i++){
	if(!(i % 100000))
	    fprintf(stderr, "Reading %d-th row\n", i);
	ret = fscanf(file, "%d", &j); // unused index to rels file
	ASSERT_ALWAYS (ret == 1);
	ret = fscanf (file, "%d", &nc);
	ASSERT_ALWAYS (ret == 1);
	if(nc == 0){
#if USE_TAB == 0
	    lengthRow(mat, i) = 0;
	    mat->data[i].val = NULL;
#else
	    mat->rows[i] = NULL;
#endif
	    mat->rem_nrows--;
	}
	else{
	    for(k = 0, ibuf = 0; k < nc; k++){
		ret = fscanf(file, PURGE_INT_FORMAT, &j);
		ASSERT_ALWAYS (ret == 1);
		ASSERT_ALWAYS (0 <= j && j < mat->ncols);
#if USE_MERGE_FAST <= 1
		if(mat->wt[j] > 0)
		    buf[ibuf++] = j;
#else
		// always store j in the right interval...!
		if((j >= jmin) && (j < jmax)){
		    buf[ibuf++] = j;
		    // this will be the weight in the current slice
		    mat->weight++; 
		    if(mat->wt[j] > 0){ // redundant test?
			mat->R[j][0]++;
			mat->R[j][mat->R[j][0]] = i;
		    }
#if DEBUG >= 1
		    if(j == 15054){
			fprintf(stderr, "i=%d, j=%d (ind/nc=%d) mat->%d\n", 
				i, j, k, mat->R[j][0]);
			fflush(stderr);
		    }
#endif
		}
#endif
	    }
	    ASSERT_ALWAYS(ibuf <= BUF_LEN);
#if USE_TAB == 0
	    lengthRow(mat, i) = ibuf;
	    mat->data[i].val = (int*) malloc (ibuf * sizeof (int));
	    memcpy(mat->data[i].val, buf, ibuf * sizeof(int));
	    // sort indices in val to ease row merges
	    qsort(mat->data[i].val, ibuf, sizeof(int), cmp);
#else
	    mat->rows[i] = (INT*) malloc ((ibuf+1) * sizeof (INT));
	    mat->rows[i][0] = ibuf;
	    memcpy(mat->rows[i]+1, buf, ibuf * sizeof(INT));
	    // sort indices in val to ease row merges
            qsort(mat->rows[i]+1, ibuf, sizeof(INT), cmp);
#endif
	}
    }
    // we need to keep informed of what really happens; this will be an upper
    // bound on the number of active columns, I guess
    mat->rem_ncols = mat->ncols;
    printStats(mat);
}

void
print_row(sparse_mat_t *mat, int i)
{
#if USE_TAB == 0
    int k;
    
    for(k = 0; k < lengthRow(mat, i); k++)
	fprintf(stderr, " %d", cell(mat, i, k));
#else
    fprintRow(stderr, mat->rows[i]);
#endif
}

void
destroyRow(sparse_mat_t *mat, int i)
{
#if USE_TAB == 0
    free(mat->data[i].val);
    mat->data[i].val = NULL;
    lengthRow(mat, i) = 0;
#else
    free(mat->rows[i]);
    mat->rows[i] = NULL;
#endif
}

void
removeWeightFromRow(sparse_mat_t *mat, int i)
{
#if USE_TAB == 0
    int k;

    for(k = 0; k < lengthRow(mat, i); k++)
	mat->wt[cell(mat, i, k)]--;
#else
    removeWeight(mat->rows, mat->wt, i);
#endif
#if TRACE_ROW >= 0
    if(i == TRACE_ROW)
	fprintf(stderr, "TRACE_ROW: removeWeightFromRow %d\n", 
		lengthRow(mat, i));
#endif
    mat->weight -= lengthRow(mat, i);
}

void
addWeightFromRow(sparse_mat_t *mat, int i)
{
#if USE_TAB == 0
    int k;

    for(k = 0; k < lengthRow(mat, i); k++)
	mat->wt[cell(mat, i, k)]++;
#else
    addWeight(mat->rows, mat->wt, i);
#endif
    mat->weight += lengthRow(mat, i);
}

int
hasCol(sparse_mat_t *mat, int i, INT j)
{
    INT k;

#if USE_TAB == 0
    for(k = 0; k < lengthRow(mat, i); k++)
	if(cell(mat, i, k) == j)
	    return 1;
#else
    for(k = 1; k <= mat->rows[i][0]; k++)
	if(cell(mat, i, k) == j)
	    return 1;
#endif
    return 0;
}

// i1 += i2; weights are taken care of outside...
void
addRowsData(rel_t *data, int i1, int i2)
{
    int k1, k2, k, len;
    INT *tmp, *tmp2;

#if DEBUG >= 2
    fprintf(stderr, "row[%d] =", i1); print_row(mat, i1);
    fprintf(stderr, "\n");
    fprintf(stderr, "row[%d] =", i2); print_row(mat, i2);
    fprintf(stderr, "\n");
#endif
    // merge row[i1] and row[i2] in i1...
    len = data[i1].len + data[i2].len;
    tmp = (INT *)malloc(len * sizeof(INT));
    k = k1 = k2 = 0;

    // loop while everybody is here
    while((k1 < data[i1].len) && (k2 < data[i2].len)){
	if(data[i1].val[k1] < data[i2].val[k2])
	    tmp[k++] = data[i1].val[k1++];
	else if(data[i1].val[k1] > data[i2].val[k2])
            tmp[k++] = data[i2].val[k2++];
	else{
#if DEBUG >= 1
	    fprintf(stderr, "ADD[%d=(%ld,%lu), %d=(%ld,%lu)]: new w[%d]=%lu\n", 
		    i1, data[i1].a, data[i1].b,
		    i2, data[i2].a, data[i2].b,
		    data[i1].val[k1], 
		    wt[data[i1].val[k1]]);
#endif
	    k1++;
	    k2++;
	}
    }
    // finish with k1
    for( ; k1 < data[i1].len; k1++)
	tmp[k++] = data[i1].val[k1];
    // finish with k2
    for( ; k2 < data[i2].len; k2++)
	tmp[k++] = data[i2].val[k2];
    // destroy and copy back
    free(data[i1].val);
    tmp2 = (INT *)malloc(k * sizeof(INT));
    memcpy(tmp2, tmp, k * sizeof(INT));
    data[i1].len = k;
    data[i1].val = tmp2;
    free(tmp);
#if DEBUG >= 2
    fprintf(stderr, "row[%d]+row[%d] =", i1, i2);
    print_row(mat, i1); fprintf(stderr, "\n");
#endif
}

// i1 += i2, mat->wt is updated at the same time.
void
addRowsWithWeight(sparse_mat_t *mat, int i1, int i2)
{
    // i1 is to disappear, replaced by a new one
    removeWeightFromRow(mat, i1);
#if USE_TAB == 0
    addRowsData(mat->data, i1, i2);
#else
    addRows(mat->rows, i1, i2, -1);
#endif
    // new row i1 has to contribute to the weight
    addWeightFromRow(mat, i1);
#if TRACE_ROW >= 0
    if(i1 == TRACE_ROW)
	fprintf(stderr, "TRACE_ROW: addRowsWithWeight i1=%d\n", i1);
    if(i2 == TRACE_ROW)
	fprintf(stderr, "TRACE_ROW: addRowsWithWeight i2=%d\n", i2);
#endif
}

// what is the weight of the sum of Ra and Rb?
int
weightSum(sparse_mat_t *mat, int i1, int i2)
{
    int k1, k2, w = 0;

#if USE_TAB == 0
    k1 = k2 = 0;
    while((k1 < lengthRow(mat, i1)) && (k2 < lengthRow(mat, i2))){
#else
    k1 = k2 = 1;
    while((k1 <= lengthRow(mat, i1)) && (k2 <= lengthRow(mat, i2))){
#endif
	if(cell(mat, i1, k1) < cell(mat, i2, k2)){
	    k1++;
	    w++;
	}
	else if(cell(mat, i1, k1) > cell(mat, i2, k2)){
	    k2++;
	    w++;
	}
	else{
	    k1++; k2++;
	}
    }
    w += (k1 > lengthRow(mat, i1) ? 0 : lengthRow(mat, i1)-k1+1);
    w += (k2 > lengthRow(mat, i2) ? 0 : lengthRow(mat, i2)-k2+1);
    return w;
}

int
findAllRowsWithGivenj(INT *ind, sparse_mat_t *mat, INT j, int nb)
{
    int i, k, r = 0;

    // TODO: special hack for nb==2???
    for(i = 0; i < mat->nrows; i++){
	if(isRowNull(mat, i))
	    continue;
#if USE_TAB == 0
	for(k = 0; k < lengthRow(mat, i)-1; k++) // trick!
	    if(cell(mat, i, k) >= j)
		break;
	if(cell(mat, i, k) == j){
	    ind[r++] = i;
	    if(r == nb)
		return 1;
	}
#else
	for(k = 1; k <= lengthRow(mat, i)-1; k++) // trick!
	    if(cell(mat, i, k) >= j)
		break;
	if(cell(mat, i, k) == j){
	    ind[r++] = i;
	    if(r == nb)
		return 1;
	}
#endif
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////
// Prim

#define QUEUE_TYPE 1 // 0 for naive, 1 for heap

#if QUEUE_TYPE == 0

// we put fQ in Q[0][0], lQ in Q[0][1].

#define isQueueEmpty(Q) (Q[0][0] == Q[0][1])

void
popQueue(int *s, int *t, int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3])
{
    int fQ = Q[0][0];

    *s = Q[fQ][0];
    *t = Q[fQ][1];
    Q[0][0]++;
}

void
printQueue(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3])
{
    int i;

    for(i = Q[0][0]; i < Q[0][1]; i++)
	fprintf(stderr, "Q[%d] = [%d, %d, %d]\n", i,Q[i][0], Q[i][1], Q[i][2]);
}

void
addEdge(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3], int u, int v, int Auv)
{
    int i, j, fQ = Q[0][0], lQ = Q[0][1];

    for(i = fQ; i < lQ; i++)
	if(Auv < Q[i][2])
	    break;
    // shift everybody: Auv >= Q[i][2]
    for(j = lQ; j > i; j--){
	Q[j][0] = Q[j-1][0];
	Q[j][1] = Q[j-1][1];
	Q[j][2] = Q[j-1][2];
    }
    Q[i][0] = u;
    Q[i][1] = v;
    Q[i][2] = Auv;
    Q[0][1]++;
}

#else
// Q[0][0] contains the number of items in Q[], so that useful part of Q
// is Q[1..Q[0][0]].

#define isQueueEmpty(Q) (Q[0][0] == 0)

void
printQueue(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3])
{
    int level = 0, imax = 1, i;

    fprintf(stderr, "L0:");
    for(i = 1; i <= Q[0][0]; i++){
	fprintf(stderr, " %d", Q[i][2]);
	if(i == imax){
	    imax = (imax<<1)+1;
	    fprintf(stderr, "\nL%d:", ++level);
	}
    }
    fprintf(stderr, "\n");
}

void
upQueue(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3], int k)
{
    int x = Q[k][0], y = Q[k][1], v = Q[k][2];

    while((k > 1) && (Q[k/2][2] >= v)){
	// we are at level > 0 and the father is >= son
	// the father replaces the son
	Q[k][0] = Q[k/2][0];
	Q[k][1] = Q[k/2][1];
	Q[k][2] = Q[k/2][2];
	k /= 2;
    }
    // we found the place of (x, y, v)
    Q[k][0] = x;
    Q[k][1] = y;
    Q[k][2] = v;
}

void
addEdge(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3], int u, int v, int Auv)
{
    Q[0][0]++;
    Q[Q[0][0]][0] = u;
    Q[Q[0][0]][1] = v;
    Q[Q[0][0]][2] = Auv;
    upQueue(Q, Q[0][0]);
    if(Q[0][0] >= MERGE_LEVEL_MAX*MERGE_LEVEL_MAX-2)
	fprintf(stderr, "size(Q) -> %d\n", Q[0][0]);
}

// Move Q[1] down.
void
downQueue(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3], int k)
{
    int x = Q[k][0], y = Q[k][1], v = Q[k][2], j;

    while(k <= Q[0][0]/2){
	// k has at least a left son
	j = 2*k;
	if(j < Q[0][0])
	    // k has a right son
	    if(Q[j][2] > Q[j+1][2])
		j++;
	// at this point, Q[j] is the largest son
	if(v <= Q[j][2])
	    break;
	else{
	    // the father takes the place of the son
	    Q[k][0] = Q[j][0];
	    Q[k][1] = Q[j][1];
	    Q[k][2] = Q[j][2];
	    k = j;
	}
    }
    // we found the place of v
    Q[k][0] = x;
    Q[k][1] = y;
    Q[k][2] = v;
}

// smallest edge (s, t) is always in Q[1]
void
popQueue(int *s, int *t, int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3])
{
    *s = Q[1][0];
    *t = Q[1][1];
    Q[1][0] = Q[Q[0][0]][0];
    Q[1][1] = Q[Q[0][0]][1];
    Q[1][2] = Q[Q[0][0]][2];
    Q[0][0]--;
    downQueue(Q, 1);
}
#endif

// Add all neighbors of u, which is already in T; all edges with v in T
// are not useful anymore.
void
addAllEdges(int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3],
	    int u, int *father, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
    int v;

    for(v = 0; v < m; v++)
	if((v != u) && (father[v] < 0))
	    addEdge(Q, u, v, A[u][v]);
}

// height[0] = 0, etc. Returns hmax.
int
minimalSpanningTreeWithPrim(int *father, int *height, 
			    int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
			    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
    int Q[MERGE_LEVEL_MAX*MERGE_LEVEL_MAX][3]; // over-conservative
    int u, s, t, i, nV, hmax = 0;

    // nodes of T
    for(i = 0; i < m; i++){
	father[i] = -1;
	height[i] = 0;
	sons[i][0] = 0; // last occupied position for a son of i
    }
    u = 0;
    father[u] = u; // for the root, isn't it?
    nV = m-1;
#if QUEUE_TYPE == 0
    Q[0][0] = 1;
    Q[0][1] = 1;
    addAllEdges(Q, u, father, A, m);
    // active part of Q is Q[fQ..lQ[
    ASSERT(Q[0][0] == 1);
    ASSERT(Q[0][1] == m);
#else
    Q[0][0] = 0;
    addAllEdges(Q, u, father, A, m);
#endif
#if DEBUG >= 1
    printQueue(Q);
#endif
    while(! isQueueEmpty(Q)){
	// while queue is non empty
	// pop queue
        popQueue(&s, &t, Q);
#if DEBUG >= 1
	fprintf(stderr, "Popping a = (%d, %d) of weight %d\n", s, t, A[s][t]);
#endif
	if(father[t] == -1){
	    // t does not close a cycle
	    // T[u] <- T[u] union (s, t)
#if DEBUG >= 1
	    fprintf(stderr, "new edge: (%d, %d) of weight %d\n",s,t,A[s][t]);
#endif
	    father[t] = s;
	    // store new son for s
	    sons[s][0]++;
	    sons[s][sons[s][0]] = t;
	    height[t] = height[s]+1;
	    if(height[t] > hmax)
		hmax = height[t];
	    nV--;
	    if(nV == 0)
		break;
	    addAllEdges(Q, t, father, A, m);
#if DEBUG >= 1
	    printQueue(Q);
#endif
	}
    }
    return hmax;
}

#if 0
void
minimalSpanningTreeWithKruskal(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
#if 0
    int *E, m2 = m*(m-1), i, j, iE = 0, lE = 3*m2;
    int nV = 0, *ancester, *father;

    E = (int *)malloc(lE * sizeof(int));
    for(i = 0; i < m; i++)
	for(j = 0; j < m; j++){
	    if(j != i){
		E[iE++] = A[i][j];
		E[iE++] = i;
		E[iE++] = j;
	    }
	}
    // sort edges
    qsort(E, m2, 3 * sizeof(int), cmp);
#if DEBUG >= 1
    for(i = 0; i < m2; i++)
	fprintf(stderr, "w(%d, %d)=%d\n", E[3*i+1], E[3*i+2], E[3*i]);
#endif
    ancester = (int *)malloc(m * sizeof(int));
    father = (int *)malloc(m * sizeof(int));
    for(i = 0; i < m; i++){
	ancester[i] = -1;
	father[i] = -1; // everybody is a root
    }
    nV = 0;
    iE = -3;
    while(1){
	iE += 3;
	if((ancester[E[iE+1]] == -1) && (ancester[E[iE+2]] == -1)){
	    // new cc
	    ancester[E[iE+1]] = E[iE+1];
	    ancester[E[iE+2]] = father[E[iE+2]] = E[iE+1];
	    nV += 2;
	}
	else if(ancester[E[iE+1]] == -1){
	    // new edge to cc of E[iE+2]
	    father[E[iE+1]] = E[iE+2];
	    ancester[E[iE+1]] = ancester[E[iE+2]];
	    nV++;
	}
	else if(ancester[E[iE+2]] == -1){
            // new edge to cc of E[iE+1]
            father[E[iE+2]] = E[iE+1];
            ancester[E[iE+2]] = ancester[E[iE+1]];
	    nV++;
        }
	else{
	    // E[iE+1] and E[iE+2] are already known, check ancesters
	    int anc1 = findAncester(E[iE+1], ancester);
	    int anc2 = findAncester(E[iE+2], ancester);
	    // if anc1 == anc2, then this would create a cycle
	    if(anc1 != anc2){
		// we have to merge the two cc's
		ancester[E[iE+2]] = ancester[E[iE+1]];
		// but how do we deal with the father's...?
	    }
	}
	if(nV == m)
	    break;
    }
    free(E);
    free(ancester);
    free(father);
#endif
}
#endif

int
minimalSpanningTree(int *father, int *height, 
		    int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		    int m, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
    return minimalSpanningTreeWithPrim(father, height, sons, A, m);
    //    minimalSpanningTreeWithKruskal(A, m);
}

void
fillRowAddMatrix(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], sparse_mat_t *mat,
                 int m, INT *ind)
{
    int i, j;

    for(i = 0; i < m; i++)
	A[i][i] = 0;
    // A[i][j] <- Weight(R[ind[i]]+R[ind[j]]);
    for(i = 0; i < m; i++)
	for(j = i+1; j < m; j++){
	    A[i][j] = weightSum(mat, ind[i], ind[j]);
	    A[j][i] = A[i][j];
#if DEBUG >= 1
	    fprintf(stderr, "A[%d][%d] = A[%d][%d] = %d\n",i,j,j,i,A[i][j]);
#endif
	}
}

// not mallocing to speed up(?).
int
findBestIndex(sparse_mat_t *mat, int m, INT *ind)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], i, j, imin, wmin, w;

    ASSERT(m <= MERGE_LEVEL_MAX);
    if(m == 2)
	return 0;
    fillRowAddMatrix(A, mat, m, ind);
    // iterate over all vertices
    imin = -1;
    wmin = 0;
    for(i = 0; i < m; i++){
	// compute the new total weight if i is used as pivot
	w = 0;
	for(j = 0; j < m; j++)
	    // note A[i][i] = 0
	    w += A[i][j];
#if DEBUG >= 1
	fprintf(stderr, "W[%d]=%d\n", i, w);
#endif
	if((imin == -1) || (w < wmin)){
	    imin = i;
	    wmin = w;
	}
    }
    return imin;
}

void
merge_m_slow(FILE *outfile, sparse_mat_t *mat, int m, int verbose)
{
    int totfindrows = 0, totfindbest = 0, totadd = 0, tt;
    int j, k, i, njproc = 0;
    INT *ind;
    int report = (verbose == 0) ? 10000 : 1000;

#if DEBUG >= 1
    fprintf(stderr, "Weight %d:", m);
#endif
    ind = (INT *)malloc(m * sizeof(INT));
    for(j = 0; j < mat->ncols; j++){
	if(mat->wt[j] != m)
	    continue;
	njproc++;
	if(!(njproc % report))
	    fprintf(stderr, "# %d columns of weight %d processed\n",njproc,m);
	// we need to find the m rows and then
	tt = cputime();
	if(!findAllRowsWithGivenj(ind, mat, j, m))
	    fprintf(stderr, "Could not find the %d required rows\n", m);
	totfindrows += (cputime()-tt);
#if DEBUG >= 1
	fprintf(stderr, " %d", j);
	fprintf(stderr, "=> the %d rows are:\n", j, m);
	for(k = 0; k < m; k++){
	    fprintf(stderr, "row[%d]=", ind[k]);
	    print_row(mat, ind[k]);
	    fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
#endif
	// try all combinations to find the smaller one
	tt = cputime();
	i = findBestIndex(mat, m, ind);
	totfindbest += (cputime()-tt);
	if(i == -1){
	    fprintf(stderr, "Sorry, could not find best index\n");
	    // or do some more clever stuff?
	    continue;
	}
#if DEBUG >= 1
	fprintf(stderr, "Minimal is i=%d\n", i);
#endif
	tt = cputime();
	for(k = 0; k < m; k++)
	    if(k != i){
		addRowsWithWeight(mat, ind[k], ind[i]);
#if DEBUG >= 1
		fprintf(stderr, "new row[%d]=", ind[k]);
		print_row(mat, ind[k]);
		fprintf(stderr, "\n");
#endif
	    }
	if(i > 0){
	    // put ind[i] in front
	    int itmp = ind[i];
            for(k = i; k > 0; k--)
                ind[k] = ind[k-1];
	    ind[0] = itmp;
	}
	reportn(outfile, ind, m);
	totadd += (cputime()-tt);
	removeWeightFromRow(mat, ind[0]);
	destroyRow(mat, ind[0]);
	mat->rem_nrows--;
	mat->rem_ncols--;
	mat->wt[j] = 0;
    }
    fprintf(stderr, "TIME: findrows=%d findbest=%d add=%d\n",
	    totfindrows, totfindbest, totadd);
    free(ind);
}

// We should have mat->wt[j] == m == nchk
void
checkCoherence(sparse_mat_t *mat, int m, int j)
{
    int nchk = 0, k;

    for(k = 1; k <= mat->R[j][0]; k++)
	if(mat->R[j][k] != -1)
	    nchk++;
    ASSERT(nchk == (mat->wt[j] >= 0 ? mat->wt[j] : -mat->wt[j]));
    if(m != -1){
	if(nchk != m){
	    fprintf(stderr, "HYPERCHECK:");
	    fprintf(stderr, "mat->R[%d][0]=%ld, m=%d\n", j,
                    (long int) mat->R[j][0], m);
	    fprintf(stderr, "Gasp: nchk=%d\n", nchk);
	}
	ASSERT(nchk == m);
    }
}

// w(j) has decreased in such a way that it can be incorporated in the
// SWAR structure. We need find all the info we need to do that. Note this
// can be costly. We must not store row index i0, since j was discovered
// when row i0 was treated.
void
incorporateColumn(sparse_mat_t *mat, INT j, int i0)
{
    int i, ni, wj;
    INT *Rj;

    wj = -mat->wt[j];
    Rj = (INT *)malloc((wj+1) * sizeof(INT));
    // find all rows in which j appears
    for(i = 0, ni = 1; i < mat->nrows; i++)
	if(!isRowNull(mat, i))
	    if((i != i0) && hasCol(mat, i, j))
#if 1
		Rj[ni++] = i;
#else
                ni++;
    fprintf(stderr, "iC: %d %d\n", j, ni);
#endif
    Rj[0] = ni-1;
    mat->R[j] = Rj;
    mat->wt[j] = wj;
    ASSERT(wj == Rj[0]);
    mat->A[j] = dclistInsert (mat->S[wj], j);
}

int
getNextj(dclist dcl)
{
    INT j;
    dclist foo;

    foo = dcl->next;
    j = foo->j;
    dcl->next = foo->next;
    if(foo->next != NULL)
	foo->next->prev = dcl;
    free(foo);
    return j;
}

void
remove_j_from_S(sparse_mat_t *mat, int j)
{
    dclist dcl = mat->A[j], foo;

    if(dcl == NULL){
	fprintf(stderr, "Column %d already removed?\n", j);
	return;
    }
#if DEBUG >= 2
    int ind = mat->wt[j];
    if(ind > mat->cwmax)
        ind = mat->cwmax+1;
    fprintf(stderr, "S[%d]_b=", ind);
    dclistPrint(stderr, mat->S[ind]->next); fprintf(stderr, "\n");
    fprintf(stderr, "dcl="); dclistPrint(stderr, dcl); fprintf(stderr, "\n");
#endif
    foo = dcl->prev;
    foo->next = dcl->next;
    if(dcl->next != NULL)
	dcl->next->prev = foo;
#if DEBUG >= 2
    fprintf(stderr, "S[%d]_a=", ind);
    dclistPrint(stderr, mat->S[ind]->next); fprintf(stderr, "\n");
#endif
#if USE_CONNECT == 0
    free(dcl);
#endif
}

void
destroyRj(sparse_mat_t *mat, int j)
{
    free(mat->R[j]);
    mat->R[j] = NULL;
}

void
remove_j_from_SWAR(sparse_mat_t *mat, int j)
{
    remove_j_from_S(mat, j);
#if USE_CONNECT
    free(mat->A[j]);
#endif
    mat->A[j] = NULL;
    mat->wt[j] = 0;
    destroyRj(mat, j);
}

// Don't touch to R[j][0]!!!!
void
remove_i_from_Rj(sparse_mat_t *mat, int i, int j)
{
    // be dumb for a while
    int k;

    if(mat->R[j] == NULL){
	fprintf(stderr, "Row %d already empty\n", j);
	return;
    }
    for(k = 1; k <= mat->R[j][0]; k++)
	if(mat->R[j][k] == i){
#if DEBUG >= 1
	    fprintf(stderr, "Removing row %d from R[%d]\n", i, j);
#endif
	    mat->R[j][k] = -1;
	    break;
	}
}

void
add_i_to_Rj(sparse_mat_t *mat, int i, int j)
{
    // be semi-dumb for a while
    int k;
    
#if DEBUG >= 1
    fprintf(stderr, "Adding row %d to R[%d]\n", i, j);
#endif
    for(k = 1; k <= mat->R[j][0]; k++)
	if(mat->R[j][k] == -1)
	    break;
    if(k <= mat->R[j][0]){
	// we have found a place where it is -1
	mat->R[j][k] = i;
    }
    else{
#if DEBUG >= 1
	fprintf(stderr, "WARNING: reallocing things in add_i_to_Rj for R[%d]\n", j);
#endif
	int l = mat->R[j][0]+2;
	mat->R[j] = (INT *)realloc(mat->R[j], l * sizeof(INT));
	mat->R[j][l-1] = i;
	mat->R[j][0] = l-1;
    }
}

int
decrS(int w)
{
#if USE_MERGE_FAST <= 1
    return w-1;
#else
    return (w >= 0 ? w-1 : w+1);
#endif
}

int
incrS(int w)
{
#if USE_MERGE_FAST <= 1
    return w+1;
#else
    return (w >= 0 ? w+1 : w-1);
#endif
}

/* remove the cell (i,j), and updates matrix correspondingly.
   Note: A[j] contains the address of the cell in S[w] where j is stored.
   
   Updates:
   - mat->wt[j] (weight of column j)
   - mat->S[w] : cell j is removed, with w = weight(j)
   - mat->S[w-1] : cell j is added
   - A[j] : points to S[w-1] instead of S[w]
*/
void
removeCellSWAR(sparse_mat_t *mat, int i, INT j)
{
    int ind;

#if TRACE_ROW >= 0
    if(i == TRACE_ROW){
	fprintf(stderr, "TRACE_ROW: removeCellSWAR i=%d j=%d\n", i, j);
    }
#endif
    // update weight
#if DEBUG >= 1
    fprintf(stderr, "removeCellSWAR: moving j=%d from S[%d] to S[%d]\n",
	    j, mat->wt[j], decrS(mat->wt[j]));
#endif
#if USE_MERGE_FAST > 1
    if(mat->wt[j] < 0){
	// if mat->wt[j] is already < 0, we don't care about
	// decreasing, updating, etc. except when > 2
# if USE_MERGE_FAST > 2
	ind = mat->wt[j] = decrS(mat->wt[j]);
	// we incorporate the column and update the data structure
	if(abs(ind) <= mat->mergelevelmax){
#  if DEBUG >= 1
	    fprintf(stderr, "WARNING: column %d becomes light at %d...!\n",
		    j, abs(ind));
#  endif
	    incorporateColumn(mat, j, i);
	}
# endif
	return;
    }
#endif
    // at this point, we should have mat->wt[j] > 0
    ind = mat->wt[j] = decrS(mat->wt[j]);
    remove_j_from_S(mat, j);
    if(mat->wt[j] > mat->cwmax)
	ind = mat->cwmax+1;
    // update A[j]
#if DEBUG >= 2
    fprintf(stderr, "S[%d]_b=", ind);
    dclistPrint(stderr, mat->S[ind]->next); fprintf(stderr, "\n");
#endif
    // TODO: replace this with a move of pointers...!
#if USE_CONNECT == 0
    mat->A[j] = dclistInsert(mat->S[ind], j);
#else
    dclistConnect(mat->S[ind], mat->A[j]);
#endif
#if DEBUG >= 2
    fprintf(stderr, "S[%d]_a=", ind);
    dclistPrint(stderr, mat->S[ind]->next); fprintf(stderr, "\n");
#endif
    // update R[j] by removing i
    remove_i_from_Rj(mat, i, j);
}

// for all j in row[i], removes j and update data
void
removeRowSWAR(sparse_mat_t *mat, int i)
{
    int k;

#if TRACE_ROW >= 0
    if(i == TRACE_ROW)
	fprintf(stderr, "TRACE_ROW: removeRowSWAR i=%d\n", i);
#endif
    mat->weight -= lengthRow(mat, i);
#if USE_TAB == 0
    for(k = 0; k < lengthRow(mat, i); k++){
#else
    for(k = 1; k <= lengthRow(mat, i); k++){
#endif
#if TRACE_COL >= 0
	if(cell(mat, i, k) == TRACE_COL){
	    fprintf(stderr, "removeRowSWAR removes %d from R_%d\n",
		    TRACE_COL, i);
	}
#endif
	removeCellSWAR(mat, i, cell(mat, i, k));
    }
}

void
addCellSWAR(sparse_mat_t *mat, int i, INT j)
{
    int ind;

#if TRACE_ROW >= 0
    if(i == TRACE_ROW)
	fprintf(stderr, "TRACE_ROW: addCellSWAR i=%d j=%d\n", i, j);
#endif
    // update weight
#if DEBUG >= 1
    fprintf(stderr, "addCellSWAR: moving j=%d from S[%d] to S[%d]\n",
	    j, mat->wt[j], incrS(mat->wt[j]));
#endif
    ind = mat->wt[j] = incrS(mat->wt[j]);
#if USE_MERGE_FAST > 1
    if(ind < 0)
	return;
#endif
    remove_j_from_S(mat, j);
    if(mat->wt[j] > mat->cwmax){
#if DEBUG >= 1
	fprintf(stderr, "WARNING: column %d is too heavy (%d)\n", j,
		mat->wt[j]);
#endif
	ind = mat->cwmax+1; // trick
    }
    // update A[j]
#if USE_CONNECT == 0
    mat->A[j] = dclistInsert(mat->S[ind], j);
#else
    dclistConnect(mat->S[ind], mat->A[j]);
#endif
    // update R[j] by adding i
    add_i_to_Rj(mat, i, j);
}

void
addRowSWAR(sparse_mat_t *mat, int i)
{
    int k;

    mat->weight += lengthRow(mat, i);
#if USE_TAB == 0
    for(k = 0; k < lengthRow(mat, i); k++)
#else
    for(k = 1; k <= lengthRow(mat, i); k++)
#endif
	addCellSWAR(mat, i, cell(mat, i, k));
}

// i1 += i2.
// len could be the real length of row[i1]+row[i2] or -1.
void
addRowsSWAR(sparse_mat_t *mat, int i1, int i2, int len)
{
#if 0 // original version
    int k1, k2, k, len, *tmp, *tmp2;

# if TEX
    fprintf(stderr, "row[%d] += row[%d]\n", i1, i2);
# endif
# if DEBUG >= 1
    fprintf(stderr, "row[%d] =", i1); print_row(mat, i1);
    fprintf(stderr, "\n");
    fprintf(stderr, "row[%d] =", i2); print_row(mat, i2);
    fprintf(stderr, "\n");
# endif
    // merge row[i1] and row[i2] in i1...
    len = mat->data[i1].len + mat->data[i2].len;
    tmp = (INT *)malloc(len * sizeof(INT));
    k = k1 = k2 = 0;

    // i1 is to disappear, replaced by a new one
    removeRowSWAR(mat, i1);
    // loop while everybody is here
    while((k1 < mat->data[i1].len) && (k2 < mat->data[i2].len)){
	if(mat->data[i1].val[k1] < mat->data[i2].val[k2]){
	    tmp[k++] = mat->data[i1].val[k1++];
	    mat->weight++;
	}
	else if(mat->data[i1].val[k1] > mat->data[i2].val[k2]){
            tmp[k++] = mat->data[i2].val[k2++];
	    mat->weight++;
        }
	else{
	    k1++;
	    k2++;
	}
    }
    // finish with k1
    for( ; k1 < mat->data[i1].len; k1++){
	tmp[k++] = mat->data[i1].val[k1];
	mat->weight++;
    }
    // finish with k2
    for( ; k2 < mat->data[i2].len; k2++){
	tmp[k++] = mat->data[i2].val[k2];
	mat->weight++;
    }
    // destroy and copy back
    free(mat->data[i1].val);
    tmp2 = (INT *)malloc(k * sizeof(INT));
    memcpy(tmp2, tmp, k * sizeof(INT));
    mat->data[i1].len = k;
    mat->data[i1].val = tmp2;
    free(tmp);
    addRowSWAR(mat, i1);
# if DEBUG >= 1
    fprintf(stderr, "row[%d]+row[%d] =", i1, i2);
    print_row(mat, i1); fprintf(stderr, "\n");
# endif
#else // cleaner one, that shares addRowsData() to prepare the next move...!
    // i1 is to disappear, replaced by a new one
    removeRowSWAR(mat, i1);
# if USE_TAB == 0
    addRowsData(mat->data, i1, i2);
# else
    // we know the length of row[i1]+row[i2]
    addRows(mat->rows, i1, i2, len);
# endif
    addRowSWAR(mat, i1);
#endif
}

// Remove j from R[i] and crunch R[i].
void
remove_j_from_row(sparse_mat_t *mat, int i, int j)
{
    int k;

#if TRACE_ROW >= 0
    if(i == TRACE_ROW){
	fprintf(stderr, "TRACE_ROW: remove_j_from_row i=%d j=%d\n", i, j);
    }
#endif
#if DEBUG >= 2
    fprintf(stderr, "row[%d]_b=", i);
    print_row(mat, i);
    fprintf(stderr, "\n");
#endif
#if USE_TAB == 0
    for(k = 0; k < lengthRow(mat, i); k++)
	if(cell(mat, i, k) == j)
	    break;
    ASSERT(k < lengthRow(mat, i));
    // crunch
    for(++k; k < lengthRow(mat, i); k++)
	cell(mat, i, k-1) = cell(mat, i, k);
#else
    for(k = 1; k <= lengthRow(mat, i); k++)
	if(cell(mat, i, k) == j)
	    break;
    ASSERT(k <= lengthRow(mat, i));
    // crunch
    for(++k; k <= lengthRow(mat, i); k++)
	cell(mat, i, k-1) = cell(mat, i, k);
#endif
    lengthRow(mat, i) -= 1;
    mat->weight--;
#if DEBUG >= 2
    fprintf(stderr, "row[%d]_a=", i);
    print_row(mat, i);
    fprintf(stderr, "\n");
#endif
}

void
removeRowDefinitely(FILE *outfile, sparse_mat_t *mat, INT i, int doreport)
{
    removeRowSWAR(mat, i);
    destroyRow(mat, i);
    if(doreport)
	report1(outfile, i);
    mat->rem_nrows--;
}

// These columns are simply removed from the current squeleton matrix, but
// not from the small matrix.
int
deleteAllColsFromStack(FILE *outfile, sparse_mat_t *mat, int iS)
{
    dclist dcl;
    INT j;
    int k, njrem = 0;

    while(1){
	dcl = mat->S[iS]->next;
	if(dcl == NULL)
	    break;
	j = dcl->j;
	njrem++;
#if DEBUG >= 1
	fprintf(stderr, "Removing column %d from S[%d]\n", j, iS);
#endif
#if TRACE_COL >= 0
	if(j == TRACE_COL)
	    fprintf(stderr, "Removing column %d from S[%d]\n", j, iS);
#endif
	// we destroy column j
	// TODO: do faster?
#if USE_MERGE_FAST <= 1
	fprintf(stderr, "FIXMEEEEEEEEE.\n");
#else
	if(iS == 0)
	    mat->rem_ncols--;
	if(iS == 1){
	    for(k = 1; k <= mat->R[j][0]; k++)
		if(mat->R[j][k] != -1){
# if TRACE_COL >= 0
		    if(j == TRACE_COL)
			fprintf(stderr, "deleteAllCols: row is %d\n",mat->R[j][k]);
# endif
		    remove_j_from_row(mat, mat->R[j][k], j);
		    removeRowDefinitely(outfile, mat, mat->R[j][k], 1);
#ifdef USE_MPI
		    mpi_send_inactive_cols(mat->R[j][k]);
#endif
		    mat->rem_ncols--;
		}
	    mat->wt[j] = 0;
	}
	k = mat->wt[j]; // make a copy of the weight
	remove_j_from_SWAR(mat, j);
	// mat->wt[j] is put to 0...
	mat->wt[j] = -k; // restore and update
#endif
    }
#if DEBUG >= 1
    if(njrem > 0)
	fprintf(stderr, "deleteAllColsFromStack[%d]: %d\n", iS, njrem);
#endif
    return njrem;
}

int
deleteHeavyColumns(FILE *outfile, sparse_mat_t *mat)
{
    return deleteAllColsFromStack(outfile, mat, mat->cwmax+1);
}

int
deleteEmptyColumns(sparse_mat_t *mat)
{
#if 0
    return deleteAllColsFromStack(mat, 0);
#else
    dclist dcl = mat->S[0], foo;
    int njrem = 0;
    INT j;

    while(dcl->next != NULL){
	foo = dcl->next;
	j = foo->j;
	dcl->next = foo->next;
	free(foo);
	njrem++;
	mat->A[j] = NULL;
	mat->wt[j] = 0;
	destroyRj(mat, j);
    }
    mat->rem_ncols -= njrem;
    return njrem;
#endif
}

int
removeSingletons(FILE *outfile, sparse_mat_t *mat)
{
#if USE_MERGE_FAST == 0
    INT i, j;

    for(j = 0; j < mat->ncols; j++)
	if(mat->wt[j] == 1){
	    // find row...
	    for(i = 0; i < mat->nrows; i++)
		if(hasCol(mat, i, j))
		    break;
	    ASSERT(i < mat->nrows);
	    removeWeightFromRow(mat, i);
	    destroyRow(mat, i);
	    report1(outfile, i); // signal to replay...!
	}
#else
    return deleteAllColsFromStack(outfile, mat, 1);
#endif
}

// try all combinations to find the smaller one
void
tryAllCombinations(FILE *outfile, sparse_mat_t *mat, int m, INT *ind)
{
    int i, k;

#ifdef USE_MPI
    // naive strategy, where the pivot is the row of smallest weight
    i = 0;
    for(k = 1; k < m; k++)
	if(lengthRow(mat, ind[k]) < lengthRow(mat, ind[i]))
	    i = k;
#else
    i = findBestIndex(mat, m, ind);
#endif
#if DEBUG >= 1
    fprintf(stderr, "Minimal is i=%d (%d %d)\n", i, ind[0], ind[1]);
#endif
    for(k = 0; k < m; k++)
	if(k != i){
	    addRowsSWAR(mat, ind[k], ind[i], -1);
#if DEBUG >= 1
	    fprintf(stderr, "new row[%d]=", ind[k]);
	    print_row(mat, ind[k]);
	    fprintf(stderr, "\n");
#endif
	}
    if(i > 0){
	// put ind[i] in front
	int itmp = ind[i];
	for(k = i; k > 0; k--)
	    ind[k] = ind[k-1];
	ind[0] = itmp;
    }
#if DEBUG >= 1
    fprintf(stderr, "=> new_ind: %d %d\n", ind[0], ind[1]);
#endif
    reportn(outfile, ind, m);
    removeRowSWAR(mat, ind[0]);
    destroyRow(mat, ind[0]);
}

void
addFatherToSonsHeight(FILE *outfile, sparse_mat_t *mat, int m, int *ind,
		      int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX],
		      int *father, int *height, int hmax)
{
    int adds[MERGE_LEVEL_MAX << 1], nadds, itab;
    INT tab[MERGE_LEVEL_MAX << 1];
    int i, h, nV = m;

    // add each father to each of his sons, but from bottom to top...!
    for(h = hmax; h > 0; h--){
	// we are currently at height h
	if(nV == 1)
	    // only the root remains!
	    break;
	nadds = 0;
	for(i = 1; i < m; i++)
	    if(height[i] == h){
		int i1 = ind[i], i2 = ind[father[i]];
#if DEBUG >= 1
		fprintf(stderr, "h = %d, i = %d: ", h, i);
		fprintf(stderr, "row[%d] += row[%d]\n", i1, i2);
#endif
		nV--;

		// we know the length of row[i1]+row[i2]...!
		addRowsSWAR(mat, i1, i2, A[i][father[i]]);

		// since we want to add without destroying, we report
		// -(i2+1) i1 // hack!!!!
#if 0 // for non compact reporting + comment reportn for adds[]
		report2(-(i2+1), i1);
#endif
		adds[nadds++] = -(i2+1);
		adds[nadds++] = i1;
	    }
	qsort(adds, nadds>>1, 2 * sizeof(int), cmp);
#if DEBUG >= 1
	fprintf(stderr, "all adds[%d]:", h);
	for(i = 0; i < nadds; i++)
	    fprintf(stderr, " %d", adds[i]);
	fprintf(stderr, "\n");
#endif
	// compress history
	tab[0] = adds[0];
	tab[1] = adds[1];
	itab = 2;
	for(i = 2; i < nadds; i += 2){
	    if(adds[i] == tab[0])
		// same father
		tab[itab++] = adds[i+1];
	    else{
		// new father
		reportn(outfile, tab, itab);
		tab[0] = adds[i];
		tab[1] = adds[i+1];
		itab = 2;
	    }
	}
	if(h == 1)
	    // we are at the root, actually
	    tab[0] = -tab[0]-1; // what a trick!
	reportn(outfile, tab, itab);
    }
}

// add u to its sons
void
addFatherToSonsRec(FILE *outfile, sparse_mat_t *mat, int m, int *ind,
		   int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX],
		   int *father,
		   int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		   int u)
{
    int k, tab[MERGE_LEVEL_MAX], itab = 0, i1, i2;

    if(sons[u][0] == 0)
	// nothing to do for a leaf...!
	return;
    i2 = ind[u];
    if(u == father[u])
	tab[itab++] = i2; // trick for the root!!!
    else
	tab[itab++] = -(i2+1); // the usual trick not to destroy row
    for(k = 1; k <= sons[u][0]; k++){
	addFatherToSonsRec(outfile, mat, m, ind, A, father, sons, sons[u][k]);
	i1 = ind[sons[u][k]];
	// add u to its son
	addRowsSWAR(mat, i1, i2, A[sons[u][k]][u]);
	tab[itab++] = i1;
    }
    reportn(outfile, tab, itab);
}

/* FIXME: Maybe height and hmax should go ? */
void
addFatherToSons(FILE *outfile, sparse_mat_t *mat, int m, int *ind,
		int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX],
		int *father,
                int *height MAYBE_UNUSED, int hmax MAYBE_UNUSED,
		int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1])
{
#if 0
    addFatherToSonsHeight(mat, m, ind, A, father, height, hmax);
#else
    addFatherToSonsRec(outfile, mat, m, ind, A, father, sons, 0);
#endif
}

void
useMinimalSpanningTree(FILE *outfile, sparse_mat_t *mat, int m,
		       INT *ind, double *tfill, double *tMST)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX];
    int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];
    int father[MERGE_LEVEL_MAX], height[MERGE_LEVEL_MAX], hmax;

    *tfill = seconds();
    fillRowAddMatrix(A, mat, m, ind);
    *tfill = seconds()-*tfill;
    *tMST = seconds();
    hmax = minimalSpanningTree(father, height, sons, m, A);
    *tMST = seconds()-*tMST;
#if DEBUG >= 1
    {
      int i;
      for (i = 0; i < m; i++)
	fprintf (stderr, "father[%d] = %d\n", i, father[i]);
      for (i = 0; i < m; i++)
        {
          int k;
          fprintf (stderr, "Sons of %d:", i);
          for (k = 1; k <= sons[i][0]; k++)
	    fprintf (stderr, " %d", sons[i][k]);
          fprintf (stderr, "\n");
        }
    }
#endif
    addFatherToSons(outfile, mat, m, ind, A, father, height, hmax, sons);
#if 0 // for non compact reporting
    // delete root!
    report1(ind[0]);
#endif
    removeRowSWAR(mat, ind[0]);
    destroyRow(mat, ind[0]);
}

void
findOptimalCombination(FILE *outfile, sparse_mat_t *mat, int m, INT j,
		       INT *ind, double *tfill, double *tMST)
{
#ifdef USE_MPI
    // tmp??!?!?!?
    *tfill = 0.0;
    *tMST = 0.0;
    tryAllCombinations(outfile, mat, m, ind);
    mpi_add_rows(mat, m, j, ind);
#else
    if(m <= 2){
	*tfill = *tMST = 0;
	tryAllCombinations(outfile, mat, m, ind);
    }
    else{
# if 0
	tryAllCombinations(outfile, mat, m, ind);
# else
	useMinimalSpanningTree(outfile, mat, m, ind, tfill, tMST);
# endif
    }
#endif
}

void
checkWeight(sparse_mat_t *mat, INT j)
{
    int i, w = 0;

    fprintf(stderr, "Rows containing %ld:", (long int) j);
    for(i = 0; i < mat->nrows; i++)
	if(!isRowNull(mat, i))
	    if(hasCol(mat, i, j)){
		fprintf(stderr, " %d", i);
		w++;
	    }
    fprintf(stderr, "\n");
    ASSERT(w == (mat->wt[j] >= 0 ? mat->wt[j] : -mat->wt[j]));
}

void
checkMatrixWeight(sparse_mat_t *mat)
{
    unsigned long w = 0;
    int i;

    for(i = 0; i < mat->nrows; i++)
	if(!isRowNull(mat, i))
	    w += lengthRow(mat, i);
    ASSERT(w == mat->weight);
}

// j has weight m.
void
mergeForColumn(FILE *outfile, double *tt, double *tfill, double *tMST,
	       sparse_mat_t *mat, int m, INT j)
{
    INT ind[MERGE_LEVEL_MAX];
    int ni, k;

#if DEBUG >= 1
    fprintf(stderr, "Treating column %d of weight %d\n", j, mat->wt[j]);
#endif
#if (DEBUG >= 1) || TEX
    fprintf(stderr, "Status before next j=%d to start\n", j);
# if TEX
    texSWAR(mat);
# else
    printSWAR(mat, mat->ncols);
# endif
#endif
    // the corresponding rows are in R[j], skipping 1st cell and -1's
#if DEBUG >= 2
    checkCoherence(mat, m, j);
#endif
    ni = 0;
    for(k = 1; k <= mat->R[j][0]; k++){
	if(mat->R[j][k] != -1){
	    ind[ni++] = mat->R[j][k];
	    if(ni == m)
		break;
	}
    }
#if DEBUG >= 1
    fprintf(stderr, " %d", j);
    fprintf(stderr, "=> the %d rows are:\n", m);
    for(k = 0; k < m; k++){
	fprintf(stderr, "row[%d]=", ind[k]);
	print_row(mat, ind[k]);
	fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
#endif
    *tt = seconds();
    findOptimalCombination(outfile, mat, m, j, ind, tfill, tMST);
    *tt = seconds()-(*tt);
    mat->rem_nrows--;
    mat->rem_ncols--;
    remove_j_from_SWAR(mat, j);
}

int
merge_m_fast(FILE *outfile, sparse_mat_t *mat, int m, int verbose)
{
    double totopt=0, tot=seconds(), tt, totfill=0, totMST=0, tfill, tMST;
    double totdel = 0;
    int njproc = 0, njrem = 0;
    INT j;
    dclist dcl = mat->S[m];
    int report = (verbose == 0) ? 10000 : 1000;

#if DEBUG >= 1
    fprintf(stderr, "Weight %d:", m);
#endif

    //    checkWeight(mat, 132461);

    while(1){
#if DEBUG >= 1
	checkMatrixWeight(mat);
	checkData(mat);
#endif
	// get next j
	dcl = mat->S[m]->next;
	if(dcl == NULL)
	    break;
	j = dcl->j;
	njproc++;
	if(!(njproc % report))
	    fprintf(stderr, "# %d columns of weight %d processed at %2.2lf\n",
		    njproc, m, seconds());
	mergeForColumn(outfile, &tt, &tfill, &tMST, mat, m, j);
#if TRACE_COL >= 0
	fprintf(stderr, "wt[%d]=%d after findOptimalCombination\n",
		TRACE_COL, mat->wt[TRACE_COL]);
#endif
	totopt += tt;
	totfill += tfill;
	totMST += tMST;
	tt = seconds();
	njrem += deleteHeavyColumns(outfile, mat);
	totdel += (seconds()-tt);
#if TRACE_COL >= 0
        fprintf(stderr, "wt[%d]=%d after deleteHeavyColumns\n",
                TRACE_COL, mat->wt[TRACE_COL]);
#endif
#ifdef USE_MPI
	// better, otherwise an incoherent state may appear between
	// master, sub-master and slaves, like the current proc continuing
	// eliminating without info from the master
	return njrem;
#endif	
    }
    tot = seconds()-tot;
    fprintf(stderr, "T=%d m=%d nj=%d", (int)seconds(), m, njproc);
    fprintf(stderr, " findopt=%2.2lf (fill=%2.2lf mst=%2.2lf) tot=%2.2lf", 
	    totopt, totfill, totMST, tot);
    fprintf(stderr, " del=%2.2lf\n", totdel);
	    
#if DEBUG >= 1
    fprintf(stderr, "Status at the end of the process\n");
    texSWAR(mat);
#endif
    return njrem;
}

int
merge_m(FILE *outfile, sparse_mat_t *mat, int m, int verbose)
{
#if USE_MERGE_FAST >= 1
  return merge_m_fast(outfile, mat, m, verbose);
#else
  return merge_m_slow(outfile, mat, m, verbose);
#endif
}

// TODO: use mergemax?
int
mergeGe2(FILE *outfile, sparse_mat_t *mat, int m, int verbose)
{
#if USE_MERGE_FAST == 0
    int j, nbm;

    // TODO: remove this and start directly merge_m?
    for(j = 0, nbm = 0; j < mat->ncols; j++)
	if(mat->wt[j] == m){
# if DEBUG >= 1
	    fprintf(stderr, "# wt[%d] = %d\n", j, m);
# endif
	    nbm++;
	}
    fprintf(stderr, "There are %d column(s) of weight %d\n", nbm, m);
    if(nbm)
#endif
	return merge_m(outfile, mat, m, verbose);
#if DEBUG >= 1
    matrix2tex(mat);
#endif
}

int
minColWeight(sparse_mat_t *mat)
{
#if USE_MERGE_FAST == 0
    int j, minw = mat->nrows;

    for(j = 0; j < mat->ncols; j++)
	if((mat->wt[j] > 0) && (mat->wt[j] < minw))
	    minw = mat->wt[j];
    return minw;
#else
    // in this case, locating minw is easy: inspect the stacks!
    int m;

    for(m = 1; m <= mat->cwmax; m++)
	if(mat->S[m]->next != NULL)
	    return m;
    return -1;
#endif
}

int
inspectRowWeight(FILE *outfile, sparse_mat_t *mat)
{
    //    double tt = seconds();
    int i, maxw = 0, nirem = 0, niremmax = 128;

    for(i = 0; i < mat->nrows; i++){
	if(!isRowNull(mat, i)){
	    if(lengthRow(mat, i) > mat->rwmax){
		if((mat->rem_nrows - mat->rem_ncols) > mat->delta){
#if DEBUG >= 1
		    fprintf(stderr, "Removing too heavy row[%d]: %d\n", 
			    i, lengthRow(mat, i));
#endif
#if TRACE_ROW >= 0
		    if(i == TRACE_ROW)
			fprintf(stderr, 
				"TRACE_ROW: removing too heavy row[%d]: %d\n", 
				i, lengthRow(mat, i));
#endif
		    removeRowDefinitely(outfile, mat, i, 1);
		    nirem++;
		    if(!(nirem % 10000))
			fprintf(stderr, "#removed_rows=%d at %2.2lf\n",
				nirem, seconds());
		    if(nirem > niremmax)
			break;
		}
	    }
	    if(!isRowNull(mat, i) && (lengthRow(mat, i) > maxw))
		maxw = lengthRow(mat, i);
	}
    }
    //    fprintf(stderr, "nirem=%d; nrows=%d; max row weight is %d (%2.2lf)\n", 
    //	    nirem, mat->rem_nrows, maxw, seconds()-tt);
    return nirem;
}

// Delete superfluous rows s.t. nrows-ncols >= keep.
// It could be fun to find rows s.t. at least one j is of small weight.
// Or components of rows with some sharing?
// We could define sf(row) = sum_{j in row} w(j) and remove the
// rows of smallest value of sf?
// For the time being, remove heaviest rows.
int
deleteSuperfluousRows(FILE *outfile, sparse_mat_t *mat, int keep, int niremmax)
{
    int nirem = 0, *tmp, ntmp, itmp, i;

    if((mat->rem_nrows - mat->rem_ncols) <= keep)
	return 0;
    ntmp = mat->rem_nrows << 1;
    tmp = (int *)malloc(ntmp * sizeof(int));
    for(i = 0, itmp = 0; i < mat->nrows; i++)
        if(!isRowNull(mat, i)){
	    tmp[itmp++] = lengthRow(mat, i);
	    tmp[itmp++] = i;
	}
    qsort(tmp, ntmp>>1, 2 * sizeof(int), cmp);
    for(i = ntmp-1; i >= 0; i -= 2){
	if((nirem >= niremmax) || (mat->rem_nrows - mat->rem_ncols) <= keep)
	    break;
	removeRowDefinitely(outfile, mat, tmp[i], 1);
	nirem++;
    }
    free(tmp);
    return nirem;
}

void
merge(FILE *outfile, sparse_mat_t *mat, int maxlevel, int verbose, int forbw)
{
    unsigned long bwcostmin = 0, oldcost = 0, cost;
    int old_nrows, old_ncols, m, mm, njrem = 0, ncost = 0, ncostmax;

    ncostmax = 20; // was 5
    m = 2;
    while(1){
	cost = ((unsigned long)mat->rem_ncols) * ((unsigned long)mat->weight);
	fprintf(stderr, "w(M)=%lu, w(M)*ncols=%lu", mat->weight, cost);
	fprintf(stderr, " w(M)/ncols=%2.2lf\n", 
		((double)mat->weight)/((double)mat->rem_ncols));
	if((bwcostmin == 0) || (cost < bwcostmin))
	    bwcostmin = cost;
	if(forbw)
	    // what a trick!!!!
	    fprintf(outfile, "BWCOST: %lu\n", cost);
	if(forbw && (oldcost != 0) && (cost > oldcost)){
	    ncost++;
	    fprintf(stderr, "WARNING: New cost > old cost (%2.6e) [%d/%d]\n",
		    ((double)cost)-((double)oldcost), ncost, ncostmax);
	    if(ncost >= ncostmax){
		int nirem;

		fprintf(stderr, "WARNING: New cost > old cost %d times",
                        ncost);
		fprintf(stderr, " in a row:");
		nirem = deleteSuperfluousRows(outfile, mat, mat->delta, 128);
		if(nirem == 0){
		    fprintf(stderr, " stopping\n");
		    break;
		}
		else{
		    fprintf(stderr, " try again after removing %d rows!\n",
			    nirem);
		    continue;
		}
	    }
	}
	else
	    ncost = 0;
	oldcost = cost;
	fprintf(stderr, "Performing merge %d at %2.2lf\n", m, seconds());
	old_nrows = mat->rem_nrows;
	old_ncols = mat->rem_ncols;
	if(m == 1)
	    njrem += removeSingletons(outfile, mat);
	else
            njrem += mergeGe2(outfile, mat, m, verbose);
#if DEBUG >= 1
	checkData(mat);
#endif
#if 0
    {
	int i, ni = 0, j0 = 191342;

	for(i = 0; i < mat->nrows; i++)
	    if((mat->data[i].val != NULL) && hasCol(mat, i, j0))
		ni++;
	printf("CHECK: %d %d\n", ni, mat->wt[j0]);
    }
#endif
	fprintf(stderr, "=> nrows=%d ncols=%d (%d) njrem=%d\n",
		mat->rem_nrows, mat->rem_ncols, 
		mat->rem_nrows - mat->rem_ncols, njrem);
	inspectRowWeight(outfile, mat);
	deleteEmptyColumns(mat);
	mm = minColWeight(mat);
	fprintf(stderr, "Min col weight = %d\n", mm);
#if M_STRATEGY == 2
	// jump to the next minimal merge level immediately
	m = mm;
	if((m > maxlevel) || (m <= 0))
	    break;
#else
#  if M_STRATEGY == 1
	if(mm < m)
	    // something new happened, anyway
	    m = mm;
	else
#  endif
	    if((old_nrows == mat->rem_nrows) && (old_ncols == mat->rem_ncols)){
		// nothing happened this time and mm > m
		m = mm;
		if((m > maxlevel) || (m <= 0))
		    break;
	    }
#endif
    }
    if(forbw){
	fprintf(outfile, "BWCOSTMIN: %lu\n", bwcostmin);
	fprintf(stderr, "Minimal bwcost found: %lu\n", bwcostmin);
    }
}

static unsigned long
my_cost(unsigned long N, unsigned long c, int forbw)
{
    if(forbw == 2){
	double K1 = .19e-9, K2 = 3.4e-05, K3 = 1.4e-10; // kinda average
	double dN = (double)N, dc = (double)c;

	return (unsigned long)((K1+K3)*dN*dc+K2*dN*log(dN)*log(dN));
    }
    else if(forbw == 3)
	return c/N;
    else if(forbw <= 1)
	return (N*c);
    return 0;
}

void
doOneMerge(FILE *outfile, sparse_mat_t *mat, int *njrem, double *totopt, double *totfill, double *totMST, double *totdel, int m, int verbose)
{
    dclist dcl;
    double tt, tfill, tMST;
    int j;

    if(m == 1){
	// do all of these
	*njrem += removeSingletons(outfile, mat);
	fprintf(stderr, "T=%d after m=1\n", (int)(seconds()));
    }
    else if(m == 2){
#if 0
	fprintf(stderr, "Performing all merges for m=%d at %2.2lf\n",
		m, seconds());
#endif
	*njrem += mergeGe2(outfile, mat, m, verbose);
    }
    else{
	//	    fprintf(stderr, "Performing one merge for m=%d\n", m);
	dcl = mat->S[m]->next;
	j = dcl->j;
	mergeForColumn(outfile, &tt, &tfill, &tMST, mat, m, j);
	*totopt += tt;
	*totfill += tfill;
	*totMST += tMST;
	tt = seconds();
	*njrem += deleteHeavyColumns(outfile, mat);
	*totdel += (seconds()-tt);
    }
}

void
mergeOneByOne(FILE *outfile, sparse_mat_t *mat, int maxlevel, int verbose, int forbw, double ratio, int coverNmax)
{
    double totopt = 0.0, totfill = 0.0, totMST = 0.0, totdel = 0.0;
    unsigned long bwcostmin = 0, oldbwcost = 0, bwcost = 0;
    int old_nrows, old_ncols, m = 2, njrem = 0, ncost = 0, ncostmax, njproc;
    int mmax = 0, target = 10000;

    fprintf(stderr, "Using mergeOneByOne\n");
    ncostmax = 20; // was 5
    njproc = 0;
    while(1){
	oldbwcost = bwcost;
	old_nrows = mat->rem_nrows;
	old_ncols = mat->rem_ncols;
	m = minColWeight(mat);
	if(m > mmax)
	    mmax = m;
	if(mmax > maxlevel){
	    fprintf(stderr, "maxlevel reached, stopping!\n");
	    break;
	}
	if(m == -1){
	    fprintf(stderr, "All stacks empty, stopping!\n");
	    break;
	}
	doOneMerge(outfile, mat, &njrem, &totopt, &totfill, &totMST, &totdel,
		   m, verbose);
	// number of columns removed
	njproc += old_ncols - mat->rem_ncols;
	deleteEmptyColumns(mat);
	if((old_nrows == mat->rem_nrows) && (old_ncols == mat->rem_ncols)){
	    // is this supposed to be activated at some point? Really?
	    if((m > maxlevel) || (m <= 0))
		break;
	}
	bwcost = my_cost((unsigned long)mat->rem_nrows,
			 (unsigned long)mat->weight, forbw);
	if(njproc >= target){ // somewhat arbitrary...!
	    int kappa = (mat->rem_nrows-mat->rem_ncols) / mat->delta, ni2rem;

	    if(kappa <= (1<<4))
		ni2rem = mat->delta/2;
	    else if(kappa <= (1<<5))
		ni2rem = mat->delta * (kappa/4);
	    else if(kappa <= (1<<10))
		ni2rem = mat->delta * (kappa/8);
	    else if(kappa <= (1<<15))
		ni2rem = mat->delta * (kappa/16);
	    else if(kappa <= (1<<20))
		ni2rem = mat->delta * (kappa/32);
	    else
		ni2rem = mat->delta * (kappa/64);
	    deleteSuperfluousRows(outfile, mat, mat->delta, ni2rem);
	    inspectRowWeight(outfile, mat);
	    fprintf(stderr, "T=%d mmax=%d nr=%d N=%d (%d)",
		    (int)seconds(),
		    mmax,
		    mat->rem_nrows, mat->rem_ncols, 
		    mat->rem_nrows - mat->rem_ncols);
	    if(forbw == 2)
		fprintf(stderr, " bw=%lu", bwcost);
	    else if(forbw == 3)
		fprintf(stderr, " cN=%lu", 
			((unsigned long)mat->rem_nrows)
			*((unsigned long)mat->weight));
	    else if(forbw <= 1)
		fprintf(stderr, " cN=%lu", bwcost);
	    fprintf(stderr, " c/N=%2.2lf\n", 
		    ((double)mat->weight)/((double)mat->rem_ncols));
	    // njrem=%d at %2.2lf\n",
	    if((forbw != 0) && (forbw != 3))
		// what a trick!!!!
		fprintf(outfile, "BWCOST: %lu\n", bwcost);
	    target = njproc + 10000;
	}
	if((bwcostmin == 0) || (bwcost < bwcostmin)){
	    bwcostmin = bwcost;
	    if((forbw != 0) && (forbw != 3))
		// what a trick!!!!
		fprintf(outfile, "BWCOST: %lu\n", bwcost);
	}
	// to be cleaned one day...
	if((forbw == 0) || (forbw == 2)){
	    double r = ((double)bwcost)/((double)bwcostmin);
	    if((r > ratio) && (mat->rem_nrows-mat->rem_ncols <= mat->delta)){
		if(forbw == 0)
		    fprintf(stderr, "cN too high, stopping [%2.2lf]\n", r);
		else
		    fprintf(stderr, "bw too high, stopping [%2.2lf]\n", r);
		break;
	    }
	}
	else if(forbw == 3){
	    if(bwcost > (unsigned long)coverNmax){
		fprintf(stderr, "c/N too high, stopping [%lu]\n", bwcost);
		break;
	    }
	}
	if((forbw != 0) && (oldbwcost != 0) && (bwcost > oldbwcost)){
	    ncost++;
#if 0
	    fprintf(stderr, "New cost > old cost (%2.6e) [%d/%d]\n",
		    ((double)bwcost)-((double)oldbwcost), ncost, ncostmax);
#endif
	    if(ncost >= ncostmax){
		int nirem;

		fprintf(stderr, "New cost > old cost %d times", ncost);
		fprintf(stderr, " in a row:");
		nirem = deleteSuperfluousRows(outfile, mat, mat->delta, 128);
		if(nirem == 0){
		    fprintf(stderr, " stopping\n");
		    break;
		}
		else{
		    fprintf(stderr, " try again after removing %d rows!\n",
			    nirem);
		    njproc += nirem; // humf: odd name for njproc...!
		    continue;
		}
	    }
	}
	else
	    ncost = 0;
    }
    deleteSuperfluousRows(outfile, mat, mat->delta, 
			  mat->rem_nrows-mat->rem_ncols+mat->delta);
    if((forbw != 0) && (forbw != 3)){
	fprintf(outfile, "BWCOSTMIN: %lu\n", bwcostmin);
	fprintf(stderr, "Minimal bwcost found: %lu\n", bwcostmin);
    }
}

void
checkmat(sparse_mat_t *mat)
{
    int i;

    for(i = 0; i < mat->nrows; i++){
	fprintf(stderr, "Row %d =", i);
	print_row(mat, i);
	fprintf(stderr, "\n");
    }
    exit(-1);
}

void
dumpSparse(FILE *ofile, sparse_mat_t *mat)
{
    long w = 0;
    int i, j, k, buf[1000], ibuf, new_nrows = 0;

    fprintf(ofile, "%d %d\n", mat->rem_nrows, mat->rem_ncols);
    for(i = 0; i < mat->nrows; i++){
	if(isRowNull(mat, i))
	    continue;
	w += lengthRow(mat, i);
	new_nrows++;
	ibuf = 0;
#if USE_TAB == 0
	for(k = 0; k < lengthRow(mat, i); k++)
#else
	for(k = 1; k <= lengthRow(mat, i); k++)
#endif
	    buf[ibuf++] = cell(mat, i, k);
	fprintf(ofile, "%d", ibuf);
	for(k = 0; k < ibuf; k++)
	    fprintf(ofile, " %d", buf[k]);
	fprintf(ofile, "\n");
    }
    fprintf(stderr, "Weight = %ld\n", w);
    ASSERT(new_nrows == mat->rem_nrows);
    // check j used
    for(j = 0; j < mat->ncols; j++)
	if(mat->wt[j] != 0){
	    fprintf(ofile, "J %d %d\n", j, abs(mat->wt[j]));
	    if(abs(mat->wt[j]) <= mat->mergelevelmax)
		fprintf(stderr, "Gasp: %d %d\n", j, mat->wt[j]);
	}
}

#ifdef USE_MPI
/*
  MPI section
 */

#include "mpi.h"

#define MPI_ERROR_TAG         -1
#define MPI_DIE_TAG            1
#define MPI_J_TAG              2
#define MPI_ADD_TAG            3
#define MPI_SEND_W_TAG         4
#define MPI_SEND_ROW_TAG       5
#define MPI_SEND_ROW_BACK_TAG  6
#define MPI_MIN_M              7
#define MPI_DO_M               8
#define MPI_M_DONE             9
#define MPI_INACTIVATE_ROWS       10

#define MPI_BUF_SIZE 1000

int *mpi_tab_j;

void
mpi_err(char *str)
{
    FILE *out = stderr;
    int mpi_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    fprintf(out, "MPI#%d# %s", mpi_rank, str);
    fflush(out);
}

void
mpi_err1(char *format, int i)
{
    FILE *out = stderr;
    int mpi_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    fprintf(out, "MPI#%d# ", mpi_rank);
    fprintf(out, format, i);
    fflush(out);
}

void
mpi_err2(char *format, int i1, int i2)
{
    FILE *out = stderr;
    int mpi_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    fprintf(out, "MPI#%d# ", mpi_rank);
    fprintf(out, format, i1, i2);
    fflush(out);
}

void
mpi_err3(char *format, int i1, int i2, int i3)
{
    FILE *out = stderr;
    int mpi_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    fprintf(out, "MPI#%d# ", mpi_rank);
    fprintf(out, format, i1, i2, i3);
    fflush(out);
}

int
mpi_get_proc_for_j(int j)
{
    int i;

    for(i = 0; j >= mpi_tab_j[i]; i++);
    // j < mpi_tab_j[i], we hope
    return i;
}

void
mpi_kill_slaves()
{
    unsigned int imsg[1];
    int nprocs, rank;
    
    printf("MPI# Stopping everybody\n");
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    imsg[0] = 0;
    for(rank = 1; rank < nprocs; ++rank)
	MPI_Send(imsg, 1, MPI_UNSIGNED, rank, MPI_DIE_TAG, MPI_COMM_WORLD);
    printf("MPI# I sent everybody the die signal\n");
}

void
mpi_check_rows(sparse_mat_t *mat, INT i, INT *tab, int ntab)
{
    int k1, k2;
#if DEBUG >= 1
    int k;
    fprintf(stderr, "Original:         ");
    for(k = 1; k <= lengthRow(mat, i); k++)
	fprintf(stderr, " %d", cell(mat, i, k));
    fprintf(stderr, "\n");
    fprintf(stderr, "Reconstructed[%d]:", ntab);
    for(k = 0; k < ntab; k++)
	fprintf(stderr, " %d", tab[k]);
    fprintf(stderr, "\n");
    for(k = 0; k < ntab; k++){
	// check that tab[k] is indeed a position in mat->rows[i]
	int ok = 0, r;
	for(r = 1; r <= lengthRow(mat, i); r++)
	    if(cell(mat, i, r) == (int)tab[k]){
		ok = 1;
		break;
	    }
	if(!ok){
	    fprintf(stderr, "PB? for i=%d jj=%u\n", i, tab[k]);
	    exit(0);
	}
    }
#endif
    k1 = 0;
    k2 = 1;
    while(k1 < ntab){
	if(cell(mat, i, k2) == tab[k1]){
	    k2++; k1++;
	}
	else if(cell(mat, i, k2) == -1)
	    k2++;
	else{
	    fprintf(stderr, "PB? for i=%d/k1=%d/k2=%d\n", i, k1, k2);
	    exit(0);
	}
    }
    for( ; k2 <= lengthRow(mat, i); k2++)
	if(cell(mat, i, k2) != -1){
	    fprintf(stderr, "Non -1 cell!!!!\n");
	    exit(0);
	}
}

#define MAX_ROW_LENGTH 1000 // who cares, really?

void
mpi_load_rows_for_j(sparse_mat_t *mat, int m, INT j)
{
    MPI_Status status;
    unsigned int buf[MPI_BUF_SIZE];
    int mpi_size, ibuf, k, nrecv, cnt, ind;
    INT i, **newrows;
    char *done;

    fprintf(stderr, "Loading m=%d rows for j=%d\n", m, j);
    buf[0] = (unsigned)m;
    buf[1] = (unsigned)j;
    ibuf = 2;
    for(k = 1; k <= mat->R[j][0]; k++)
        if(mat->R[j][k] != -1)
	    buf[ibuf++] = (unsigned)mat->R[j][k];
    // at this point, we should have ibuf-2 == m...!
    if((ibuf-2) != m){
	fprintf(stderr, "#!# ibuf-2 != m in mpi_load_rows_for_j\n");
	exit(0);
    }
#if DEBUG >= 1
    fprintf(stderr, "Ready to send R[%d]=%d // #i=%d\n", j, 
	    mat->R[j][0], ibuf-2);
    fflush(stderr);
#endif
    // asking for the rows referenced by mat->R[j]
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    done = (char *)malloc(mpi_size * sizeof(char));
    for(k = 1; k < mpi_size; k++){
	MPI_Send(buf, ibuf, MPI_UNSIGNED, k, MPI_SEND_ROW_TAG, MPI_COMM_WORLD);
	done[k] = 0;
    }
    // feed newrows
    newrows = (INT **)malloc(m * sizeof(INT *));
    for(k = 0; k < m; k++){
	newrows[k] = (INT *)malloc(MAX_ROW_LENGTH * sizeof(INT));
	newrows[k][0] = (INT)buf[k+2];
	newrows[k][1] = 1; // last index fed
    }
    nrecv = 0;
    // waiting for the rows to come back from the slaves
    while(1){
# if DEBUG >= 1
	fprintf(stderr, "[m_l_r] Master waiting... nrecv=%d\n", nrecv);
	fflush(stderr);
# endif
	MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
# if DEBUG >= 1
	fprintf(stderr, "Something received: %d\n", status.MPI_TAG);
	fflush(stderr);
# endif
	if(status.MPI_TAG == MPI_SEND_ROW_BACK_TAG){
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
# if DEBUG >= 1
	    fprintf(stderr, "[m_l_r] source=%u:", status.MPI_SOURCE);
	    for(k = 0; k < cnt; k++)
		fprintf(stderr, " %u", buf[k]);
	    fprintf(stderr, "\n");
# endif
	    // buf = [m, j, i, j_1, ..., j_r] where the indices j_s
	    // are positions in mat->rows[i]
	    i = (INT)buf[2];
	    if(isRowNull(mat, i)){
		fprintf(stderr, "Can it really happen??? i=%d j=%u cnt=%d\n", 
			i, buf[1], cnt);
		fflush(stderr);
	    }
	    // feed the corresponding row
	    for(ind = 0; ind < m; ind++)
		if(newrows[ind][0] == i)
		    break;
	    for(k = 3; k < cnt; k++){
		newrows[ind][1]++;
		if(newrows[ind][1] >= MAX_ROW_LENGTH){
		    fprintf(stderr, "#!# newrows.length exceeded, sorry!\n");
		    exit(0);
		}
		newrows[ind][newrows[ind][1]] = (INT)buf[k];
	    }
	    done[status.MPI_SOURCE] += 1;
	    if(done[status.MPI_SOURCE] == m){
		nrecv++;
		if(nrecv == (mpi_size-1))
		    break;
	    }
	}
    }
    for(k = 0; k < m; k++){
	qsort(newrows[k]+2, newrows[k][1]-1, sizeof(INT), cmp);
	mpi_check_rows(mat, newrows[k][0], newrows[k]+2, newrows[k][1]-1);
    }
    free(done);
    for(k = 0; k < m; k++)
	free(newrows[k]);
    free(newrows);
}

void
mpi_send_inactive_cols(int i)
{
    unsigned int buf[MPI_BUF_SIZE];
    int mpi_size, mpi_rank, k;

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    buf[0] = mpi_index;
    buf[1] = i;
    // TODO_MPI: include the master?
    for(k = 1; k < mpi_size; k++)
	if(k != mpi_rank)
	    MPI_Send(buf, 2, MPI_UNSIGNED, k, 
		     MPI_INACTIVATE_ROWS, MPI_COMM_WORLD);
}

void
mpi_inactivate_rows(FILE *outfile, sparse_mat_t *mat, unsigned int *tab, int ntab)
{
    int k, i;

    for(k = 1; k < ntab; k++){
	i = (int)tab[k];
	mpi_err2("index=%d Row[%d] has to be inactivated\n", tab[0], i);
	if(isRowNull(mat, i))
	    mpi_err1("Row[%d] already nulled!!!\n", i);
	else
	    removeRowDefinitely(outfile, mat, i, 0);
    }
}

int
mpi_get_number_of_active_colums(sparse_mat_t *mat, INT jmin, INT jmax)
{
    INT j;
    int nb = 0;

    for(j = jmin; j < jmax; j++)
	nb += (mat->wt[j] != 0);
    return nb;
}

// Sending rows to add; get new weight in return; rows are sent to all
// procs except the current one and the master. At the end, the new weights
// are sent to the master for update.
// the slave *has* to collect the weights itself, since otherwise the
// master would be receiving a lot of suborders in a probably random order.
void
mpi_add_rows(sparse_mat_t *mat, int m, INT j, INT *ind)
{
    MPI_Status status;
    unsigned int buf[MPI_BUF_SIZE];
    int new_weight[MERGE_LEVEL_MAX];
    int mpi_rank, mpi_size, ibuf, cnt, k, nrecv, finished;
    char *done;
    
# if DEBUG >= 1
    fprintf(stderr, "Column %d is for processor %d (m=%d)\n",
	    j, mpi_get_proc_for_j(j), m);
# endif
    ibuf = 0;
    buf[ibuf++] = mpi_index;
    buf[ibuf++] = m;
    buf[ibuf++] = j;
    for(k = 0; k < m; k++){
	buf[ibuf++] = (unsigned)ind[k];
	buf[ibuf++] = -1; // TODO_MPI: some real length??
    }
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    done = (char *)malloc(mpi_size * sizeof(char));
    // get current proc rank
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    done[mpi_rank] = 2;
    for(k = 1; k < mpi_size; k++){
	if(k != mpi_rank){
	    MPI_Send(buf, ibuf, MPI_UNSIGNED, k, MPI_ADD_TAG, MPI_COMM_WORLD);
	    done[k] = 1;
	}
    }
    nrecv = 0;
    for(k = 0; k < m; k++)
	// this is the local weight
	new_weight[k] = (isRowNull(mat, ind[k]) ? 0 : lengthRow(mat, ind[k]));
    finished = 0;
    while(!finished){
	// we loop until we receive all new weights and then we exit
# if DEBUG >= 1
	fprintf(stdout, "[m_a_r] Proc %d waiting... nrecv=%d\n", 
		mpi_rank, nrecv);
# endif
	MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	switch(status.MPI_TAG){
	case MPI_SEND_W_TAG:
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
# if DEBUG >= 1
	    fprintf(stdout, "[m_a_r] source=%u:", status.MPI_SOURCE);
	    for(k = 0; k < cnt; k++)
		fprintf(stdout, " %u", buf[k]);
	    fprintf(stdout, "\n");
# endif
	    for(k = 5; k < cnt; k += 2)
		new_weight[(k>>1)-1] += (int)buf[k+1];
	    if(done[status.MPI_SOURCE] == 1){
		done[status.MPI_SOURCE] = 2;
		nrecv++;
		finished = (nrecv == (mpi_size-2));
	    }
	    else
		fprintf(stderr, "MPI??\n");
	    break;
	default:
	    // should not happen!!!!
	    mpi_err1("tag received: %d\n", status.MPI_TAG);
	}
    }
    new_weight[0] = 0; // this is the destroyed row
    ibuf = 0;
    buf[ibuf++] = mpi_index;
    buf[ibuf++] = m;
    buf[ibuf++] = j;
    for(k = 0; k < m; k++){
	buf[ibuf++] = ind[k];
	buf[ibuf++] = new_weight[k];
    }
    MPI_Send(buf, ibuf, MPI_UNSIGNED, 0, MPI_SEND_W_TAG, MPI_COMM_WORLD);
#if DEBUG >= 1
    mpi_err("Exiting from mpi_add_rows\n");
#endif
    free(done);
}

void
mpi_slave(int mpi_rank, sparse_mat_t *mat, FILE *purgedfile, char *purgedname,
	  char *outname)
{
    char *str;
    FILE *outfile;
    unsigned int buf[MPI_BUF_SIZE];
    MPI_Status status;
    double totopt = 0.0, totfill = 0.0, totMST = 0.0, totdel = 0.0;
    int jmin, jmax, cnt, i, k, m, verbose = 1, njrem = 0, i0;

    str = (char *)malloc(strlen(outname)+3);
    if(is_gzip(outname)){
	printf("Case of gz file NYI\n");
	exit(0);
    }
    else
	sprintf(str, "%s%03d", outname, mpi_rank);
    fprintf(stderr, "Slave outfile will be %s\n", str);
    outfile = gzip_open(str, "w");
    // loop forever
    while(1){
#if DEBUG >= 1
	mpi_err("Waiting...\n");
#endif
	MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
		 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#if DEBUG >= 1
	mpi_err2("Received %d from %d\n", status.MPI_TAG, status.MPI_SOURCE);
#endif
	// FIXME: sometimes, only 0 has the right to send things
	switch(status.MPI_TAG){
	case MPI_DIE_TAG:
	    mpi_err("I was asked to die...\n");
	    break;
	case MPI_MIN_M:
	    // buf = [index]
	    if(buf[0] >= mpi_index){
		mpi_index = buf[0];
#if DEBUG >= 1
		mpi_err1("New mpi_index = %d\n", mpi_index);
#endif
	    }
	    // buf <- [index, m_min]
	    m = minColWeight(mat);
	    buf[1] = m;
	    buf[2] = mpi_get_number_of_active_colums(mat, jmin, jmax);
	    MPI_Send(buf, 3, MPI_UNSIGNED, 0, MPI_MIN_M, MPI_COMM_WORLD);
	    break;
	case MPI_DO_M:
	    // buf = [index, m]
	    m = (int)buf[1];
	    // the local part
	    mpi_err1("doOneMerge(%d)\n", m);
	    doOneMerge(outfile, mat, &njrem, &totopt, &totfill, &totMST,
		       &totdel, m, verbose);
	    deleteEmptyColumns(mat);
	    fflush(outfile);
#if DEBUG >= 1
	    mpi_err2("nrows=%d ncols=%d\n", mat->rem_nrows, mat->rem_ncols);
#endif
	    MPI_Send(buf, 1, MPI_UNSIGNED, 0, MPI_M_DONE, MPI_COMM_WORLD);
	    break;
	case MPI_J_TAG:
	    jmin = (int)buf[0];
	    jmax = (int)buf[1];
	    mpi_err2("A slave does not need to do too much things...\n... but has to treat C[%d..%d[\n", jmin, jmax);
	    initWeightFromFile(mat, purgedfile, jmin, jmax);
	    gzip_close(purgedfile, purgedname);
	    fillSWAR(mat, jmin, jmax);
	    purgedfile = gzip_open(purgedname, "r");
	    readmat(mat, purgedfile, jmin, jmax);
	    gzip_close(purgedfile, purgedname);
	    break;
	case MPI_ADD_TAG:
	    // buf = [index, m, j, i1, l1, i2, l2, ..., ir, lr]
	    // TODO_MPI: receive length of r1+r2???? Perhaps pbs with MST
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
#if DEBUG >= 1
	    fprintf(stderr, "From proc=%d, I received %d values:",
		    status.MPI_SOURCE, cnt);
	    fprintf(stderr, " index=%d m=%d j=%d", buf[0], buf[1], buf[2]);
	    for(k = 3; k < cnt; k++)
		fprintf(stderr, " %d", (int)buf[k]);
	    fprintf(stderr, "\n");
#endif
	    // if buf[3] is non null, we must add it to all buf[k > 3]
	    i0 = (int)buf[3];
	    // humffffff
	    if(!isRowNull(mat, i0)){
		for(k = 5; k < cnt; k += 2){
		    i = (int)buf[k];
#if DEBUG >= 1
		    mpi_err2("R[%d] += R[%d]\n", i, i0);
#endif
		    if(isRowNull(mat, i)){
			fprintf(stderr, "Row %d is null -- %p!!!!!\n", 
				i, mat->rows[i]);
			// just copy
			mat->rows[i] = copyRow(mat->rows[3]);
			addRowSWAR(mat, i);
		    }
		    else
			addRowsSWAR(mat, i, i0, (int)buf[k+1]);
		}
		// i0 is no longer needed...!
		removeRowSWAR(mat, i0);
		destroyRow(mat, i0);
	    }
	    // now, send back the new weights of the rows (skipping 2...)
	    buf[4] = 0;
	    for(i = 5; i < cnt; i += 2){
		if(isRowNull(mat, (int)buf[i]))
		    buf[i+1] = 0;
		else
		    buf[i+1] = (unsigned)lengthRow(mat, (int)buf[i]);
#if DEBUG >= 1
		fprintf(stderr, "index=%u, m=%u j=%u", buf[0], buf[1], buf[2]);
		fprintf(stderr, " i=%u => w=%u\n", buf[i], buf[i+1]);
#endif
	    }
	    // back to the current sub-master
	    // buf = [index, m, j, i1, 0, i2, new_l2, ..., ir, newlr]
	    MPI_Send(buf, cnt, MPI_UNSIGNED, status.MPI_SOURCE,
		     MPI_SEND_W_TAG,MPI_COMM_WORLD);
	    break;
	case MPI_INACTIVATE_ROWS:
	    MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
	    mpi_inactivate_rows(outfile, mat, buf, cnt);
	    break;
	default:
	    mpi_err1("Unknown tag: %d\n", status.MPI_TAG);
	    break;
	}
	if(status.MPI_TAG == MPI_DIE_TAG)
	    break;
    }
    gzip_close(outfile, str);
}

void
mpi_start_slaves(int mpi_size, sparse_mat_t *mat)
{
    unsigned int buf[MPI_BUF_SIZE];
    int jmin, jmax, kappa;
    int rk;

    mpi_tab_j = (int *)malloc(mpi_size * sizeof(int));
#if 0
    kappa = mat->ncols/(mpi_size-1);
#endif
    for(rk = 1; rk < mpi_size; rk++){
#if 0
	// let kappa = ncols/(mpi_size-1)
	// slave_j for j < mpi_size-1 will have to treat [jmin..jmax[
	// jmin = (j-1)*kappa, jmax = j*kappa-1;
	// for j = mpi_size-1, this will be jmax = ncols
	jmin = (rk-1) * kappa;
	jmax = (rk == (mpi_size-1) ? mat->ncols : rk * kappa);
#else
	jmin = (rk == 1 ? 0 : jmax);
	jmax = jmin;
	kappa = 0;
	for(jmax = jmin; jmax < mat->ncols; jmax++){
	    kappa += abs(mat->wt[jmax]);
# if DEBUG >= 1
	    fprintf(stderr, "wt=%d kappa=%d\n", mat->wt[jmax], kappa);
# endif
	    if(((unsigned long)kappa) > (mat->weight / ((unsigned long)((mpi_size-1)))))
		break;
	}
#endif
	mpi_err3("Weight[%d..%d[ = %d\n", jmin, jmax, kappa);
	mpi_tab_j[rk-1] = jmin;
	buf[0] = (unsigned)jmin;
	buf[1] = (unsigned)jmax;
        MPI_Send(buf, 2, MPI_UNSIGNED, rk, MPI_J_TAG, MPI_COMM_WORLD);
    }
    mpi_tab_j[mpi_size-1] = jmax;
    fprintf(stderr, "# mpi_tab_j:");
    for(rk = 0; rk < mpi_size; rk++)
	fprintf(stderr, " %d", mpi_tab_j[rk]);
    fprintf(stderr, "\n");
}

int
mpi_get_minimal_m_proc(int *m, int *proc, int mpi_size, unsigned int index)
{
    MPI_Status status;
    unsigned int buf[MPI_BUF_SIZE];
    int k, nrecv, curr_ncols = 0;

    *m = MERGE_LEVEL_MAX, nrecv;
    buf[0] = index;
    for(k = 1; k < mpi_size; k++)
	MPI_Send(buf, 1, MPI_UNSIGNED, k, MPI_MIN_M, MPI_COMM_WORLD);
    *proc = -1;
    nrecv = 0;
    while(1){
	MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);
#if DEBUG >= 1
	fprintf(stderr, "source=%u: index=%u m_min=%u\n",
		status.MPI_SOURCE, buf[0], buf[1]);
#endif
	if(status.MPI_TAG == MPI_MIN_M){
	    if(buf[0] != index){
		fprintf(stderr, 
			"Should not happen: bad index %u instead of %u\n", 
			buf[0], index);
	    }
	    else
		nrecv++;
	    if(*m > (int)buf[1]){
		*proc = status.MPI_SOURCE;
		*m = (int)buf[1];
	    }
	    curr_ncols += (int)buf[2];
	    if(nrecv == (mpi_size-1))
		break;
	}
    }
    return curr_ncols;
}

void
mpi_master(FILE *outfile, sparse_mat_t *mat, int mpi_size)
{
    MPI_Status status;
    unsigned int buf[MPI_BUF_SIZE];
    int ibuf, m, proc, cnt, k, done, maxm = 0;
    int *row_weight, curr_ncols;

#if 0
    mat->mergelevelmax = 2;
#endif
    // let's initialize a new simpler array with the weights of the rows
    row_weight = (int *)malloc(mat->nrows * sizeof(int));
    for(k = 0; k < mat->nrows; k++)
	if(isRowNull(mat, k))
	    row_weight[k] = 0;
	else
	    row_weight[k] = lengthRow(mat, k);
    while(1){
	mpi_index++;
	curr_ncols = mpi_get_minimal_m_proc(&m, &proc, mpi_size, mpi_index);
	if((mpi_index == 1) || ((mpi_index % 10000) == 0))
	    fprintf(stderr, 
		    "Round %d (%2.2lf): maxm=%d nrows=%d ncols=%d[%d] (%d)\n",
		    mpi_index, seconds(), maxm, 
		    mat->rem_nrows, mat->rem_ncols, curr_ncols,
		    mat->rem_nrows - mat->rem_ncols);
	// look for minimal m for index
	//	mat->rem_ncols = curr_ncols;
#if DEBUG >= 1	
	mpi_err2("Minimal m is %d for proc=%d\n", m, proc);
#endif
	if(m > mat->mergelevelmax)
	    break;
	if(m > maxm)
	    maxm = m;
	// ask the proc to execute one step
	buf[0] = mpi_index;
	buf[1] = m;
	MPI_Send(buf, 2, MPI_UNSIGNED, proc, MPI_DO_M, MPI_COMM_WORLD);
	while(1){
#if DEBUG >= 1
	    mpi_err("Waiting...\n");
#endif
	    done = 0;
	    MPI_Recv(buf, MPI_BUF_SIZE, MPI_UNSIGNED, MPI_ANY_SOURCE,
		     MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	    switch(status.MPI_TAG){
	    case MPI_M_DONE:
#if DEBUG >= 1
		mpi_err("We stop here...\n");
#endif
		done = 1;
		break;
	    case MPI_SEND_W_TAG:
		// new rows received after additions
		MPI_Get_count(&status, MPI_UNSIGNED, &cnt);
		// buf = [index, m, j, i1, new_l1, i2, new_l2, ..., ir, new_lr]
		for(k = 3; k < cnt; k += 2){
#if DEBUG >= 1
		    mpi_err3("New row %d has now weight %d instead of %d\n",
			     buf[k], buf[k+1], row_weight[buf[k]]);
#endif
		    row_weight[buf[k]] = buf[k+1];
		    if(buf[k+1] == 0){
			mat->rem_nrows--;
			mat->rem_ncols--;
		    }
		}
		ibuf = 1;
		if((mat->rem_nrows - mat->rem_ncols) > mat->delta){
		    // look for too heavy rows
		    for(k = 3; k < cnt; k += 2){
			if(row_weight[buf[k]] > mat->rwmax){
#if DEBUG >= 1
			    mpi_err3("Row %d is too heavy (%d/%d)\n",
				     buf[k], row_weight[buf[k]], mat->rwmax);
#endif
			    buf[ibuf++] = buf[k]; // ohhhhhhhhhh!
			    row_weight[buf[k]] = 0;
			    report1(outfile, (int)buf[k]);
			    mat->rem_nrows--;
			}
		    }
		}
		if(ibuf > 1)
		    for(k = 1; k < mpi_size; k++){
			mpi_index++;
			buf[0] = mpi_index;
			MPI_Send(buf, ibuf, MPI_UNSIGNED, k, 
				 MPI_INACTIVATE_ROWS, MPI_COMM_WORLD);
		    }
		break;
	    default:
		fprintf(stderr, "Received tag was %d\n", status.MPI_TAG);
	    }
	    if(done)
		break;
	}
    }
    fprintf(stderr, "Finally (%d -- %2.2lf): nrows=%d ncols=%d\n",
		    mpi_index, seconds(), mat->rem_nrows, mat->rem_ncols);
}

void
mpi_start_master(char *outname, int mpi_size, sparse_mat_t *mat)
{
    FILE *outfile;
    char *str;

    str = (char *)malloc(strlen(outname)+3);
    if(is_gzip(outname)){
	printf("Case of gz file NYI\n");
	exit(0);
    }
    else
	sprintf(str, "%s%03d", outname, 0);
    fprintf(stderr, "Master outfile will be %s\n", str);
    mpi_index = 0;
    outfile = gzip_open(str, "w");
    report2(outfile, mat->nrows, mat->ncols);
    mpi_start_slaves(mpi_size, mat);
    mpi_master(outfile, mat, mpi_size);
    gzip_close(outfile, str);
    free(str);
}
#endif

int
main(int argc, char *argv[])
{
    FILE *purgedfile, *outfile;
    sparse_mat_t mat;
    char *purgedname = NULL, *hisname = NULL, *outname = NULL;
    int nrows, ncols;
    int cwmax = 20, rwmax = 1000000, maxlevel = 2, iprune = 0, keep = 128;
    int verbose = 0; /* default verbose level */
    double tt;
    double kprune = 1.0; /* prune keeps kprune * (initial excess) */
    double ratio = 1.1; /* bound on cN_new/cN to stop the computation */
    int i, forbw = 0, coverNmax = 0;
#if USE_MPI
    int mpi_rank, mpi_size; // 0 <= mpi_rank < mpi_size
#endif
    
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    fprintf(stderr, "Hello, world, I am %d of %d\n", mpi_rank, mpi_size);
    if(mpi_rank == 0){
	// I don't want this to appear many times
#endif
    fprintf (stderr, "%s.r%s", argv[0], REV);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");
#ifdef USE_MPI
    }
#endif

#if TEX
    fprintf(stderr, "\\begin{verbatim}\n");
#endif    
    while(argc > 1 && argv[1][0] == '-'){
        if (argc > 2 && strcmp (argv[1], "-mat") == 0){
	    purgedname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-rebuild") == 0){
	    hisname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-cwmax") == 0){
	    cwmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-rwmax") == 0){
	    rwmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-maxlevel") == 0){
	    maxlevel = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-prune") == 0){
	    kprune = strtod(argv[2], NULL);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-keep") == 0){
	    keep = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 1 && strcmp (argv[1], "-v") == 0){
            verbose ++;
	    argc -= 1;
	    argv += 1;
	}
	else if (argc > 2 && strcmp (argv[1], "-forbw") == 0){
	    forbw = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-ratio") == 0){
	    ratio = strtod(argv[2], NULL);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-coverNmax") == 0){
	    coverNmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-out") == 0){
	    outname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else 
	    break;
    }
#ifdef USE_MPI    
    if(mpi_rank != 0){
	char mpitmp[1024];
	sprintf(mpitmp, "%s%03d.err", outname, mpi_rank);
	freopen(mpitmp, "w", stderr);
    }
#endif

    purgedfile = gzip_open(purgedname, "r");
    ASSERT_ALWAYS(purgedfile != NULL);
    fscanf(purgedfile, "%d %d", &nrows, &ncols);

    mat.nrows = nrows;
    mat.ncols = ncols;
    mat.delta = keep; // FIXME: change the name of delta
    
    tt = seconds();
    initMat(&mat, cwmax, rwmax, maxlevel);
    fprintf(stderr, "Time for initMat: %2.2lf\n", seconds()-tt);

#ifdef USE_MPI
    if(mpi_rank != 0){
	mpi_slave(mpi_rank, &mat, purgedfile, purgedname, outname);
	// TODO: clean the mat data structure for a slave also
    }
    else{ // for the master only...!
	// TODO: actually, the master does not need to do other things
	// like reading a piece of the matrix...
#endif    
    tt = seconds();
    initWeightFromFile(&mat, purgedfile, 0, mat.ncols);
    fprintf(stderr, "Time for initWeightFromFile: %2.2lf\n", seconds()-tt);
    gzip_close(purgedfile, purgedname);
    
    tt = seconds();
    fillSWAR(&mat, 0, mat.ncols);
    fprintf(stderr, "Time for fillSWAR: %2.2lf\n", seconds()-tt);
#ifdef USE_MPI
    }
#endif
    
#ifdef USE_MPI
    if(mpi_rank == 0){
#endif
    purgedfile = gzip_open(purgedname, "r");
    ASSERT_ALWAYS(purgedfile != NULL);
    readmat(&mat, purgedfile, 0, mat.ncols);
    gzip_close(purgedfile, purgedname);
#if DEBUG >= 3
    checkmat(&mat);
#endif

#ifdef USE_MPI
    mpi_start_master(outname, mpi_size, &mat);
    mpi_kill_slaves(mpi_size);
    MPI_Finalize();
    return 0;
#endif

    outfile = gzip_open(outname, "w");
    report2(outfile, mat.nrows, mat.ncols);

    /* iprune is the excess we want at the end of prune */
    iprune = (mat.nrows-mat.ncols) * kprune;
    if (iprune < mat.delta) /* ensures iprune >= DELTA */
      iprune = mat.delta;
    /* only call prune if the current excess is larger than iprune */
    if (iprune < mat.nrows - mat.ncols){
	double tt = seconds();
	prune(outfile, &mat, iprune);
	fprintf(stderr, "Pruning: nrows=%d ncols=%d %2.2lf\n",
		mat.rem_nrows, mat.rem_ncols, seconds()-tt);
    }

    // ouhhhhhhhhhh
    // we do not have a clear idea of which function to minimize...!
#if 0
    forbw = 0;
    fprintf(stderr, "WARNING: forcing forbw=0...!!!!\n");
#endif

#if M_STRATEGY <= 2
    merge(&mat, maxlevel, verbose, forbw);
#else
    mergeOneByOne(outfile, &mat, maxlevel, verbose, forbw, ratio, coverNmax);
#endif
    gzip_close(outfile, outname);
    fprintf(stderr, "Final matrix has w(M)=%lu, ncols*w(M)=%lu\n",
	    mat.weight, ((unsigned long)mat.rem_ncols) * mat.weight);
#if TEX
    fprintf(stderr, "\\end{verbatim}\n");
#endif
#if DEBUG >= 1
    {
	FILE *ofile = fopen("toto", "w");
	dumpSparse(ofile, &mat);
	fclose(ofile);
    }
#endif
    closeSWAR(/*&mat*/);
#ifdef USE_MPI
    mpi_kill_slaves(mpi_size);
    }
    MPI_Finalize();
#endif
    return 0;
}
