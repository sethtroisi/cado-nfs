// TODO: use unified compact lists...!
// TODO: reintroduce mat->weight

/* If one wants to change the R data structure, please check the diff of
   revision 1568, which clearly identifies the places where R is used. */

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
#include "dclist.h"
#include "sparse_mat.h"
#include "report.h"
#include "swar.h"

#ifdef USE_MARKOWITZ
#include "markowitz.h"
#endif
#include "merge_mono.h"
#include "prune.h"

#define STR_LEN_MAX 1024

#define DEBUG 0
#define TEX 0

#define USE_CONNECT 0

typedef struct {
  int len;    /* number of non-zero entries */
  INT *val;   /* array of size len containing column indices of non-zero
                 entries */
} rel_t;

int 
cmp(const void *p, const void *q) {
    int x = *((int *)p);
    int y = *((int *)q);
    return (x <= y ? -1 : 1);
}

// TODO: ncols could be a new renumbered ncols...
void
initMat(sparse_mat_t *mat, INT jmin, INT jmax)
{
#if USE_MERGE_FAST
    // mat->cwmax+2 to prevent bangs
    dclist *S = (dclist *)malloc((mat->cwmax+2) * sizeof(dclist));
    dclist *A = (dclist *)malloc((jmax-jmin) * sizeof(dclist));
    INT **R;
    int k;

    ASSERT_ALWAYS(S != NULL);
    ASSERT_ALWAYS(A != NULL);
    // S[0] has a meaning at least for temporary reasons
    for(k = 0; k <= mat->cwmax+1; k++)
	S[k] = dclistCreate(-1);
    mat->S = S;
    mat->A = A;
#endif

    mat->jmin = jmin;
    mat->jmax = jmax;
#if USE_TAB == 0
    mat->data = (rel_t*) malloc (mat->nrows * sizeof (rel_t));
    ASSERT_ALWAYS(mat->data != NULL);
#else
    mat->rows = (INT **)malloc(mat->nrows * sizeof (INT *));
    ASSERT_ALWAYS(mat->rows != NULL);
#endif
    mat->wt = (int*) malloc ((mat->jmax-mat->jmin) * sizeof (int));
    memset(mat->wt, 0, (mat->jmax-mat->jmin) * sizeof (int));
#if USE_MERGE_FAST
    R = (INT **)malloc((mat->jmax-mat->jmin) * sizeof(INT *));
    ASSERT_ALWAYS(R != NULL);
    mat->R = R;
#endif
}

void
checkData(sparse_mat_t *mat)
{
    int j, nbj = 0;

    for(j = 0; j < mat->ncols; j++)
	if(mat->wt[GETJ(mat, j)])
	    nbj++;
    if(mat->rem_ncols != nbj){
	fprintf(stderr, "rem_ncols=%d nbj=%d\n", mat->rem_ncols, nbj);
	exit(1);
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
	fprintf(stderr, "& %d", mat->wt[GETJ(mat, j)]);
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
	fprintf(stderr, "& %d", mat->wt[GETJ(mat, j)]);
    fprintf(stderr, "\\\\\\hline\n\\end{array}$$\n");
#ifndef USE_COMPACT_R
    fprintf(stderr, "\\begin{verbatim}\n");
    for(j = 0; j < mat->ncols; j++){
	if(mat->R[GETJ(mat, j)] == NULL)
	    continue;
	fprintf(stderr, "  R[%d]=", j);
	for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
          fprintf(stderr, " %ld", (long int) mat->R[GETJ(mat, j)][k]);
	fprintf(stderr, "\n");
    }
#endif
}

void
printStats(sparse_mat_t *mat)
{
    int w;

    for(w = 0; w <= mat->cwmax; w++)
	fprintf(stderr, "I found %d primes of weight %d\n",
		dclistLength(mat->S[w]->next), w);
}

/* compute the weight mat->wt[j] of each column j, for the set of relations
   in file purgedfile.
*/
void
initWeightFromFile(sparse_mat_t *mat, FILE *purgedfile)
{
    INT j, jmin = mat->jmin, jmax = mat->jmax;
    int i, k, ret, nc;

    memset(mat->wt, 0, (mat->jmax - mat->jmin) * sizeof(int));
    for(i = 0; i < mat->nrows; i++){
	ret = fscanf(purgedfile, "%d", &k); // unused index to rels file
	ASSERT_ALWAYS (ret == 1);
	ret = fscanf (purgedfile, "%d", &nc);
	ASSERT_ALWAYS (ret == 1);
	for(k = 0; k < nc; k++){
	    ret = fscanf(purgedfile, PURGE_INT_FORMAT, &j);
	    ASSERT_ALWAYS (ret == 1);
	    if((j >= jmin) && (j < jmax))
		mat->wt[GETJ(mat, j)]++;
	}
    }
}

#define BUF_LEN 100

/* Reads a matrix file, and puts in mat->wt[j], 0 <= j < ncols, the
   weight of column j (adapted from matsort.c).

   We skip columns that are too heavy.

   mat->wt is already filled and put to - values when needed. So don't touch.
*/
int
readmat(sparse_mat_t *mat, FILE *file)
{
    int ret;
    int i, k;
    int nc, buf[BUF_LEN];
#if (USE_MERGE_FAST < 3) && !defined(USE_MPI)
    int nh = 0;
#endif
    INT ibuf, j, jmin = mat->jmin, jmax = mat->jmax;

    ret = fscanf (file, "%d %d", &(mat->nrows), &(mat->ncols));
    ASSERT_ALWAYS (ret == 2);
    
    fprintf(stderr, "Reading matrix of %d rows and %d columns: excess is %d\n",
	    mat->nrows, mat->ncols, mat->nrows - mat->ncols);
    mat->rem_nrows = mat->nrows;
    mat->weight = 0;

    for (i = 0; i < mat->nrows; i++){
	if(!(i % 100000))
	    fprintf(stderr, "Reading %d-th row\n", i);
	ret = fscanf(file, "%d", &nc); // unused index to rels file
	if(ret != 1){
	    fprintf(stderr, "Pb1 in readmat\n");
	    return 0;
	}
	ret = fscanf (file, "%d", &nc);
	if(ret != 1){
	    fprintf(stderr, "Pb2 in readmat\n");
	    return 0;
	}
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
#if (USE_MERGE_FAST < 3) && !defined(USE_MPI)
	    int nb_heavy_j = 0;
#endif
	    for(k = 0, ibuf = 0; k < nc; k++){
		ret = fscanf(file, PURGE_INT_FORMAT, &j);
		if(ret != 1){
		    fprintf(stderr, "Pb3 in readmat: k=%d\n", k);
		    return 0;
		}
		ASSERT_ALWAYS (0 <= j && j < mat->ncols);
#if (USE_MERGE_FAST < 3) && !defined(USE_MPI)
		if(mat->wt[GETJ(mat, j)] < 0)
		    nb_heavy_j++;
#endif
#if USE_MERGE_FAST <= 1
		if(mat->wt[GETJ(mat, j)] > 0)
		    buf[ibuf++] = j;
#else
		// always store j in the right interval...!
		if((j >= jmin) && (j < jmax)){
		    buf[ibuf++] = j;
		    // this will be the weight in the current slice
		    mat->weight++; 
		    if(mat->wt[GETJ(mat, j)] > 0){ // redundant test?
#ifndef USE_COMPACT_R
			mat->R[GETJ(mat, j)][0]++;
			mat->R[GETJ(mat, j)][mat->R[GETJ(mat, j)][0]] = i;
#else
			fprintf(stderr, "R: NYI in readmat\n");
			exit(1);
#endif
		    }
# if DEBUG >= 1
		    if(j == 15054){
			fprintf(stderr, "i=%d, j=%d (ind/nc=%d) mat->%d\n", 
				i, j, k, mat->R[GETJ(mat, j)][0]);
			fflush(stderr);
		    }
# endif
		}
#endif
	    }
#if (USE_MERGE_FAST < 3) && !defined(USE_MPI)
	    if(nb_heavy_j == nc){
		// all the columns are heavy and thus will never participate
		mat->rows[i] = NULL;
		nh++;
		continue;
	    }
#endif
	    ASSERT_ALWAYS(ibuf <= BUF_LEN);
	    // TODO: do not store rows not having at least one light
	    // column, but do not decrease mat->nrows!
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
#if (USE_MERGE_FAST < 3) && !defined(USE_MPI)
    fprintf(stderr, "Number of heavy rows: %d\n", nh);
#endif
    // we need to keep informed of what really happens; this will be an upper
    // bound on the number of active columns, I guess
    mat->rem_ncols = mat->ncols;
    printStats(mat);
    return 1;
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
	mat->wt[GETJ(mat, cell(mat, i, k))]--;
#else
    INT k;

    for(k = 1; k <= mat->rows[i][0]; k++)
	mat->wt[GETJ(mat, mat->rows[i][k])]--;
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
	mat->wt[GETJ(mat, cell(mat, i, k))]++;
#else
    INT k;

    for(k = 1; k <= mat->rows[i][0]; k++)
	mat->wt[GETJ(mat, mat->rows[i][k])]++;
#endif
    mat->weight += lengthRow(mat, i);
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

// what is the weight of the sum of Ra and Rb? Works even in the partial
// scenario of MPI.
int
weightSum(sparse_mat_t *mat, int i1, int i2)
{
    int k1, k2, w = 0, len1, len2;

    len1 = (isRowNull(mat, i1) ? 0 : lengthRow(mat, i1));
    len2 = (isRowNull(mat, i2) ? 0 : lengthRow(mat, i2));
    if((len1 == 0) || (len2 == 0))
	fprintf(stderr, "i1=%d i2=%d len1=%d len2=%d\n", i1, i2, len1, len2);
#if USE_TAB == 0
    k1 = k2 = 0;
    while((k1 < len1) && (k2 < len2)){
#else
    k1 = k2 = 1;
    while((k1 <= len1) && (k2 <= len2)){
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
    w += (k1 > len1 ? 0 : len1-k1+1);
    w += (k2 > len2 ? 0 : len2-k2+1);
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
merge_m_slow(report_t *rep, sparse_mat_t *mat, int m, int verbose)
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
	if(mat->wt[GETJ(mat, j)] != m)
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
	reportn(rep, ind, m);
	totadd += (cputime()-tt);
	removeWeightFromRow(mat, ind[0]);
	destroyRow(mat, ind[0]);
	mat->rem_nrows--;
	mat->rem_ncols--;
	mat->wt[GETJ(mat, j)] = 0;
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

#ifndef USE_COMPACT_R
    for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
	if(mat->R[GETJ(mat, j)][k] != -1)
	    nchk++;
    ASSERT(nchk == (mat->wt[GETJ(mat, j)] >= 0 ? mat->wt[GETJ(mat, j)] : -mat->wt[GETJ(mat, j)]));
    if(m != -1){
	if(nchk != m){
	    fprintf(stderr, "HYPERCHECK:");
	    fprintf(stderr, "mat->R[%d][0]=%ld, m=%d\n", j,
                    (long int) mat->R[GETJ(mat, j)][0], m);
	    fprintf(stderr, "Gasp: nchk=%d\n", nchk);
	}
	ASSERT(nchk == m);
    }
#else
    fprintf(stderr, "R: NYI in checkCoherence\n");
    exit(1);
#endif
}

void
destroyRj(sparse_mat_t *mat, int j)
{
#ifndef USE_COMPACT_R
    free(mat->R[GETJ(mat, j)]);
    mat->R[GETJ(mat, j)] = NULL;
#else
    fprintf(stderr, "R: NYI in destroyRj\n");
    exit(1);
#endif
}

// Don't touch to R[j][0]!!!!
void
remove_i_from_Rj(sparse_mat_t *mat, int i, int j)
{
    // be dumb for a while
    int k;

#ifndef USE_COMPACT_R
    if(mat->R[GETJ(mat, j)] == NULL){
	fprintf(stderr, "Row %d already empty\n", j);
	return;
    }
    for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
	if(mat->R[GETJ(mat, j)][k] == i){
#if DEBUG >= 1
	    fprintf(stderr, "Removing row %d from R[%d]\n", i, j);
#endif
	    mat->R[GETJ(mat, j)][k] = -1;
	    break;
	}
#else
    fprintf(stderr, "R: NYI in remove_i_from_Rj\n");
    exit(1);
#endif
}

// cell M[i, j] is incorporated in the data structure. It is used
// later on in cases where i is not already present in row[j].
void
add_i_to_Rj(sparse_mat_t *mat, int i, int j)
{
    // be semi-dumb for a while
    int k;
    
#if DEBUG >= 1
    fprintf(stderr, "Adding row %d to R[%d]\n", i, j);
#endif
#ifndef USE_COMPACT_R
    for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
	if(mat->R[GETJ(mat, j)][k] == -1)
	    break;
    if(k <= mat->R[GETJ(mat, j)][0]){
	// we have found a place where it is -1
	mat->R[GETJ(mat, j)][k] = i;
    }
    else{
#if DEBUG >= 1
	fprintf(stderr, "WARNING: reallocing things in add_i_to_Rj for R[%d]\n", j);
#endif
	int l = mat->R[GETJ(mat, j)][0]+2;
	mat->R[GETJ(mat, j)] = (INT *)realloc(mat->R[GETJ(mat, j)], l * sizeof(INT));
	mat->R[GETJ(mat, j)][l-1] = i;
	mat->R[GETJ(mat, j)][0] = l-1;
    }
#else
    fprintf(stderr, "R: NYI in add_i_to_Rj\n");
    exit(1);
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

//////////////////////////////////////////////////////////////////////
// making things independent of the real data structure used
//////////////////////////////////////////////////////////////////////

// The cell [i, j] must be incorporated to the data structure, at least
// if j is not too heavy, etc. 
// TODO: this should be shared and redistributed
// between swar.c and the present file.
void
addCellAndUpdate(sparse_mat_t *mat, int i, INT j)
{
#if TRACE_ROW >= 0
    if(i == TRACE_ROW)
	fprintf(stderr, "TRACE_ROW: addCellSWAR i=%d j=%d\n", i, j);
#endif
    if(addColSWAR(mat, j) >= 0)
	// update R[j] by adding i
	add_i_to_Rj(mat, i, j);
}

void
removeCellAndUpdate(sparse_mat_t *mat, int i, INT j)
{
    removeCellSWAR(mat, i, j);
}

void
removeColumnAndUpdate(sparse_mat_t *mat, int j)
{
    remove_j_from_SWAR(mat, j);
}

// These columns are simply removed from the current squeleton matrix, but
// not from the small matrix.
int
deleteAllColsFromStack(report_t *rep, sparse_mat_t *mat, int iS)
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
#ifndef USE_COMPACT_R
	    for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
		if(mat->R[GETJ(mat, j)][k] != -1){
# if TRACE_COL >= 0
		    if(j == TRACE_COL)
			fprintf(stderr, "deleteAllCols: row is %d\n",mat->R[GETJ(mat, j)][k]);
# endif
		    remove_j_from_row(mat, mat->R[GETJ(mat, j)][k], j);
		    removeRowDefinitely(rep, mat, mat->R[GETJ(mat, j)][k]);
		    mat->rem_ncols--;
		}
	    mat->wt[GETJ(mat, j)] = 0;
#else
	    fprintf(stderr, "R: NYI in deleteAllColsFromStack\n");
	    exit(1);
#endif
	}
	k = mat->wt[GETJ(mat, j)]; // make a copy of the weight
	removeColumnAndUpdate(mat, j);
	// mat->wt[j] was put to 0...
	mat->wt[GETJ(mat, j)] = -k; // restore and update
#endif
    }
#if DEBUG >= 1
    if(njrem > 0)
	fprintf(stderr, "deleteAllColsFromStack[%d]: %d\n", iS, njrem);
#endif
    return njrem;
}

//////////////////////////////////////////////////////////////////////
// now, these guys are generic...!

// for all j in row[i], removes j and update data
void
removeRowAndUpdate(sparse_mat_t *mat, int i)
{
    int k;

#if TRACE_ROW >= 0
    if(i == TRACE_ROW)
	fprintf(stderr, "TRACE_ROW: removeRowAndUpdate i=%d\n", i);
#endif
    mat->weight -= lengthRow(mat, i);
    for(k = 1; k <= lengthRow(mat, i); k++){
#if TRACE_COL >= 0
	if(cell(mat, i, k) == TRACE_COL){
	    fprintf(stderr, "removeRowAndUpdate removes %d from R_%d\n",
		    TRACE_COL, i);
	}
#endif
	removeCellAndUpdate(mat, i, cell(mat, i, k));
    }
}

// All entries M[i, j] are added to the structure.
void
addOneRowAndUpdate(sparse_mat_t *mat, int i)
{
    int k;

    mat->weight += lengthRow(mat, i);
    for(k = 1; k <= lengthRow(mat, i); k++)
	addCellAndUpdate(mat, i, cell(mat, i, k));
}

// realize mat[i1] += mat[i2] and update the data structure.
// len could be the real length of row[i1]+row[i2] or -1.
void
addRowsAndUpdate(sparse_mat_t *mat, int i1, int i2, int len)
{
    // cleaner one, that shares addRowsData() to prepare the next move...!
    // i1 is to disappear, replaced by a new one
    removeRowAndUpdate(mat, i1);
    // we know the length of row[i1]+row[i2]
    addRows(mat->rows, i1, i2, len);
    addOneRowAndUpdate(mat, i1);
}

int
removeSingletons(report_t *rep, sparse_mat_t *mat)
{
#if USE_MERGE_FAST == 0
    INT i, j;

    for(j = 0; j < mat->ncols; j++)
	if(mat->wt[GETJ(mat, j)] == 1){
	    // find row...
	    for(i = 0; i < mat->nrows; i++)
		if(hasCol(mat->rows, i, j))
		    break;
	    ASSERT(i < mat->nrows);
	    removeWeightFromRow(mat, i);
	    destroyRow(mat, i);
	    report1(rep, i); // signal to replay...!
	}
#else
    return deleteAllColsFromStack(rep, mat, 1);
#endif
}

int
deleteHeavyColumns(report_t *rep, sparse_mat_t *mat)
{
    return deleteAllColsFromStack(rep, mat, mat->cwmax+1);
}

void
removeRowDefinitely(report_t *rep, sparse_mat_t *mat, INT i)
{
    removeRowAndUpdate(mat, i);
    destroyRow(mat, i);
    report1(rep, i);
    mat->rem_nrows--;
}

// try all combinations to find the smaller one
void
tryAllCombinations(report_t *rep, sparse_mat_t *mat, int m, INT *ind)
{
    int i, k;

    i = findBestIndex(mat, m, ind);
#if DEBUG >= 1
    fprintf(stderr, "Minimal is i=%d (%d %d)\n", i, ind[0], ind[1]);
#endif
    for(k = 0; k < m; k++)
	if(k != i){
	    addRowsAndUpdate(mat, ind[k], ind[i], -1);
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
    reportn(rep, ind, m);
    removeRowAndUpdate(mat, ind[0]);
    destroyRow(mat, ind[0]);
}

// add u to its sons; we save the history in the history array, so that
// we can report all at once and prepare for MPI.
int
addFatherToSonsRec(int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		   sparse_mat_t *mat, int m, int *ind,
		   int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX],
		   int *father,
		   int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		   int u, int level0)
{
    int k, itab = 1, i1, i2, level = level0;

    if(sons[u][0] == 0)
	// nothing to do for a leaf...!
	return -1;
    i2 = ind[u];
    if(u == father[u])
	history[level0][itab++] = i2; // trick for the root!!!
    else
	// the usual trick not to destroy row
	history[level0][itab++] = -(i2+1);
    for(k = 1; k <= sons[u][0]; k++){
	i1 = addFatherToSonsRec(history, mat, m, ind, A, 
				father, sons, sons[u][k], level+1);
	if(i1 != -1)
	    level = i1;
	i1 = ind[sons[u][k]];
	// add u to its son
	addRowsAndUpdate(mat, i1, i2, A[sons[u][k]][u]);
	history[level0][itab++] = i1;
    }
    history[level0][0] = itab-1;
    return level;
}

int
addFatherToSons(int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		sparse_mat_t *mat, int m, int *ind,
		int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX],
		int *father,
                int *height MAYBE_UNUSED, int hmax MAYBE_UNUSED,
		int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1])
{
    return addFatherToSonsRec(history, mat, m, ind, A, father, sons, 0, 0);
}

void
MSTWithA(report_t *rep, sparse_mat_t *mat, int m, INT *ind, double *tMST,
	 int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
    int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];
    int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];
    int father[MERGE_LEVEL_MAX], height[MERGE_LEVEL_MAX], hmax, i;

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
    hmax = addFatherToSons(history, mat, m, ind, A, father, height, hmax,sons);
    for(i = hmax; i >= 0; i--)
#if 0
	reporthis(rep, history, i);
#else
        reportn(rep, history[i]+1, history[i][0]);
#endif
    removeRowAndUpdate(mat, ind[0]);
    destroyRow(mat, ind[0]);
}

void
useMinimalSpanningTree(report_t *rep, sparse_mat_t *mat, int m,
		       INT *ind, double *tfill, double *tMST)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX];

    *tfill = seconds();
    fillRowAddMatrix(A, mat, m, ind);
    *tfill = seconds()-*tfill;
    MSTWithA(rep, mat, m, ind, tMST, A);
}

void
findOptimalCombination(report_t *rep, sparse_mat_t *mat, int m,
		       INT *ind, double *tfill, double *tMST, int useMST)
{
    if((m <= 2) || !useMST){
	*tfill = *tMST = 0;
	tryAllCombinations(rep, mat, m, ind);
    }
    else
	useMinimalSpanningTree(rep, mat, m, ind, tfill, tMST);
}

void
checkWeight(sparse_mat_t *mat, INT j)
{
    int i, w = 0;

    fprintf(stderr, "Rows containing %ld:", (long int) j);
    for(i = 0; i < mat->nrows; i++)
	if(!isRowNull(mat, i))
	    if(hasCol(mat->rows, i, j)){
		fprintf(stderr, " %d", i);
		w++;
	    }
    fprintf(stderr, "\n");
    ASSERT(w == (mat->wt[GETJ(mat, j)] >= 0 ? mat->wt[GETJ(mat, j)] : -mat->wt[GETJ(mat, j)]));
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

// j has weight m, which should coherent with mat->wt[j] == m
void
mergeForColumn(report_t *rep, double *tt, double *tfill, double *tMST,
	       sparse_mat_t *mat, int m, INT j, int useMST)
{
    INT ind[MERGE_LEVEL_MAX];
    int ni, k;

#ifdef USE_MARKOWITZ
    // let's be cautious...
    if(mat->wt[j] != m){
	fprintf(stderr, "GASP: wt[%d]=%d != %d\n", j, mat->wt[j], m);
    }
#endif

#if DEBUG >= 1
    fprintf(stderr, "Treating column %d of weight %d\n", j, mat->wt[GETJ(mat, j)]);
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
#ifndef USE_COMPACT_R
    for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++){
	if(mat->R[GETJ(mat, j)][k] != -1){
	    ind[ni++] = mat->R[GETJ(mat, j)][k];
	    if(ni == m)
		break;
	}
    }
#else
    fprintf(stderr, "R: NYI in mergeForColumn\n");
    exit(1);
#endif
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
    findOptimalCombination(rep, mat, m, ind, tfill, tMST, useMST);
    *tt = seconds()-(*tt);
    mat->rem_nrows--;
    mat->rem_ncols--;
    removeColumnAndUpdate(mat, j);
}

// maxdo is 0 if we want to perform a non-bounded number of operations; 
// an integer >= 1 otherwise.
// Default for the monoproc version is 0.
int
merge_m_fast(report_t *rep, sparse_mat_t *mat, int m, int maxdo, int verbose)
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
	mergeForColumn(rep, &tt, &tfill, &tMST, mat, m, j, 1);
#if TRACE_COL >= 0
	fprintf(stderr, "wt[%d]=%d after findOptimalCombination\n",
		TRACE_COL, mat->wt[GETJ(mat, TRACE_COL)]);
#endif
	totopt += tt;
	totfill += tfill;
	totMST += tMST;
	tt = seconds();
	njrem += deleteHeavyColumns(rep, mat);
	totdel += (seconds()-tt);
#if TRACE_COL >= 0
        fprintf(stderr, "wt[%d]=%d after deleteHeavyColumns\n",
                TRACE_COL, mat->wt[GETJ(mat, TRACE_COL)]);
#endif
	maxdo--;
	if(maxdo == 0)
	    // what a trick: maxdo = 0 => never == 0 until infinity
	    return njrem;
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
merge_m(report_t *rep, sparse_mat_t *mat, int m, int maxdo, int verbose)
{
#if USE_MERGE_FAST >= 1
    return merge_m_fast(rep, mat, m, maxdo, verbose);
#else
  return merge_m_slow(rep, mat, m, verbose);
#endif
}

// TODO: use mergemax?
int
mergeGe2(report_t *rep, sparse_mat_t *mat, int m, int maxdo, int verbose)
{
#if USE_MERGE_FAST == 0
    int j, nbm;

    // TODO: remove this and start directly merge_m?
    for(j = 0, nbm = 0; j < mat->ncols; j++)
	if(mat->wt[GETJ(mat, j)] == m){
# if DEBUG >= 1
	    fprintf(stderr, "# wt[%d] = %d\n", j, m);
# endif
	    nbm++;
	}
    fprintf(stderr, "There are %d column(s) of weight %d\n", nbm, m);
    if(nbm)
#endif
	return merge_m(rep, mat, m, maxdo, verbose);
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
	if((mat->wt[GETJ(mat, j)] > 0) && (mat->wt[GETJ(mat, j)] < minw))
	    minw = mat->wt[GETJ(mat, j)];
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

#if (USE_MERGE_FAST < 3) && !defined(USE_MPI)
int
isRowToBeDeleted(sparse_mat_t *mat, int i)
{
    int heavy = 0;

    int k;

    heavy = 1;
    for(k = 1; k <= lengthRow(mat, i); k++)
	if((cell(mat,i,k) != -1) && (mat->wt[GETJ(mat, cell(mat,i,k))] >= 0)){
	    heavy = 0;
	    break;
	}
    return heavy;
}
#endif

int
inspectRowWeight(report_t *rep, sparse_mat_t *mat)
{
    int i, nirem = 0, niremmax = 128;
#if (USE_MERGE_FAST < 3) && !defined(USE_MPI)
    int useless = 0;
#endif

    for(i = 0; i < mat->nrows; i++){
	if((mat->rem_nrows - mat->rem_ncols) <= mat->delta)
	    return nirem;
	if(!isRowNull(mat, i)){
	    if(lengthRow(mat, i) > mat->rwmax){
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
		removeRowDefinitely(rep, mat, i);
		nirem++;
		if(!(nirem % 10000))
		    fprintf(stderr, "#removed_rows=%d at %2.2lf\n",
			    nirem, seconds());
		if(nirem > niremmax)
		    break;
	    }
	}
    }
#if (USE_MERGE_FAST < 3) && !defined(USE_MPI)
    // only activated rarely...!
    for(i = 0; i < mat->nrows; i++)
	if(!isRowNull(mat, i) && isRowToBeDeleted(mat, i)){
# if DEBUG >= 1
	    fprintf(stderr, "Row %d is no longer useful in merge\n", i);
# endif
	    removeRowAndUpdate(mat, i);
	    destroyRow(mat, i);
	    useless++;
	}
    fprintf(stderr, "#useless rows=%d\n", useless);
#endif
    return nirem;
}

// It could be fun to find rows s.t. at least one j is of small weight.
// Or components of rows with some sharing?
// We could define sf(row) = sum_{j in row} w(j) and remove the
// rows of smallest value of sf?
// For the time being, remove heaviest rows.
int
deleteScore(sparse_mat_t *mat, INT i)
{
#if 1
    // plain weight to remove heaviest rows
    return lengthRow(mat, i);
#endif
#if 0
    // -plain weight to remove lightest rows
    return -lengthRow(mat, i);
#endif
#if 0
    // not using rows with too many heavy column part
    int k, s = 0;

    for(k = 1; k <= lengthRow(mat, i); k++)
	s += abs(mat->wt[GETJ(mat, mat->rows[i][k])]);
    return s;
#endif
#if 0
    // not using rows with too many light columns
    int k, s = 0;

    for(k = 1; k <= lengthRow(mat, i); k++)
	s += abs(mat->wt[GETJ(mat, mat->rows[i][k])]);
    return -s;
#endif
#if 0
    // return the weight of the lightest column
    int k, s = abs(mat->wt[GETJ(mat, mat->rows[i][1])]);

    for(k = 2; k <= lengthRow(mat, i); k++)
	if(abs(mat->wt[GETJ(mat, mat->rows[i][k])]) > s)
	    s = abs(mat->wt[GETJ(mat, mat->rows[i][k])]);
    return s;
#endif
}

// locate columns of weight 3 and delete one of the rows, but just one!
int
findSuperfluousRowsFor2(int *tmp, int ntmp, sparse_mat_t *mat)
{
    dclist dcl;
    INT *Rj;
    int itmp = 0, j, k, kmin;

    for(dcl = mat->S[3]->next; dcl != NULL; dcl = dcl->next){
	j = dcl->j;
#ifndef USE_COMPACT_R
	Rj = mat->R[GETJ(mat, j)];
	// be semi-stupid
	for(kmin = 1; Rj[kmin] == -1; kmin++); // should stop at some point...!
	for(k = kmin+1; k <= mat->R[GETJ(mat, j)][0]; k++)
	    if(Rj[k] != -1)
		if(lengthRow(mat, Rj[k]) > lengthRow(mat, Rj[kmin]))
		    kmin = k;
	tmp[itmp++] = 0; // horrible trick not to change things elsewhere...!
	tmp[itmp++] = Rj[kmin];
	if(itmp >= ntmp)
	    // should not happen, but...
	    break;
#else
	fprintf(stderr, "R: NYI in findSuperfluousRowsFor2\n");
	exit(1);
#endif
    }
    return itmp;
}

// this guy is linear => total is quadratic, be careful!
int
findSuperfluousRows(int *tmp, int ntmp, sparse_mat_t *mat, int m)
{
    int i, itmp;

    if(m == 2){
	itmp = findSuperfluousRowsFor2(tmp, ntmp, mat);
	if(itmp != 0)
	    return itmp;
    }
    for(i = 0, itmp = 0; i < mat->nrows; i++)
        if(!isRowNull(mat, i)){
	    tmp[itmp++] = deleteScore(mat, i);
	    tmp[itmp++] = i;
	}
    // rows with largest score will be at the end
    qsort(tmp, itmp>>1, 2 * sizeof(int), cmp);
    return itmp;
}

// Delete superfluous rows s.t. nrows-ncols >= keep.
// Let's say we are at "step" m.
int
deleteSuperfluousRows(report_t *rep, sparse_mat_t *mat, int keep, int niremmax, int m)
{
    int nirem = 0, *tmp, ntmp, i;

    if((mat->rem_nrows - mat->rem_ncols) <= keep)
	return 0;
    ntmp = mat->rem_nrows << 1;
    tmp = (int *)malloc(ntmp * sizeof(int));
    ntmp = findSuperfluousRows(tmp, ntmp, mat, m);
    fprintf(stderr, "Number of removable rows[%d]: %d\n", m, ntmp>>1);
    // remove rows with largest score
    for(i = ntmp-1; i >= 0; i -= 2){
	if((nirem >= niremmax) || (mat->rem_nrows - mat->rem_ncols) <= keep)
	    break;
	if(mat->rows[tmp[i]] != NULL){
	    removeRowDefinitely(rep, mat, tmp[i]);
	    nirem++;
	}
    }
    free(tmp);
    return nirem;
}

void
merge(report_t *rep, sparse_mat_t *mat, int maxlevel, int verbose, int forbw)
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
	    fprintf(rep->outfile, "BWCOST: %lu\n", cost);
	if(forbw && (oldcost != 0) && (cost > oldcost)){
	    ncost++;
	    fprintf(stderr, "WARNING: New cost > old cost (%2.6e) [%d/%d]\n",
		    ((double)cost)-((double)oldcost), ncost, ncostmax);
	    if(ncost >= ncostmax){
		int nirem;

		fprintf(stderr, "WARNING: New cost > old cost %d times",
                        ncost);
		fprintf(stderr, " in a row:");
		nirem = deleteSuperfluousRows(rep, mat, mat->delta, 128, m);
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
	    njrem += removeSingletons(rep, mat);
	else
            njrem += mergeGe2(rep, mat, m, 0, verbose);
#if DEBUG >= 1
	checkData(mat);
#endif
#if 0
    {
	int i, ni = 0, j0 = 191342;

	for(i = 0; i < mat->nrows; i++)
	    if((mat->data[i].val != NULL) && hasCol(mat->rows, i, j0))
		ni++;
	printf("CHECK: %d %d\n", ni, mat->wt[GETJ(mat, j0)]);
    }
#endif
	fprintf(stderr, "=> nrows=%d ncols=%d (%d) njrem=%d\n",
		mat->rem_nrows, mat->rem_ncols, 
		mat->rem_nrows - mat->rem_ncols, njrem);
	inspectRowWeight(rep, mat);
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
	fprintf(rep->outfile, "BWCOSTMIN: %lu\n", bwcostmin);
	fprintf(stderr, "Minimal bwcost found: %lu\n", bwcostmin);
    }
}

static double
my_cost (unsigned long N, unsigned long c, int forbw)
{
    double dN = (double)N, dc = (double)c;
  
    if(forbw == 2){
	double K1 = .19e-9, K2 = 3.4e-05, K3 = 1.4e-10; // kinda average

	return ((K1+K3)*dN*dc+K2*dN*log(dN)*log(dN));
    }
    else if(forbw == 3)
      return dc / dN;
    else if(forbw <= 1)
	return dN * dc;
    return 0;
}

void
mergeForColumn2(report_t *rep, sparse_mat_t *mat, int *njrem, 
		double *totopt, double *totfill, double *totMST, 
		double *totdel, int useMST, INT j)
{
    double tt, tfill, tMST;

    mergeForColumn(rep, &tt, &tfill, &tMST, mat, mat->wt[j], j, useMST);
    *totopt += tt;
    *totfill += tfill;
    *totMST += tMST;
    tt = seconds();
    *njrem += deleteHeavyColumns(rep, mat);
    *totdel += (seconds()-tt);
}

// njrem is the number of columns removed.
void
doOneMerge(report_t *rep, sparse_mat_t *mat, int *njrem, double *totopt, double *totfill, double *totMST, double *totdel, int m, int maxdo, int useMST, int verbose)
{
    dclist dcl;
    INT j;

    if(m == 1){
	// do all of these
	*njrem += removeSingletons(rep, mat);
	fprintf(stderr, "T=%d after m=1\n", (int)(seconds()));
    }
    else if(m == 2){
#if 0
	fprintf(stderr, "Performing %d merges for m=%d at %2.2lf\n",
		maxdo, m, seconds());
#endif
	*njrem += mergeGe2(rep, mat, m, maxdo, verbose);
    }
    else{
	//	    fprintf(stderr, "Performing one merge for m=%d\n", m);
	dcl = mat->S[m]->next;
	j = dcl->j;
	mergeForColumn2(rep, mat, njrem, totopt, totfill, totMST, totdel,
			useMST, j);
    }
}

int
number_of_superfluous_rows(sparse_mat_t *mat)
{
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
    return ni2rem;
}

void
mergeOneByOne(report_t *rep, sparse_mat_t *mat, int maxlevel, int verbose, int forbw, double ratio, int coverNmax)
{
    double totopt = 0.0, totfill = 0.0, totMST = 0.0, totdel = 0.0;
    double bwcostmin = 0.0, oldbwcost = 0.0, bwcost = 0.0;
    int old_nrows, old_ncols, m = 2, njrem = 0, ncost = 0, ncostmax, njproc;
    int mmax = 0, target = 10000, ni2rem;
#ifdef USE_MARKOWITZ
    INT j, mkz;
#endif

    fprintf(stderr, "Using mergeOneByOne\n");
    ncostmax = 20; // was 5
    njproc = 0;
    while(1){
	oldbwcost = bwcost;
	old_nrows = mat->rem_nrows;
	old_ncols = mat->rem_ncols;
#ifndef USE_MARKOWITZ
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
	doOneMerge(rep, mat, &njrem, &totopt, &totfill, &totMST, &totdel,
		   m, 0, 1, verbose);
#else
	MkzPopQueue(&j, &mkz, mat->MKZQ, mat->MKZA);
	fprintf(stderr, "I pop j=%d mkz=%d\n", j, mkz); 
	mergeForColumn2(rep, mat, &njrem, &totopt, &totfill, &totMST, &totdel,
			1, j);
#endif
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
	    ni2rem = number_of_superfluous_rows(mat);
	    deleteSuperfluousRows(rep, mat, mat->delta, ni2rem, m);
	    inspectRowWeight(rep, mat);
	    fprintf(stderr, "T=%d mmax=%d N=%d nc=%d (%d)",
		    (int)seconds(),
		    mmax,
		    mat->rem_nrows, mat->rem_ncols, 
		    mat->rem_nrows - mat->rem_ncols);
	    if(forbw == 2)
		fprintf(stderr, " bw=%1.2e", bwcost);
	    else if(forbw == 3)
		fprintf(stderr, " cN=%1.2e",
			((double) mat->rem_nrows) * ((double) mat->weight));
	    else if(forbw <= 1)
		fprintf(stderr, " cN=%1.2e", bwcost);
	    fprintf(stderr, " c/N=%2.2lf\n", 
		    ((double)mat->weight)/((double)mat->rem_ncols));
	    // njrem=%d at %2.2lf\n",
	    if((forbw != 0) && (forbw != 3))
		// what a trick!!!!
		fprintf(rep->outfile, "BWCOST: %1.0f\n", bwcost);
	    target = njproc + 10000;
	}
	if((bwcostmin == 0) || (bwcost < bwcostmin)){
	    bwcostmin = bwcost;
	    if((forbw != 0) && (forbw != 3))
		// what a trick!!!!
		fprintf(rep->outfile, "BWCOST: %1.0f\n", bwcost);
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
		fprintf(stderr, "c/N too high, stopping [%1.2e]\n", bwcost);
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
		nirem = deleteSuperfluousRows(rep, mat, mat->delta, 128, m);
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
    deleteSuperfluousRows(rep, mat, mat->delta, 
			  mat->rem_nrows-mat->rem_ncols+mat->delta, -1);
    if((forbw != 0) && (forbw != 3)){
	fprintf(rep->outfile, "BWCOSTMIN: %1.0f\n", bwcostmin);
	fprintf(stderr, "Minimal bwcost found: %1.0f\n", bwcostmin);
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
	if(mat->wt[GETJ(mat, j)] != 0){
	    fprintf(ofile, "J %d %d\n", j, abs(mat->wt[GETJ(mat, j)]));
	    if(abs(mat->wt[GETJ(mat, j)]) <= mat->mergelevelmax)
		fprintf(stderr, "Gasp: %d %d\n", j, mat->wt[GETJ(mat, j)]);
	}
}

//////////////////////////////////////////////////////////////////////
// 
// Resume section: very much inspired by replay.c, of course...!
//
//////////////////////////////////////////////////////////////////////

// A line is "i i1 ... ik".
int
indicesFromString(INT *ind, char *str)
{
    int ni = 0;
    char *tok = strtok(str, " ");

    while(tok != NULL){
	ind[ni++] = atoi(tok);
	tok = strtok(NULL, " ");
    }
    return ni;
}

// A line is "i i1 ... ik".
// If i >= 0 then
//     row[i] is to be added to rows i1...ik and destroyed at the end of
//     the process.
//     Works also is i is alone (hence: destroyed row).
// If i < 0 then
//     row[-i-1] is to be added to rows i1...ik and NOT destroyed.
//
void
doAllAdds(report_t *rep, sparse_mat_t *mat, char *str)
{
    INT ind[MERGE_LEVEL_MAX], i0;
    int ni, sg, k;

    ni = indicesFromString(ind, str);
    if(ind[0] < 0){
	i0 = -ind[0]-1;
	sg = -1;
    }
    else{
	i0 = ind[0];
	sg = 1;
    }
    for(k = 1; k < ni; k++)
	addRowsAndUpdate(mat, ind[k], i0, -1);
    reportn(rep, ind, ni);
    if(sg > 0){
	// when ni == 1, then the corresponding row was too heavy and dropped
	// unless m = 1 was used, in which case S[1] will contain j and
	// can dispensed of next time (?) In the case of S[0], it'll be done
	// later on also (?)
	removeRowAndUpdate(mat, i0);
	destroyRow(mat, i0);
	mat->rem_nrows--;
	// the number of active j's is recomputed, anyway, later on
    }
}

// resumename is a file of the type mergehis.
void
resume(report_t *rep, sparse_mat_t *mat, char *resumename)
{
    FILE *resumefile = fopen(resumename, "r");
    char str[STR_LEN_MAX];
    unsigned long addread = 0;
    int nactivej;
    INT j;

    fprintf(stderr, "Resuming computations from %s\n", resumename);
    // skip first line containing nrows ncols
    fgets(str, STR_LEN_MAX, resumefile);
    fprintf(stderr, "Reading row additions\n");
    while(fgets(str, STR_LEN_MAX, resumefile)){
	addread++;
	if((addread % 100000) == 0)
	    fprintf(stderr, "%lu lines read at %2.2lf\n", addread, seconds());
	if(str[strlen(str)-1] != '\n'){
	    fprintf(stderr, "Gasp: not a complete a line!");
	    fprintf(stderr, " I stop reading and go to the next phase\n");
	    break;
	}
	if(strncmp(str, "BWCOST", 6) != 0)
	    doAllAdds(rep, mat, str);
    }
    fclose(resumefile);
    for(j = mat->jmin; j < mat->jmax; j++)
	if((mat->wt[GETJ(mat, j)] == 0) && (mat->A[GETJ(mat, j)] != NULL))
	    // be sure j was removed...
	    removeColumnAndUpdate(mat, j);
    nactivej = 0;
    for(j = mat->jmin; j < mat->jmax; j++)
	if(mat->wt[GETJ(mat, j)] != 0)
	    nactivej++;
    mat->rem_ncols = nactivej;
    fprintf(stderr, "At the end of resume, we have");
    fprintf(stderr, " nrows=%d ncols=%d (%d)\n",
	    mat->rem_nrows, mat->rem_ncols, mat->rem_nrows-mat->rem_ncols);
}
