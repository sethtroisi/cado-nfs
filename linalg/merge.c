/* 
 * Program: merge
 * Author : F. Morain
 * Purpose: merging relations
 * 
 * Algorithm: digest and interpolation from Cavallar.
 *
 */

// TODO: use compact lists...!
// TODO: keep track of the number of columns dropped at first sight
// TODO: merge SWAR into sparse?

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>
#include <string.h>

#include "utils/utils.h"

#define DEBUG 0
#define TEX 0

#define USE_USED 0
#if USE_USED >= 1
char *used;
#endif

#define MERGE_LEVEL_MAX 30

// 0 means dummy filtering using slow methods
// 1 means:
//  * fast with optimized-though-memory-consuming data structure;
//  * heavy columns are killed.
// 2 means:
//  * fast as above;
//  * heavy columns are present, but not in S[]; this yields more acurate
//    weights.
#define USE_MERGE_FAST 2

typedef struct {
  int len;    /* number of non-zero entries */
  int *val;   /* array of size len containing column indices of non-zero
                 entries */
} rel_t;

typedef struct {
  int nrows;
  int ncols;
  int rem_nrows;     /* number of remaining rows */
  int rem_ncols;     /* number of remaining columns */
  rel_t *data;
  int *wt; /* weight of prime j, <= 1 for a deleted prime */
} sparse_mat_t;

static int 
cmp(const void *p, const void *q) {
    int x = *((int *)p);
    int y = *((int *)q);
    return (x <= y ? -1 : 1);
}

// doubly chained lists
typedef struct dclist{
    int j;
    struct dclist *prev, *next;
} *dclist;

typedef struct{
    dclist *S, *A;
    int *Wj;
    int **R;
    int cwmax, mergelevelmax;
} swar_t;

dclist
dclistCreate(int j)
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
	fprintf(file, " -> %d", dcl->j);
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
	fprintf(file, "%d\\\\\n", dcl->j);
	dcl = dcl->next;
    }
    fprintf(stderr, "\\end{array}\n");
}

dclist
dclistInsert(dclist dcl, int j)
{
    dclist newdcl = dclistCreate(j);

    newdcl->next = dcl->next;
    newdcl->prev = dcl;
    if(newdcl->next != NULL)
	newdcl->next->prev = newdcl;
    dcl->next = newdcl;
    return newdcl;
}

// TODO: ncols could be a new renumbered ncols...
void
initSWAR(swar_t *SWAR, int ncols, int cwmax, int mergelevelmax)
{
#if USE_MERGE_FAST
    // cwmax+2 to prevent bangs
    dclist *S = (dclist *)malloc((cwmax+2) * sizeof(dclist));
    dclist *A = (dclist *)malloc(ncols * sizeof(dclist));
    int *Wj = (int *)malloc(ncols * sizeof(int));
    int **R;
    int j;
#endif
    SWAR->cwmax = cwmax;
    SWAR->mergelevelmax = mergelevelmax;
#if USE_MERGE_FAST
    // S[0] has a meaning at least for temporary reasons
    for(j = 0; j <= cwmax+1; j++)
	S[j] = dclistCreate(-1);
    R = (int **)malloc(ncols * sizeof(int *));
    SWAR->S = S;
    SWAR->A = A;
    SWAR->R = R;
    SWAR->Wj = Wj;
#endif
}

void
fillSWAR(swar_t *SWAR, int *usecol, int ncols)
{
    int j, *Rj;

    for(j = 0; j < ncols; j++){
	if(usecol[j] <= SWAR->cwmax){
	    SWAR->Wj[j] = usecol[j]; // TODO: simplify this!
	    SWAR->A[j] = dclistInsert(SWAR->S[usecol[j]], j);
#  if DEBUG >= 1
	    fprintf(stderr, "Inserting %d in S[%d]:", j, usecol[j]);
	    dclistPrint(stderr, SWAR->S[usecol[j]]->next);
	    fprintf(stderr, "\n");
#  endif
	    Rj = (int *)malloc((usecol[j]+1) * sizeof(int));
	    Rj[0] = 0; // last index used
	    SWAR->R[j] = Rj;
	}
	else{
#if USE_MERGE_FAST <= 1
	    SWAR->Wj[j] = -1;
#else
	    SWAR->Wj[j] = -usecol[j]; // trick!!!
#endif
	    SWAR->A[j] = NULL; // TODO: renumber j's?????
	    SWAR->R[j] = NULL;
	}
    }
}

// TODO
void
closeSWAR(swar_t *SWAR)
{
}

void
printStats(swar_t *SWAR)
{
    int w;

    for(w = 0; w <= SWAR->cwmax; w++)
	fprintf(stderr, "I found %d primes of weight %d\n",
		dclistLength(SWAR->S[w]->next), w);
}

void
checkData(sparse_mat_t *mat, swar_t *SWAR)
{
    int j;

    fprintf(stderr, "Checking final data\n");
    for(j = 0; j < mat->ncols; j++)
	if(SWAR->Wj[j] > 0)
	    fprintf(stderr, "Wj[%d] = %d\n", j, SWAR->Wj[j]);
}

// dump for debugging reasons
void
printSWAR(swar_t *SWAR, int ncols)
{
    int j, k, w;

    fprintf(stderr, "===== S is\n");
    for(w = 0; w <= SWAR->cwmax; w++){
	fprintf(stderr, "  S[%d] ", w);
	dclistPrint(stderr, SWAR->S[w]->next);
	fprintf(stderr, "\n");
    }
    fprintf(stderr, "===== Wj is\n");
    for(j = 0; j < ncols; j++)
	fprintf(stderr, "  Wj[%d]=%d\n", j, SWAR->Wj[j]);
    fprintf(stderr, "===== R is\n");
    for(j = 0; j < ncols; j++){
	if(SWAR->R[j] == NULL)
	    continue;
	fprintf(stderr, "  R[%d]=", j);
	for(k = 1; k <= SWAR->R[j][0]; k++)
	    fprintf(stderr, " %d", SWAR->R[j][k]);
	fprintf(stderr, "\n");
    }
}

void
matrix2tex(sparse_mat_t *mat, int *W)
{
    int i, j, k, *tab;

    tab = (int *)malloc(mat->ncols * sizeof(int));
    fprintf(stderr, "$$\\begin{array}{l|");
    for(j = 0; j < mat->ncols; j++)
	fprintf(stderr, "c");
    fprintf(stderr, "}\n");
    for(i = 0; i < mat->nrows; i++){
	memset(tab, 0, mat->ncols * sizeof(int));
	for(k = 0; k < mat->data[i].len; k++)
	    tab[mat->data[i].val[k]] = 1;
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
texSWAR(sparse_mat_t *mat, swar_t *SWAR)
{
    int j, k, w;

    fprintf(stderr, "\\end{verbatim}\n");
    matrix2tex(mat, SWAR->Wj);
    fprintf(stderr, "$$\\begin{array}{|");
    for(w = 0; w <= SWAR->cwmax; w++)
	fprintf(stderr, "r|");
    fprintf(stderr, "|r|}\\hline\n");
    for(w = 0; w <= SWAR->cwmax+1; w++){
	fprintf(stderr, "S[%d]", w);
	if(w < SWAR->cwmax+1)
	    fprintf(stderr, "&");
    }
    fprintf(stderr, "\\\\\\hline\n");
    for(w = 0; w <= SWAR->cwmax+1; w++){
	dclistTex(stderr, SWAR->S[w]->next);
	if(w < SWAR->cwmax+1)
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
	fprintf(stderr, "& %d", SWAR->Wj[j]);
    fprintf(stderr, "\\\\\\hline\n\\end{array}$$\n");
    fprintf(stderr, "\\begin{verbatim}\n");
    for(j = 0; j < mat->ncols; j++){
	if(SWAR->R[j] == NULL)
	    continue;
	fprintf(stderr, "  R[%d]=", j);
	for(k = 1; k <= SWAR->R[j][0]; k++)
	    fprintf(stderr, " %d", SWAR->R[j][k]);
	fprintf(stderr, "\n");
    }
}

// terrific hack: everybody on the same line
// the first is to be destroyed in replay!!!
void
reportn(int *ind, int n)
{
    int i;

#if DEBUG >= 1
    fprintf(stderr, "Reporting for n=%d\n", n);
#endif
    for(i = 0; i < n; i++){
	printf("%d", ind[i]);
	if(i < n-1)
	    printf(" ");
    }
    printf("\n");
}

void
report1(int i)
{
    reportn(&i, 1);
}

void
report2(int i1, int i2)
{
    int tmp[2];
    tmp[0] = i1;
    tmp[1] = i2;
    reportn(tmp, 2);
}

void
inspectWeight(int *usecol, FILE *purgedfile, int nrows, int ncols, int cwmax)
{
    int i, j, ret, x, nc;

    memset(usecol, 0, ncols * sizeof(int));
    for(i = 0; i < nrows; i++){
	ret = fscanf(purgedfile, "%d", &j); // unused index to rels file
	assert (ret == 1);
	ret = fscanf (purgedfile, "%d", &nc);
	assert (ret == 1);
	for(j = 0; j < nc; j++){
	    ret = fscanf(purgedfile, "%d", &x);
	    assert (ret == 1);
	    usecol[x]++;
	}
    }
}

/* Reads a matrix file, and puts in mat->wt[j], 0 <= j < ncols, the
   weight of column j (adapted from matsort.c).

   We skip columns that are too heavy, that is columns for which usecol[j]==0.

*/
void
readmat(swar_t *SWAR, FILE *file, sparse_mat_t *mat, int *usecol)
{
    int ret;
    int i, j;
    int nc, x, buf[100], ibuf;

    ret = fscanf (file, "%d %d", &(mat->nrows), &(mat->ncols));
    assert (ret == 2);
    
    fprintf(stderr, "Reading matrix of %d rows and %d columns: excess is %d\n",
	    mat->nrows, mat->ncols, mat->nrows - mat->ncols);
    mat->rem_nrows = mat->nrows;
    
    mat->wt = (int*) malloc (mat->ncols * sizeof (int));
    for (j = 0; j < mat->ncols; j++)
	mat->wt[j] = 0;
    
    mat->data = (rel_t*) malloc (mat->nrows * sizeof (rel_t));
    
    for (i = 0; i < mat->nrows; i++){
	ret = fscanf(file, "%d", &j); // unused index to rels file
	assert (ret == 1);
	ret = fscanf (file, "%d", &nc);
	assert (ret == 1);
	if(nc == 0){
	    mat->data[i].len = 0;
	    mat->data[i].val = NULL;
	    mat->rem_nrows--;
	}
	else{
	    for(j = 0, ibuf = 0; j < nc; j++){
		ret = fscanf(file, "%d", &x);
#if DEBUG >= 1
		fprintf(stderr, "i = %d, j = %d, x = %d\n", i, j, x);
#endif
		assert (ret == 1);
		assert (0 <= x && x < mat->ncols);
#if USE_MERGE_FAST <= 1
		if(usecol[x] <= SWAR->cwmax){
		    buf[ibuf++] = x;
		    mat->wt[x]++;
		}
#else
		// always store x
		buf[ibuf++] = x;
		mat->wt[x]++;
#endif
		if(usecol[x] <= SWAR->cwmax){
		    SWAR->R[x][0]++;
		    SWAR->R[x][SWAR->R[x][0]] = i;
		}
	    }
	    mat->data[i].len = ibuf;
	    mat->data[i].val = (int*) malloc (ibuf * sizeof (int));
	    memcpy(mat->data[i].val, buf, ibuf * sizeof(int));
	    // sort indices in val to ease row merges
	    qsort(mat->data[i].val, ibuf, sizeof(int), cmp);
	}
    }
    // we need to keep informed of what really happens; this will be an upper
    // bound on the number of active columns, I guess
    mat->rem_ncols = mat->ncols;
    printStats(SWAR);
}

void
print_row(sparse_mat_t *mat, int i)
{
    int k;
    
    for(k = 0; k < mat->data[i].len; k++)
	fprintf(stderr, " %d", mat->data[i].val[k]);
}

void
destroyRow(sparse_mat_t *mat, int i)
{
    free(mat->data[i].val);
    mat->data[i].val = NULL;
    mat->data[i].len = 0;
}

void
removeWeight(sparse_mat_t *mat, int i)
{
    int k;

    for(k = 0; k < mat->data[i].len; k++)
	mat->wt[mat->data[i].val[k]]--;
}

int
hasCol(sparse_mat_t *mat, int i, int j)
{
    int k;

    for(k = 0; k < mat->data[i].len; k++)
	if(mat->data[i].val[k] == j)
	    return 1;
    return 0;
}

// i1 += i2, mat->wt is updated at the same time.
void
addRows(sparse_mat_t *mat, int i1, int i2)
{
    int k1, k2, k, len, *tmp, *tmp2;

#if DEBUG >= 2
    fprintf(stderr, "row[%d] =", i1); print_row(mat, i1);
    fprintf(stderr, "\n");
    fprintf(stderr, "row[%d] =", i2); print_row(mat, i2);
    fprintf(stderr, "\n");
#endif
    // merge row[i1] and row[i2] in i1...
    len = mat->data[i1].len + mat->data[i2].len;
    tmp = (int *)malloc(len * sizeof(tmp));
    k = k1 = k2 = 0;

    // TODO: how do we update ad????

    // i1 is to disappear, replaced by a new one
    removeWeight(mat, i1);
    // loop while everybody is here
    while((k1 < mat->data[i1].len) && (k2 < mat->data[i2].len)){
	if(mat->data[i1].val[k1] < mat->data[i2].val[k2]){
	    tmp[k++] = mat->data[i1].val[k1++];
	    mat->wt[mat->data[i1].val[k1-1]]++;
	}
	else if(mat->data[i1].val[k1] > mat->data[i2].val[k2]){
            tmp[k++] = mat->data[i2].val[k2++];
	    mat->wt[mat->data[i2].val[k2-1]]++;
	}
	else{
#if DEBUG >= 1
	    fprintf(stderr, "ADD[%d=(%ld,%lu), %d=(%ld,%lu)]: new w[%d]=%lu\n", 
		    i1, mat->data[i1].a, mat->data[i1].b,
		    i2, mat->data[i2].a, mat->data[i2].b,
		    mat->data[i1].val[k1], 
		    mat->wt[mat->data[i1].val[k1]]);
#endif
	    k1++;
	    k2++;
	}
    }
    // finish with k1
    for( ; k1 < mat->data[i1].len; k1++){
	tmp[k++] = mat->data[i1].val[k1];
	mat->wt[mat->data[i1].val[k1]]++;
    }
    // finish with k2
    for( ; k2 < mat->data[i2].len; k2++){
	tmp[k++] = mat->data[i2].val[k2];
	mat->wt[mat->data[i2].val[k2]]++;
    }
    // destroy and copy back
    free(mat->data[i1].val);
    tmp2 = (int *)malloc(k * sizeof(int));
    memcpy(tmp2, tmp, k * sizeof(int));
    mat->data[i1].len = k;
    mat->data[i1].val = tmp2;
    free(tmp);
#if DEBUG >= 2
    fprintf(stderr, "row[%d]+row[%d] =", i1, i2);
    print_row(mat, i1); fprintf(stderr, "\n");
#endif
}

void
merge2rows(sparse_mat_t *mat, int j)
{
    int i1, i2, ok = 0, k1, k2;

    for(i1 = 0; i1 < mat->nrows; i1++){
	if(mat->data[i1].val == NULL)
	    continue;
	for(k1 = 0; k1 < mat->data[i1].len; k1++)
	    if(mat->data[i1].val[k1] == j){
		ok = 1;
		break;
	    }
	if(ok)
	    break;
    }
    if(!ok){
	printf("Pb in merge2rows: no row containing %d\n", j);
	return;
    }
    ok = 0;
    for(i2 = i1+1; i2 < mat->nrows; i2++){
	if(mat->data[i2].val == NULL)
	    continue;
	for(k2 = 0; k2 < mat->data[i2].len; k2++)
            if(mat->data[i2].val[k2] == j){
                ok = 1;
                break;
            }
	if(ok)
	    break;
    }
    if(!ok){
	printf("Pb in merge2rows: no 2nd row containing %d\n", j);
	return;
    }
#if USE_USED >= 1
    if(used[i1])
	fprintf(stderr, "R%d used %d times\n", i1, (int)used[i1]);
    if(used[i2])
	fprintf(stderr, "R%d used %d times\n", i2, (int)used[i2]);
    used[i1]++;
    used[i2] = -1;
#endif
#if DEBUG >= 1
    fprintf(stderr, "Merging rows: %d += %d [j=%d]\n", i1, i2, j);
#endif
    printf("%d %d\n", i2, i1); // note the order!!!
    addRows(mat, i1, i2);
    removeWeight(mat, i2);
    destroyRow(mat, i2);
    mat->rem_nrows--;
    mat->rem_ncols--;
#if DEBUG >= 1
    matrix2tex(mat);
#endif
}

void
merge2(sparse_mat_t *mat, int nb_merge_max)
{
    int j, nb_merge = 0;

#if DEBUG >= 1
    matrix2tex(mat);
#endif
    for(j = 0; j < mat->ncols; j++){
	if(nb_merge == nb_merge_max){
	    fprintf(stderr, "Warning: nb_merge has reached its maximum:");
	    fprintf(stderr, " %d\n", nb_merge_max);
	    break; // useful for debugging and more...!
	}
	if(mat->wt[j] == 2){
#if DEGUB >= 1
	    fprintf(stderr, "j = %d has weight 2\n", j);
#endif
	    nb_merge++;
	    merge2rows(mat, j);
	    if(!(nb_merge % 1000))
		fprintf(stderr, "nrows=%d ncols=%d\n",
			mat->rem_nrows, mat->rem_ncols);
	}
    }
}

// what is the weight of the sum of Ra and Rb?
int
weightSum(sparse_mat_t *mat, int i1, int i2)
{
    int k1 = 0, k2 = 0, w = 0;

    while((k1 < mat->data[i1].len) && (k2 < mat->data[i2].len)){
	if(mat->data[i1].val[k1] < mat->data[i2].val[k2]){
	    k1++;
	    w++;
	}
	else if(mat->data[i1].val[k1] > mat->data[i2].val[k2]){
	    k2++;
	    w++;
	}
	else{
	    k1++; k2++;
	}
    }
    w += (k1 > mat->data[i1].len ? 0 : mat->data[i1].len-k1+1);
    w += (k2 > mat->data[i2].len ? 0 : mat->data[i2].len-k2+1);
    return w;
}

int
findAllRowsWithGivenj(int *ind, sparse_mat_t *mat, int j, int nb)
{
    int i, k, r = 0;

    // TODO: special hack for nb==2???
    for(i = 0; i < mat->nrows; i++){
	if(mat->data[i].val == NULL)
	    continue;
	for(k = 0; k < mat->data[i].len-1; k++) // trick!
	    if(mat->data[i].val[k] >= j)
		break;
	if(mat->data[i].val[k] == j){
	    ind[r++] = i;
	    if(r == nb)
		return 1;
	}
    }
    return 0;
}

int
findBestIndex3(sparse_mat_t *mat, int *ind)
{
    int w01, w02, w12, wtot, w0, w1, w2, k, i;

    w01 = weightSum(mat, ind[0], ind[1]);
    w02 = weightSum(mat, ind[0], ind[2]);
    w12 = weightSum(mat, ind[1], ind[2]);
#if DEBUG >= 1
    fprintf(stderr,"W(0+1)=%d\n", w01);
    fprintf(stderr,"W(0+2)=%d\n", w02);
    fprintf(stderr,"W(1+2)=%d\n", w12);
#endif
    // first scenario r1 += r0, r2 += r0, r0 disappears
    wtot = 0;
    for(k = 0; k < 3; k++)
	wtot += mat->data[ind[k]].len;
    w0 = w01+w02-wtot;
    w1 = w01+w12-wtot;
    w2 = w02+w12-wtot;
#if DEBUG >= 1
    fprintf(stderr, "Using r0 = r[%d]: %d\n", ind[0], w0);
    fprintf(stderr, "Using r1 = r[%d]: %d\n", ind[1], w1);
    fprintf(stderr, "Using r2 = r[%d]: %d\n", ind[2], w2);
#endif
    if(w0 < w1)
	i = (w0 < w2 ? 0 : 2);
    else
	// w0 >= w1
	i = (w1 < w2 ? 1 : 2);
    return i;
}

//////////////////////////////////////////////////////////////////////
// Prim

void
popQueue(int *s, int *t, int **Q, int *fQ, int *lQ)
{
    *s = Q[*fQ][0];
    *t = Q[*fQ][1];
    *fQ += 1;
}

void
printQueue(int **Q, int fQ, int lQ)
{
    int i;

    for(i = fQ; i < lQ; i++)
	fprintf(stderr, "Q[%d] = [%d, %d, %d]\n", i,Q[i][0], Q[i][1], Q[i][2]);
}

void
addEdge(int **Q, int *fQ, int *lQ, int u, int v, int Auv)
{
    int i, j;

    for(i = *fQ; i < *lQ; i++)
	if(Auv < Q[i][2])
	    break;
    // shift everybody: Auv >= Q[i][2]
    for(j = *lQ; j > i; j--){
	Q[j][0] = Q[j-1][0];
	Q[j][1] = Q[j-1][1];
	Q[j][2] = Q[j-1][2];
    }
    Q[i][0] = u;
    Q[i][1] = v;
    Q[i][2] = Auv;
    *lQ += 1;
}

// Add all neighbors of u, which is already in T; all edges with v in T
// are not useful anymore.
void
addAllEdges(int **Q, int *fQ, int *lQ, int u, int *father, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
    int v;

    for(v = 0; v < m; v++)
	if((v != u) && (father[v] < 0))
	    addEdge(Q, fQ, lQ, u, v, A[u][v]);
}

// height[0] = 0, etc. Returns hmax.
int
minimalSpanningTreeWithPrim(int *father, int *height, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
    int **Q, m2 = m*m, u, fQ = 0, lQ = 0, s, t, i, nV, hmax = 0;

    // over-conservative
    Q = (int **)malloc(m2 * sizeof(int *));
    for(i = 0; i < m2; i++)
	Q[i] = (int *)malloc(3 * sizeof(int));
    // nodes of T
    for(i = 0; i < m; i++){
	father[i] = -1;
	height[i] = 0;
    }
    u = 0;
    father[u] = u; // for the root, isn't it?
    nV = m-1;
    addAllEdges(Q, &fQ, &lQ, u, father, A, m);
#if DEBUG >= 1
    printQueue(Q, fQ, lQ);
#endif
    // active part of Q is Q[fQ..lQ[
    assert(fQ == 0);
    assert(lQ == (m-1));
    while(fQ != lQ){
	// while queue is non empty
	// pop queue
	popQueue(&s, &t, Q, &fQ, &lQ);
#if DEBUG >= 1
	fprintf(stderr, "Popping a = (%d, %d)\n", s, t);
#endif
	if(father[t] == -1){
	    // t does not close a cycle
	    // T[u] <- T[u] union (s, t)
#if DEBUG >= 1
	    fprintf(stderr, "new edge: (%d, %d)\n", s, t);
#endif
	    father[t] = s;
	    height[t] = height[s]+1;
	    if(height[t] > hmax)
		hmax = height[t];
	    nV--;
	    if(nV == 0)
		break;
	    addAllEdges(Q, &fQ, &lQ, t, father, A, m);
#if DEBUG >= 1
	    printQueue(Q, fQ, lQ);
#endif
	}
    }
    for(i = 0; i < m2; i++)
	free(Q[i]);
    free(Q);
    return hmax;
}

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

int
minimalSpanningTree(int *father, int *height, sparse_mat_t *mat, int m, swar_t *SWAR, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
    return minimalSpanningTreeWithPrim(father, height, A, m);
    //    minimalSpanningTreeWithKruskal(A, m);
}

void
fillRowAddMatrix(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], sparse_mat_t *mat, int m, int *ind)
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
findBestIndex(sparse_mat_t *mat, int m, int *ind)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], i, j, imin, wmin, w;

    assert(m <= MERGE_LEVEL_MAX);
    if(m == 2)
	return 0;
#if 0 // obsolete...
    if(m == 3)
	return findBestIndex3(mat, ind);
#endif
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
merge_m_slow(sparse_mat_t *mat, int m)
{
    int totfindrows = 0, totfindbest = 0, totadd = 0, tt;
    int *ind, j, k, i, njproc = 0;

#if DEBUG >= 1
    fprintf(stderr, "Weight %d:", m);
#endif
    ind = (int *)malloc(m * sizeof(int));
    for(j = 0; j < mat->ncols; j++){
	if(mat->wt[j] != m)
	    continue;
	njproc++;
	if(!(njproc % 1000))
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
		addRows(mat, ind[k], ind[i]);
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
	reportn(ind, m);
	totadd += (cputime()-tt);
	removeWeight(mat, ind[0]);
	destroyRow(mat, ind[0]);
	mat->rem_nrows--;
	mat->rem_ncols--;
	mat->wt[j] = 0;
    }
    fprintf(stderr, "TIME: findrows=%d findbest=%d add=%d\n",
	    totfindrows, totfindbest, totadd);
    free(ind);
}

// w(j) has decreased in such a way that it can be incorporated in the
// SWAR structure. We need find all the info we need to do that. Note this
// can be costly. We must not store row index i0, since j was discovered
// when row i was treated.
void
incorporateColumn(swar_t *SWAR, sparse_mat_t *mat, int j, int i0)
{
    int i, ni, *Rj, wj;

    wj = -SWAR->Wj[j]-1;
    Rj = (int *)malloc((wj+1) * sizeof(int));
    // find all rows in which j appears
    for(i = 0, ni = 1; i < mat->nrows; i++)
	if(mat->data[i].val != NULL)
	    if((i != i0) && hasCol(mat, i, j))
#if 1
		Rj[ni++] = i;
#else
                ni++;
    fprintf(stderr, "iC: %d %d\n", j, ni);
#endif
    Rj[0] = ni-1;
    SWAR->R[j] = Rj;
    SWAR->Wj[j] = wj;
    SWAR->A[j] = dclistInsert(SWAR->S[wj], j);
}

int
getNextj(dclist dcl)
{
    int j;
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
remove_j_from_S(swar_t *SWAR, int j)
{
    dclist dcl = SWAR->A[j], foo;

#if DEBUG >= 2
    int ind = SWAR->Wj[j];
    if(ind > SWAR->cwmax)
        ind = SWAR->cwmax+1;
    fprintf(stderr, "S[%d]_b=", ind);
    dclistPrint(stderr, SWAR->S[ind]->next); fprintf(stderr, "\n");
    fprintf(stderr, "dcl="); dclistPrint(stderr, dcl); fprintf(stderr, "\n");
#endif
    foo = dcl->prev;
    foo->next = dcl->next;
    if(dcl->next != NULL)
	dcl->next->prev = foo;
#if DEBUG >= 2
    fprintf(stderr, "S[%d]_a=", ind);
    dclistPrint(stderr, SWAR->S[ind]->next); fprintf(stderr, "\n");
#endif
    free(dcl);
}

void
destroyRj(swar_t *SWAR, int j)
{
    free(SWAR->R[j]);
    SWAR->R[j] = NULL;
}

void
remove_j_from_SWAR(swar_t *SWAR, int j)
{
    remove_j_from_S(SWAR, j);
    SWAR->A[j] = NULL;
    SWAR->Wj[j] = -1;
    destroyRj(SWAR, j);
}

// Don't touch to R[j][0]!!!!
void
remove_i_from_Rj(swar_t *SWAR, int i, int j)
{
    // be dumb for a while
    int k;
    
    for(k = 1; k <= SWAR->R[j][0]; k++)
	if(SWAR->R[j][k] == i){
#if DEBUG >= 1
	    fprintf(stderr, "Removing row %d from R[%d]\n", i, j);
#endif
	    SWAR->R[j][k] = -1;
	    break;
	}
}

void
add_i_to_Rj(swar_t *SWAR, int i, int j)
{
    // be semi-dumb for a while
    int k;
    
#if DEBUG >= 1
    fprintf(stderr, "Adding row %d to R[%d]\n", i, j);
#endif
    for(k = 1; k <= SWAR->R[j][0]; k++)
	if(SWAR->R[j][k] == -1)
	    break;
    if(k <= SWAR->R[j][0]){
	// we have found a place where it is -1
	SWAR->R[j][k] = i;
    }
    else{
#if DEBUG >= 1
	fprintf(stderr, "WARNING: reallocing things in add_i_to_Rj for R[%d]\n", j);
#endif
	int l = SWAR->R[j][0]+2;
	SWAR->R[j] = (int *)realloc(SWAR->R[j], l * sizeof(int));
	SWAR->R[j][l-1] = i;
	SWAR->R[j][0] = l-1;
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

// A[j] contains the address where j is stored
void
removeCellFast(sparse_mat_t *mat, int i, int j, swar_t *SWAR)
{
    int ind;

    // update weight
#if DEBUG >= 1
    fprintf(stderr, "Moving j=%d from S[%d] to S[%d]\n",
	    j, SWAR->Wj[j], SWAR->Wj[j]-1);
#endif
    ind = SWAR->Wj[j] = decrS(SWAR->Wj[j]);
#if USE_MERGE_FAST > 1
    if(ind < 0){
	// what if abs(weight) becomes below cwmax?
	// we should incorporate the column and update the data structure, no?
	if(abs(ind) <= SWAR->mergelevelmax){
	    fprintf(stderr, "WARNING: column %d becomes light at %d...!\n",
		    j, abs(ind));
	    incorporateColumn(SWAR, mat, j, i);
	}
	return;
    }
#endif
    remove_j_from_S(SWAR, j);
    if(SWAR->Wj[j] > SWAR->cwmax)
	ind = SWAR->cwmax+1;
    // update A[j]
#if DEBUG >= 2
    fprintf(stderr, "S[%d]_b=", ind);
    dclistPrint(stderr, SWAR->S[ind]->next); fprintf(stderr, "\n");
#endif
    SWAR->A[j] = dclistInsert(SWAR->S[ind], j);
#if DEBUG >= 2
    fprintf(stderr, "S[%d]_a=", ind);
    dclistPrint(stderr, SWAR->S[ind]->next); fprintf(stderr, "\n");
#endif
    // update R[j] by removing i
    remove_i_from_Rj(SWAR, i, j);
}

// for all j in row[i], removes j and update data
void
removeRowFast(sparse_mat_t *mat, int i, swar_t *SWAR)
{
    int k;

    for(k = 0; k < mat->data[i].len; k++)
	removeCellFast(mat, i, mat->data[i].val[k], SWAR);
}

void
addCellFast(sparse_mat_t *mat, int i, int j, swar_t *SWAR)
{
    int ind;

    // update weight
#if DEBUG >= 1
    fprintf(stderr, "Moving j=%d from S[%d] to S[%d]\n",
	    j, SWAR->Wj[j], SWAR->Wj[j]+1);
#endif
    ind = SWAR->Wj[j] = incrS(SWAR->Wj[j]);
#if USE_MERGE_FAST > 1
    if(ind < 0)
	return;
#endif
    remove_j_from_S(SWAR, j);
    if(SWAR->Wj[j] > SWAR->cwmax){
#if DEBUG >= 1
	fprintf(stderr, "WARNING: column %d is too heavy (%d)\n", j,
		SWAR->Wj[j]);
#endif
	ind = SWAR->cwmax+1; // trick
    }
    // update A[j]
    SWAR->A[j] = dclistInsert(SWAR->S[ind], j);
    // update R[j] by adding i
    add_i_to_Rj(SWAR, i, j);
}

void
addRowFast(sparse_mat_t *mat, int i, swar_t *SWAR)
{
    int k;

    for(k = 0; k < mat->data[i].len; k++)
	addCellFast(mat, i, mat->data[i].val[k], SWAR);
}

// i1 += i2.
void
addRowsFast(sparse_mat_t *mat, int i1, int i2, swar_t *SWAR)
{
    int k1, k2, k, len, *tmp, *tmp2;

#if TEX
    fprintf(stderr, "row[%d] += row[%d]\n", i1, i2);
#endif
#if DEBUG >= 1
    fprintf(stderr, "row[%d] =", i1); print_row(mat, i1);
    fprintf(stderr, "\n");
    fprintf(stderr, "row[%d] =", i2); print_row(mat, i2);
    fprintf(stderr, "\n");
#endif
    // merge row[i1] and row[i2] in i1...
    len = mat->data[i1].len + mat->data[i2].len;
    tmp = (int *)malloc(len * sizeof(tmp));
    k = k1 = k2 = 0;

    // i1 is to disappear, replaced by a new one
    removeRowFast(mat, i1, SWAR);
    // loop while everybody is here
    while((k1 < mat->data[i1].len) && (k2 < mat->data[i2].len)){
	if(mat->data[i1].val[k1] < mat->data[i2].val[k2])
	    tmp[k++] = mat->data[i1].val[k1++];
	else if(mat->data[i1].val[k1] > mat->data[i2].val[k2])
            tmp[k++] = mat->data[i2].val[k2++];
	else{
	    k1++;
	    k2++;
	}
    }
    // finish with k1
    for( ; k1 < mat->data[i1].len; k1++)
	tmp[k++] = mat->data[i1].val[k1];
    // finish with k2
    for( ; k2 < mat->data[i2].len; k2++)
	tmp[k++] = mat->data[i2].val[k2];
    // destroy and copy back
    free(mat->data[i1].val);
    tmp2 = (int *)malloc(k * sizeof(int));
    memcpy(tmp2, tmp, k * sizeof(int));
    mat->data[i1].len = k;
    mat->data[i1].val = tmp2;
    free(tmp);
    addRowFast(mat, i1, SWAR);
#if DEBUG >= 1
    fprintf(stderr, "row[%d]+row[%d] =", i1, i2);
    print_row(mat, i1); fprintf(stderr, "\n");
#endif
}

void
remove_j_from_row(sparse_mat_t *mat, int i, int j)
{
    int k;

#if DEBUG >= 2
    fprintf(stderr, "row[%d]_b=", i);
    print_row(mat, i);
    fprintf(stderr, "\n");
#endif
    for(k = 0; k < mat->data[i].len; k++)
	if(mat->data[i].val[k] == j)
	    break;
    // crunch
    for(++k; k < mat->data[i].len; k++)
	mat->data[i].val[k-1] = mat->data[i].val[k];
    mat->data[i].len -= 1;
#if DEBUG >= 2
    fprintf(stderr, "row[%d]_a=", i);
    print_row(mat, i);
    fprintf(stderr, "\n");
#endif
}

// These columns are simply removed from the current squeleton matrix, but
// not from the small matrix.
int
deleteAllColsFromStack(sparse_mat_t *mat, swar_t *SWAR, int iS)
{
    dclist dcl;
    int j, k, njrem = 0;

    while(1){
	dcl = SWAR->S[iS]->next;
	if(dcl == NULL)
	    break;
	j = dcl->j;
	njrem++;
#if DEBUG >= 1
	fprintf(stderr, "Removing column %d from S[%d]\n", j, iS);
#endif
	// we destroy column j
	// TODO: do faster?
#if USE_MERGE_FAST > 1
	if(iS == 1)
#endif
	for(k = 1; k <= SWAR->R[j][0]; k++)
	    if(SWAR->R[j][k] != -1)
		remove_j_from_row(mat, SWAR->R[j][k], j);
#if USE_MERGE_FAST > 1
	k = SWAR->Wj[j]; // make a copy of the weight
#endif
	remove_j_from_SWAR(SWAR, j);
	// SWAR->Wj[j] is put to -1...
#if USE_MERGE_FAST > 1
	SWAR->Wj[j] = -k-1; // restore and update
#endif
    }
#if DEBUG >= 1
    if(njrem > 0)
	fprintf(stderr, "deleteAllColsFromStack[%d]: %d\n", iS, njrem);
#endif
    return njrem;
}

int
deleteHeavyColumns(sparse_mat_t *mat, swar_t *SWAR)
{
    return deleteAllColsFromStack(mat, SWAR, SWAR->cwmax+1);
}

int
removeSingletons(sparse_mat_t *mat, swar_t *SWAR)
{
#if USE_MERGE_FAST == 0
    int i, j;

    for(j = 0; j < mat->ncols; j++)
	if(mat->wt[j] == 1){
	    // find row...
	    for(i = 0; i < mat->nrows; i++)
		if(hasCol(mat, i, j))
		    break;
	    assert(i < mat->nrows);
	    removeWeight(mat, i);
	    destroyRow(mat, i);
	    report1(i); // signal to replay...!
	}
#else
    return deleteAllColsFromStack(mat, SWAR, 1);
#endif
}

// try all combinations to find the smaller one
void
tryAllCombinations(sparse_mat_t *mat, int m, swar_t *SWAR, int *ind)
{
    int i, k;
    
    i = findBestIndex(mat, m, ind);
#if DEBUG >= 1
    fprintf(stderr, "Minimal is i=%d\n", i);
#endif
    for(k = 0; k < m; k++)
	if(k != i){
	    addRowsFast(mat, ind[k], ind[i], SWAR);
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
    reportn(ind, m);
    removeRowFast(mat, ind[0], SWAR);
    destroyRow(mat, ind[0]);
}

void
useMinimalSpanningTree(sparse_mat_t *mat, int m, swar_t *SWAR, int *ind)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX];
    int i, father[MERGE_LEVEL_MAX], height[MERGE_LEVEL_MAX], hmax, nV = m, h;
    int adds[MERGE_LEVEL_MAX << 1], nadds, tab[MERGE_LEVEL_MAX << 1], itab;

    fillRowAddMatrix(A, mat, m, ind);
    hmax = minimalSpanningTree(father, height, mat, m, SWAR, A);
#if DEBUG >= 1
    for(i = 0; i < m; i++)
	fprintf(stderr, "father[%d] = %d\n", i, father[i]);
#endif
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

		addRowsFast(mat, i1, i2, SWAR);

		// TODO: use the fact we know w(i1, i2) = A[][] to compute
		// the length of row[i1]+row[i2] directly

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
		reportn(tab, itab);
		tab[0] = adds[i];
		tab[1] = adds[i+1];
		itab = 2;
	    }
	}
	if(h == 1)
	    // we are at the root, actually
	    tab[0] = -tab[0]-1; // what a trick!
	reportn(tab, itab);
    }
#if 0 // for non compact reporting
    // delete root!
    report1(ind[0]);
#endif
    removeRowFast(mat, ind[0], SWAR);
    destroyRow(mat, ind[0]);
}

void
findOptimalCombination(sparse_mat_t *mat, int m, swar_t *SWAR, int *ind)
{
    if(m <= 3)
	tryAllCombinations(mat, m, SWAR, ind);
    else{
#if 0
	tryAllCombinations(mat, m, SWAR, ind);
#else
	useMinimalSpanningTree(mat, m, SWAR, ind);
#endif
    }
}

int
merge_m_fast(sparse_mat_t *mat, int m, swar_t *SWAR)
{
    int totfindrows = 0, totfindbest = 0, totadd = 0;
    int *ind, j, k, njproc = 0, ni, njrem = 0;
    dclist dcl = SWAR->S[m];

#if DEBUG >= 1
    fprintf(stderr, "Weight %d:", m);
#endif
    ind = (int *)malloc(m * sizeof(int));
    while(1){
	// get next j
	dcl = SWAR->S[m]->next;
	if(dcl == NULL)
	    break;
	j = dcl->j;
	njproc++;
	if(!(njproc % 1000))
	    fprintf(stderr, "# %d columns of weight %d processed\n",njproc,m);
#if DEBUG >= 1
	fprintf(stderr, "Treating %d-th column %d\n", njproc, j);
#endif
#if (DEBUG >= 1) || TEX
	fprintf(stderr, "Status before next j=%d to start\n", j);
# if TEX
	texSWAR(mat, SWAR);
# else
	printSWAR(SWAR, mat->ncols);
# endif
#endif
	// the corresponding rows are in R[j], skipping 1st cell and -1's
	ni = 0;
	for(k = 1; k <= SWAR->R[j][0]; k++){
	    if(SWAR->R[j][k] != -1){
		ind[ni++] = SWAR->R[j][k];
		if(ni == m)
		    break;
	    }
	}
#if DEBUG >= 1
	fprintf(stderr, " %d", j);
	fprintf(stderr, "=> the %d rows are:\n", j);
	for(k = 0; k < m; k++){
	    fprintf(stderr, "row[%d]=", ind[k]);
	    print_row(mat, ind[k]);
	    fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
#endif
	findOptimalCombination(mat, m, SWAR, ind);
	mat->rem_nrows--;
	mat->rem_ncols--;
	remove_j_from_SWAR(SWAR, j);
	njrem += deleteHeavyColumns(mat, SWAR);
    }
    fprintf(stderr, "TIME: findrows=%d findbest=%d add=%d\n",
	    totfindrows, totfindbest, totadd);
    free(ind);
#if DEBUG >= 1
    fprintf(stderr, "Status at the end of the process\n");
    texSWAR(mat, SWAR);
#endif
    return njrem;
}

int
merge_m(sparse_mat_t *mat, int m, swar_t *SWAR)
{
#if USE_MERGE_FAST >= 1
    return merge_m_fast(mat, m, SWAR);
#else
    return merge_m_slow(mat, m);
#endif
}

// TODO: use mergemax?
int
mergeGe2(sparse_mat_t *mat, int m, int nb_merge_max, swar_t *SWAR)
{
#if USE_MERGE_FAST == 0
    int j, nbm;

    // TODO: remove this and start directly merge_m?
    for(j = 0, nbm = 0; j < mat->ncols; j++)
	if(mat->wt[j] == m){
# if DEBUG >= 0
	    fprintf(stderr, "# wt[%d] = %d\n", j, m);
# endif
	    nbm++;
	}
    fprintf(stderr, "There are %d column(s) of weight %d\n", nbm, m);
    if(nbm)
#endif
	return merge_m(mat, m, SWAR);
#if DEBUG >= 1
    matrix2tex(mat);
#endif
}

int
minColWeight(sparse_mat_t *mat, swar_t *SWAR)
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

    for(m = 1; m <= SWAR->cwmax; m++)
	if(SWAR->S[m]->next != NULL)
	    return m;
    return -1;
#endif
}

void
inspectRowWeight(sparse_mat_t *mat, int rwmax)
{
    int i, maxw = 0;

    for(i = 0; i < mat->nrows; i++)
	if(mat->data[i].val != NULL){
	    if(mat->data[i].len > rwmax){
		fprintf(stderr, "Too heavy row[%d]: %d\n", i, mat->data[i].len);
		// TODO: destroy row...!
	    }
	    if(mat->data[i].len > maxw)
		maxw = mat->data[i].len;
	}
    fprintf(stderr, "Max row weight is %d\n", maxw);
}

void
merge(sparse_mat_t *mat, int nb_merge_max, int maxlevel, int rwmax, swar_t *SWAR)
{
    int old_nrows, old_ncols, m, mm, njrem = 0;

    report2(mat->nrows, mat->ncols);
    m = 2;
    while(1){
	fprintf(stderr, "Performing merge %d\n", m);
	old_nrows = mat->rem_nrows;
	old_ncols = mat->rem_ncols;
	if(m == 1)
	    njrem += removeSingletons(mat, SWAR);
#if 0
	else if(m == 2)
	    merge2(mat, nb_merge_max);
#endif
	else
	    njrem += mergeGe2(mat, m, nb_merge_max, SWAR);
#if 0
    {
	int i, ni = 0;
	for(i = 0; i < mat->nrows; i++)
	    if((mat->data[i].val != NULL) && hasCol(mat, i, 301))
		ni++;
	printf("CHECK: %d %d\n", ni, SWAR->Wj[301]);
    }
#endif
	fprintf(stderr, "=> nrows=%d ncols=%d njrem=%d\n",
		mat->rem_nrows, mat->rem_ncols, njrem);
	inspectRowWeight(mat, rwmax);
	mm = minColWeight(mat, SWAR);
	fprintf(stderr, "Min col weight = %d\n", mm);
	if((old_nrows == mat->rem_nrows) && (old_ncols == mat->rem_ncols)){
	    // nothing happened this time and mm > m
	    m = mm;
	    if((m > maxlevel) || (m <= 0))
		break;
	}
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

// not activated anymore...!
void
dumpSparse(FILE *ofile, sparse_mat_t *mat, int *code)
{
    int i, j, k, buf[1000], ibuf, new_nrows = 0;

    fprintf(ofile, "%d %d\n", mat->rem_nrows, mat->rem_ncols);
    for(i = 0; i < mat->nrows; i++){
	if(mat->data[i].val == NULL)
	    continue;
	new_nrows++;
	ibuf = 0;
	for(k = 0; k < mat->data[i].len; k++){
	    j = mat->data[i].val[k];
	    if(code[j])
		buf[ibuf++] = code[j]-1;
	}
	fprintf(ofile, "%d", ibuf);
	for(k = 0; k < ibuf; k++)
	    fprintf(ofile, " %d", buf[k]);
	fprintf(ofile, "\n");
    }
    assert(new_nrows == mat->rem_nrows);
}

int
main (int argc, char *argv[])
{
  FILE *purgedfile;
  sparse_mat_t mat;
  char *purgedname = NULL, *hisname = NULL;
  int nb_merge_max = 0, *usecol, nrows, ncols;
  int cwmax = 20, rwmax = 1000000, maxlevel = 2;
  swar_t SWAR;

#if TEX
  fprintf(stderr, "\\begin{verbatim}\n");
#endif    
  fprintf (stderr, "%s rev. %s\n", argv[0], REV);

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 2 && strcmp (argv[1], "-merge") == 0)
	{
	  nb_merge_max = atoi(argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-mat") == 0)
	{
	  purgedname = argv[2];
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-rebuild") == 0)
	{
	  hisname = argv[2];
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-cwmax") == 0)
	{
	  cwmax = atoi(argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-rwmax") == 0)
	{
	  rwmax = atoi(argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-maxlevel") == 0)
	{
	  maxlevel = atoi(argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else 
	break;
    }

  purgedfile = fopen(purgedname, "r");
  assert(purgedfile != NULL);
  fscanf(purgedfile, "%d %d", &nrows, &ncols);
  usecol = (int *)malloc(ncols *sizeof(int));
  inspectWeight(usecol, purgedfile, nrows, ncols, cwmax);

  initSWAR(&SWAR, ncols, cwmax, maxlevel);
  fillSWAR(&SWAR, usecol, ncols);
    
  rewind(purgedfile);
  readmat(&SWAR, purgedfile, &mat, usecol);
  fclose(purgedfile);
#if DEBUG >= 3
  checkmat(&mat);
#endif
#if USE_USED >= 1
  used = (char *)malloc(mat.nrows * sizeof(char));
  memset(used, 0, mat.nrows * sizeof(char));
#endif
  merge(&mat, nb_merge_max, maxlevel, rwmax, &SWAR);
#if TEX
  fprintf(stderr, "\\end{verbatim}\n");
#endif
  //  checkData(&mat, &SWAR);
  closeSWAR(&SWAR);
  return 0;
}
