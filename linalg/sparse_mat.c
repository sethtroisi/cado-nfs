#include <string.h>

#include "utils/utils.h"
#include "files.h"
#include "merge_opts.h"
#include "sparse.h"
#ifndef USE_MARKOWITZ
# include "dclist.h"
#endif
#include "sparse_mat.h"
#ifndef USE_MARKOWITZ
# include "swar.h"
#endif

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

// TODO: ncols could be a new renumbered ncols...
void
initMat(sparse_mat_t *mat, INT jmin, INT jmax)
{
#if USE_MERGE_FAST
    INT **R;
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
Sparse2Tex(sparse_mat_t *mat)
{
    int j;
    INT k;

    fprintf(stderr, "\\end{verbatim}\n");
    matrix2tex(mat);
#ifndef USE_MARKOWITZ
    texSWAR(mat);
#endif
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

#if USE_TAB == 0
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
#endif

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

void
freeRj(sparse_mat_t *mat, int j)
{
#ifndef USE_COMPACT_R
    free(mat->R[GETJ(mat, j)]);
    mat->R[GETJ(mat, j)] = NULL;
#else
    fprintf(stderr, "R: NYI in freeRj\n");
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

