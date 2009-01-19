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

#define DEBUG 0

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
    mat->rows = (INT **)malloc(mat->nrows * sizeof (INT *));
    ASSERT_ALWAYS(mat->rows != NULL);
    mat->wt = (int*) malloc ((mat->jmax-mat->jmin) * sizeof (int));
    memset(mat->wt, 0, (mat->jmax-mat->jmin) * sizeof (int));
    mat->wburried = (int*) malloc (mat->nrows * sizeof (int));
    memset(mat->wburried, 0, mat->nrows * sizeof (int));
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
initWeightFromFile(sparse_mat_t *mat, FILE *purgedfile, int skipfirst)
{
    INT j, jmin = mat->jmin, jmax = mat->jmax;
    int i, k, ret, nc;

    memset(mat->wt, 0, (mat->jmax - mat->jmin) * sizeof(int));
    for(i = 0; i < mat->nrows; i++){
	if(skipfirst != 0){
	    ret = fscanf(purgedfile, "%d", &k); // unused index to rels file
	    ASSERT_ALWAYS (ret == 1);
	}
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

// builds Rj[j] for light j's.
void
fillmat(sparse_mat_t *mat)
{
    INT j, *Rj, jmin = mat->jmin, jmax = mat->jmax;

    for(j = jmin; j < jmax; j++){
#  if DEBUG >= 2
	fprintf(stderr, "Treating column %d\n", j);
#  endif
	if(mat->wt[GETJ(mat, j)] <= mat->cwmax){
#ifndef USE_MARKOWITZ
	    mat->A[GETJ(mat, j)] = dclistInsert(mat->S[mat->wt[GETJ(mat, j)]], j);
#  if DEBUG >= 2
	    fprintf(stderr, "Inserting %d in S[%d]:", j, mat->wt[GETJ(mat, j)]);
	    dclistPrint(stderr, mat->S[mat->wt[GETJ(mat, j)]]->next);
	    fprintf(stderr, "\n");
#  endif
#endif // USE_MARKOWITZ
#ifndef USE_COMPACT_R
	    Rj = (INT *)malloc((mat->wt[GETJ(mat, j)]+1) * sizeof(INT));
	    Rj[0] = 0; // last index used
	    mat->R[GETJ(mat, j)] = Rj;
#else
	    fprintf(stderr, "R: NYI in fillSWAR\n");
	    exit(1);
#endif
	}
	else{
#if USE_MERGE_FAST <= 1
	    mat->wt[GETJ(mat, j)] = -1;
#else
	    mat->wt[GETJ(mat, j)] = -mat->wt[GETJ(mat, j)]; // trick!!!
#endif
#ifndef USE_MARKOWITZ
	    mat->A[GETJ(mat, j)] = NULL; // TODO: renumber j's?????
#endif
#ifndef USE_COMPACT_R
	    mat->R[GETJ(mat, j)] = NULL;
#else
            fprintf(stderr, "R: NYI2 in fillSWAR\n");
            exit(1);
#endif
	}
    }
}

int 
cmp(const void *p, const void *q) {
    int x = *((int *)p);
    int y = *((int *)q);
    return (x <= y ? -1 : 1);
}

/* Reads a matrix file, and puts in mat->wt[j], 0 <= j < ncols, the
   weight of column j (adapted from matsort.c).

   We skip columns that are too heavy.

   mat->wt is already filled and put to - values when needed. So don't touch.

   mat->weight is correct on exit; it will be only approximately true
   when skipheavycols will be activated.
*/
int
readmat(sparse_mat_t *mat, FILE *file, int skipfirst, int skipheavycols)
{
    int ret;
    int i, k;
    int nc, lbuf = 100, *buf;
#if (USE_MERGE_FAST < 3) && !defined(USE_MPI)
    int nh = 0;
#endif
    INT ibuf, j, jmin = mat->jmin, jmax = mat->jmax;
    char *tooheavy = NULL;
    int nburried = 0;

    buf = (int *)malloc(lbuf * sizeof(int));
    ret = fscanf (file, "%d %d", &i, &k); // already set up...!
    ASSERT_ALWAYS (ret == 2);
    
    fprintf(stderr, "Reading matrix of %d rows and %d columns: excess is %d\n",
	    mat->nrows, mat->ncols, mat->nrows - mat->ncols);
    mat->rem_nrows = mat->nrows;
    mat->weight = 0;

    if(skipheavycols != 0){
	int bmin = mat->nrows, bmax = 0, wc, wmax;

	wmax = 5 * mat->cwmax;
	// heavy columns already have wt < 0
	tooheavy = (char *)malloc(mat->ncols * sizeof(char));
	memset(tooheavy, 0, mat->ncols * sizeof(char));
	for(j = 0; j < mat->ncols; j++){
	    if((wc = abs(mat->wt[j])) > wmax){
#if DEBUG >= 1 
		fprintf(stderr, "Burrying j=%d (wt=%d)\n", j, wc);
#endif
		tooheavy[j] = 1;
		nburried++;
		if(wc > bmax) bmax = wc;
		if(wc < bmin) bmin = wc;
	    }
	}
	fprintf(stderr, "# Number of burried columns is %d", nburried);
	fprintf(stderr, " (min=%d max=%d)\n", bmin, bmax);
    }

    for (i = 0; i < mat->nrows; i++){
	if(!(i % 100000))
	    fprintf(stderr, "Reading %d-th row\n", i);
	if(skipfirst != 0){
	    ret = fscanf(file, "%d", &nc); // unused index to rels file
	    if(ret != 1){
		fprintf(stderr, "Pb1 in readmat\n");
		return 0;
	    }
	}
	ret = fscanf (file, "%d", &nc);
	if(ret != 1){
	    fprintf(stderr, "Pb2 in readmat\n");
	    return 0;
	}
	if(nc == 0){
	    mat->rows[i] = NULL;
	    mat->rem_nrows--;
	}
	else{
#if (USE_MERGE_FAST < 3) && !defined(USE_MPI)
	    int nb_heavy_j = 0;
#endif
	    if(nc > lbuf){
		lbuf <<= 1;
		fprintf(stderr, "Warning: doubling lbuf in readmat;");
		fprintf(stderr, " new value is %d\n", lbuf);
		free(buf);
		buf = (int *)malloc(lbuf * sizeof(int));
	    }
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
		    // this will be the weight in the current slice
		    mat->weight++; 
		    if((tooheavy != NULL) && (tooheavy[j] != 0))
			mat->wburried[i] += 1;
		    else{
			buf[ibuf++] = j;
			if(mat->wt[GETJ(mat, j)] > 0){ // redundant test?
#ifndef USE_COMPACT_R
			    mat->R[GETJ(mat, j)][0]++;
			    mat->R[GETJ(mat, j)][mat->R[GETJ(mat, j)][0]] = i;
#else
			    fprintf(stderr, "R: NYI in readmat\n");
			    exit(1);
#endif
			}
		    }
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
	    // TODO: do not store rows not having at least one light
	    // column, but do not decrease mat->nrows!
	    mat->rows[i] = (INT*) malloc ((ibuf+1) * sizeof (INT));
	    mat->rows[i][0] = ibuf;
	    memcpy(mat->rows[i]+1, buf, ibuf * sizeof(INT));
	    // sort indices in val to ease row merges
            qsort(mat->rows[i]+1, ibuf, sizeof(INT), cmp);
            /* check all indices are distinct, otherwise this is a bug of
               purge */
            for (k = 1; k < ibuf; k++)
              if (mat->rows[i][k] == mat->rows[i][k+1])
                {
                  fprintf (stderr, "Error, duplicate ideal %x in row %i\n",
                           mat->rows[i][k], i);
                  exit (1);
                }
	}
    }
#if (USE_MERGE_FAST < 3) && !defined(USE_MPI)
    fprintf(stderr, "Number of heavy rows: %d\n", nh);
#endif
    // we need to keep informed of what really happens; this will be an upper
    // bound on the number of active columns, I guess
    mat->rem_ncols = mat->ncols;
#ifndef USE_MARKOWITZ
    printStatsSWAR(mat);
#endif
    free(buf);
    if(tooheavy != NULL)
	free(tooheavy);
    return 1;
}

void
print_row(sparse_mat_t *mat, int i)
{
    fprintRow(stderr, mat->rows[i]);
}

void
destroyRow(sparse_mat_t *mat, int i)
{
    free(mat->rows[i]);
    mat->rows[i] = NULL;
}

// we do not use/touch wburried[i]
void
removeWeightFromRow(sparse_mat_t *mat, int i)
{
    INT k;

    for(k = 1; k <= mat->rows[i][0]; k++)
	mat->wt[GETJ(mat, mat->rows[i][k])]--;
#if TRACE_ROW >= 0
    if(i == TRACE_ROW)
	fprintf(stderr, "TRACE_ROW: removeWeightFromRow %d\n", 
		lengthRow(mat, i));
#endif
    mat->weight -= lengthRow(mat, i);
}

// we do not use/touch wburried[i]
void
addWeightFromRow(sparse_mat_t *mat, int i)
{
    INT k;

    for(k = 1; k <= mat->rows[i][0]; k++)
	mat->wt[GETJ(mat, mat->rows[i][k])]++;
    mat->weight += lengthRow(mat, i);
}

// i1 += i2, mat->wt is updated at the same time.
void
addRowsWithWeight(sparse_mat_t *mat, int i1, int i2)
{
    // i1 is to disappear, replaced by a new one
    removeWeightFromRow(mat, i1);
    addRows(mat->rows, i1, i2, -1);
    mat->wburried[i1] += mat->wburried[i2];
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
// New: uses the estimated wburried stuff. TODO: what about MPI???
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
    return w + mat->wburried[i1] + mat->wburried[i2];
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

void
fillTabWithRowsForGivenj(INT *ind, sparse_mat_t *mat, INT j)
{
    int ni = 0, k, i;

    for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
	if((i = mat->R[GETJ(mat, j)][k]) != -1){
	    ind[ni++] = i;
	    if(ni == mat->wt[GETJ(mat, j)])
		break;
	}
}
