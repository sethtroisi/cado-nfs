#include <string.h>

#include "utils.h"
#include "merge_opts.h"
#include "sparse.h"
#ifndef USE_MARKOWITZ
# include "dclist.h"
#endif
#include "filter_matrix.h"
#ifndef USE_MARKOWITZ
# include "swar.h"
#endif

#define DEBUG 0

int
decrS (int w)
{
  return (w >= 0 ? w - 1 : w + 1);
}

int
incrS (int w)
{
  return (w >= 0 ? w + 1 : w - 1);
}

/* initialize the sparse matrix mat,
   where we consider only columns from jmin to jmax-1 */
// TODO: ncols could be a new renumbered ncols...
void
initMat (filter_matrix_t *mat, int32_t jmin, int32_t jmax)
{
    int32_t **R;

    mat->jmin = jmin;
    mat->jmax = jmax;
    mat->rows = (int32_t **) malloc (mat->nrows * sizeof (int32_t *));
    ASSERT_ALWAYS (mat->rows != NULL);
    mat->wt = (int*) malloc ((mat->jmax - mat->jmin) * sizeof (int));
    memset (mat->wt, 0, (mat->jmax - mat->jmin) * sizeof (int));
    mat->wburied = (int*) malloc (mat->nrows * sizeof (int));
    memset (mat->wburied, 0, mat->nrows * sizeof (int));
    R = (int32_t **) malloc ((mat->jmax - mat->jmin) * sizeof(int32_t *));
    ASSERT_ALWAYS(R != NULL);
    mat->R = R;
}

void
clearMat (filter_matrix_t *mat)
{
  int32_t i, j;

  for (i = 0; i < mat->nrows; i++)
    free (mat->rows[i]);
  free (mat->rows);
  free (mat->wt);
  free (mat->wburied);
  for (j = 0; j < mat->jmax - mat->jmin; j++)
    free (mat->R[j]);
  free (mat->R);
}

void
checkData(filter_matrix_t *mat)
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


#define MAX_LENGTH 1024

/* compute the weight mat->wt[j] of each column j, for the set of relations
   in file purgedfile.
   Assume mat->wt is initialized to 0 (as done by initMat).
   Assume we have already read the 1st line of purgedfile (nrows, ncols).
   If skipfirst is non-zero, skip the 1st entry of each line.
*/
void
filter_matrix_read_weights(filter_matrix_t * mat, purgedfile_stream_ptr ps)
{
    int32_t jmin = mat->jmin, jmax = mat->jmax;
    for (; purgedfile_stream_get(ps,NULL) >= 0;) {
	for (int k = 0; k < ps->nc; k++) {
	    int32_t j = ps->cols[k];
	    if (LIKELY((j >= jmin) && (j < jmax)))
		mat->wt[GETJ(mat, j)]++;
	}
    }
}

/* initialize Rj[j] for light columns, i.e., for those of weight <= cwmax */
void
fillmat (filter_matrix_t *mat)
{
  int32_t j, *Rj, jmin = mat->jmin, jmax = mat->jmax, wj;

    for(j = jmin; j < jmax; j++){
        wj = mat->wt[GETJ(mat, j)];
	if (wj <= mat->cwmax){
#ifndef USE_MARKOWITZ
            mat->A[GETJ(mat, j)] = dclistInsert(mat->S[wj], j);
#endif
#ifndef USE_COMPACT_R
	    Rj = (int32_t *) malloc((wj + 1) * sizeof(int32_t));
	    Rj[0] = 0; /* last index used */
	    mat->R[GETJ(mat, j)] = Rj;
#else
	    fprintf(stderr, "R: NYI in fillSWAR\n");
	    exit(1);
#endif
	}
	else{ /* weight is larger than cwmax */
	    mat->wt[GETJ(mat, j)] = -mat->wt[GETJ(mat, j)]; // trick!!!
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

/* Reads a matrix file.

   We skip columns that are too heavy.

   mat->wt is already filled:
   * if mat->wt[j] > 0, it is the weight of column j that we consider
   * if mat->wt[j] < 0, we don't consider column j

   mat->weight is correct on exit; it will be only approximately true
   when skipheavycols will be activated.
*/
int
filter_matrix_read (filter_matrix_t *mat, purgedfile_stream_ptr ps, int verbose)
{
    int lbuf = 100, *buf;
#if !defined(USE_MPI)
    int nh = 0;
#endif
    int32_t ibuf, j, jmin = mat->jmin, jmax = mat->jmax;
    char *tooheavy = NULL;

    mat->nburied = 0;
    buf = (int *) malloc(lbuf * sizeof(int));

    fprintf(stderr, "Reading matrix of %d rows and %d columns: excess is %d\n",
            mat->nrows, mat->ncols, mat->nrows - mat->ncols);
    mat->rem_nrows = mat->nrows;
    mat->weight = 0;

    int bmin = mat->nrows, bmax = 0, wmax;

    /* don't consider heavy columns of density > 1/2000, this roughly
       corresponds to primes < 2000 or ideals of norm < 2000 */
    wmax = -(mat->nrows / 2000);
    /* heavy columns already have wt < 0 */
    tooheavy = (char *) malloc (mat->ncols * sizeof(char));
    memset (tooheavy, 0, mat->ncols * sizeof(char));
    for(j = 0; j < mat->ncols; j++){
        if(mat->wt[j] < wmax){
            int wc = -mat->wt[j];
#if DEBUG >= 1
            fprintf(stderr, "Burying j=%d (wt=%d)\n", j, wc);
#endif
            tooheavy[j] = 1;
            mat->nburied += 1;
            if(wc > bmax) bmax = wc;
            if(wc < bmin) bmin = wc;
        }
    }
    fprintf(stderr, "# Number of buried columns is %d", mat->nburied);
    fprintf(stderr, " (min=%d max=%d)\n", bmin, bmax);

    for (int i = 0 ; purgedfile_stream_get(ps, NULL) >= 0 ; i++) {
        ASSERT_ALWAYS(i < mat->nrows);
        if (purgedfile_stream_disp_progress_now_p(ps))
            fprintf(stderr, "read %d rows in %.1f s (%.1f MB/s, %.1f rows/s)\n",
                    ps->rrows, ps->dt, ps->mb_s, ps->rows_s);

        if(ps->nc == 0) {
            mat->rows[i] = NULL;
            mat->rem_nrows--;
        } else{
#if !defined(USE_MPI)
            int nb_heavy_j = 0;
#endif
            if(ps->nc > lbuf){
                lbuf <<= 1;
                fprintf(stderr, "Warning: doubling lbuf in readmat;");
                fprintf(stderr, " new value is %d\n", lbuf);
                free(buf);
                buf = (int *)malloc(lbuf * sizeof(int));
            }
            ibuf = 0;
            for(int k = 0 ; k < ps->nc; k++) {
                int32_t j = ps->cols[k];
                ASSERT_ALWAYS (0 <= j && j < mat->ncols);
#if !defined(USE_MPI)
                if(mat->wt[GETJ(mat, j)] < 0)
                    nb_heavy_j++;
#endif
                // always store j in the right interval...!
                if((j >= jmin) && (j < jmax)){
                    // this will be the weight in the current slice
                    mat->weight++; 
                    if ((tooheavy != NULL) && (tooheavy[j] != 0)) {
                        mat->wburied[i] += 1;
                    } else {
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
            }
#if !defined(USE_MPI)
            if(nb_heavy_j == ps->nc){
                // all the columns are heavy and thus will never participate
                mat->rows[i] = NULL;
                nh++;
                continue;
            }
#endif
            // TODO: do not store rows not having at least one light
            // column, but do not decrease mat->nrows!
            mat->rows[i] = (int32_t*) malloc ((ibuf+1) * sizeof (int32_t));
            mat->rows[i][0] = ibuf;
            memcpy(mat->rows[i]+1, buf, ibuf * sizeof(int32_t));
            // sort indices in val to ease row merges
            qsort(mat->rows[i]+1, ibuf, sizeof(int32_t), cmp);
            /* check all indices are distinct, otherwise this is a bug of
               purge */
            for (int k = 1; k < ibuf; k++)
                if (mat->rows[i][k] == mat->rows[i][k+1])
                {
                    fprintf (stderr, "Error, duplicate ideal %x in row %i\n",
                            mat->rows[i][k], i);
                    exit (1);
                }
        }
    }
    fprintf (stderr, "\n"); /* to keep last output on screen */
#if !defined(USE_MPI)
    fprintf(stderr, "Number of heavy rows: %d\n", nh);
#endif
    // we need to keep informed of what really happens; this will be an upper
    // bound on the number of active columns, I guess
    mat->rem_ncols = mat->ncols;
#ifndef USE_MARKOWITZ
    if (verbose)
        printStatsSWAR (mat);
#endif
    free (buf);
    if (tooheavy != NULL)
        free (tooheavy);
    return 1;
}

void
print_row(filter_matrix_t *mat, int i)
{
    fprintRow(stderr, mat->rows[i]);
}

void
destroyRow(filter_matrix_t *mat, int i)
{
    free(mat->rows[i]);
    mat->rows[i] = NULL;
}

void
freeRj(filter_matrix_t *mat, int j)
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
remove_i_from_Rj(filter_matrix_t *mat, int i, int j)
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
add_i_to_Rj(filter_matrix_t *mat, int i, int j)
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
	mat->R[GETJ(mat, j)] = (int32_t *)realloc(mat->R[GETJ(mat, j)], l * sizeof(int32_t));
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
remove_j_from_row(filter_matrix_t *mat, int i, int j)
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
// New: uses the estimated wburied stuff. TODO: what about MPI???
int
weightSum(filter_matrix_t *mat, int i1, int i2)
{
    int k1, k2, w, len1, len2;

    w = mat->wburied[i1] + mat->wburied[i2];
    if(w > mat->nburied)
	w = mat->nburied;
    len1 = (isRowNull(mat, i1) ? 0 : lengthRow(mat, i1));
    len2 = (isRowNull(mat, i2) ? 0 : lengthRow(mat, i2));
    if((len1 == 0) || (len2 == 0))
	fprintf(stderr, "i1=%d i2=%d len1=%d len2=%d\n", i1, i2, len1, len2);
    k1 = k2 = 1;
    while((k1 <= len1) && (k2 <= len2)){
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

void
fillTabWithRowsForGivenj(int32_t *ind, filter_matrix_t *mat, int32_t j)
{
    int ni = 0, k, i;

    for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
	if((i = mat->R[GETJ(mat, j)][k]) != -1){
	    ind[ni++] = i;
	    if(ni == mat->wt[GETJ(mat, j)])
		break;
	}
}
