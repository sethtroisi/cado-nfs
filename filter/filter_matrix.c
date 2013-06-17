#include "cado.h"

#include <string.h>
#include <stdio.h>

#include "portability.h"
#include "utils.h"
#include "merge_opts.h"
#include "filter_matrix.h"
#include "sparse.h"

#define DEBUG 0

unsigned int weight_ffs (int e)
{
  if (e == 0)
      return 0;
  else
      return 1; /* Should depend on e, for now jsut constant*/
}

int
decrS (int w)
{
  return (w >= 0 ? w - 1 : w + 1);
}

/* increment the absolute value of a signed number */
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
    mat->rows = (typerow_t **) malloc (mat->nrows * sizeof (typerow_t *));
    ASSERT_ALWAYS (mat->rows != NULL);
    mat->wt = (int*) malloc ((mat->jmax - mat->jmin) * sizeof (int));
    memset (mat->wt, 0, (mat->jmax - mat->jmin) * sizeof (int));
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
*/
void
filter_matrix_read_weights(filter_matrix_t * mat, purgedfile_stream_ptr ps)
{
    int32_t jmin = mat->jmin, jmax = mat->jmax;
    for (; purgedfile_stream_get(ps,NULL) >= 0;) 
      {
#ifdef FOR_FFS
            int32_t previousj = -1;
#endif
        for (int k = 0; k < ps->nc; k++) 
          {
            int32_t j = ps->cols[k];
            if (LIKELY((j >= jmin) && (j < jmax)))
              {
#ifdef FOR_FFS
                /* For FFS we do not want to count multiplicity here*/
                if (j != previousj)
                    mat->wt[GETJ(mat, j)]++;
                previousj = j;
#else
                mat->wt[GETJ(mat, j)]++;
#endif
              }
          }
      }
#ifdef FOR_FFS
    int32_t j = mat->ncols - 1;
    mat->wt[GETJ(mat, j)] = mat->nrows;
#endif
}

/* initialize Rj[j] for light columns, i.e., for those of weight <= cwmax */
void
fillmat (filter_matrix_t *mat)
{
  int32_t j, *Rj, jmin = mat->jmin, jmax = mat->jmax, wj;

    for(j = jmin; j < jmax; j++){
        wj = mat->wt[GETJ(mat, j)];
	if (wj <= mat->cwmax){
	    Rj = (int32_t *) malloc((wj + 1) * sizeof(int32_t));
	    Rj[0] = 0; /* last index used */
	    mat->R[GETJ(mat, j)] = Rj;
	}
	else{ /* weight is larger than cwmax */
            mat->wt[GETJ(mat, j)] = -mat->wt[GETJ(mat, j)]; // trick!!!
	    mat->R[GETJ(mat, j)] = NULL;
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
filter_matrix_read (filter_matrix_t *mat, purgedfile_stream_ptr ps, int skip)
{
    int lbuf = 100; 
    typerow_t *buf;
#if !defined(USE_MPI)
    int nh = 0;
#endif
    int32_t ibuf, j, k, jmin = mat->jmin, jmax = mat->jmax;
    char *tooheavy = NULL;
    int *wskip, wc;

    mat->nburied = 0;
    buf = (typerow_t *) malloc(lbuf * sizeof(typerow_t));

    fprintf(stderr, "Reading matrix of %d rows and %d columns: excess is %d\n",
            mat->nrows, mat->ncols, mat->nrows - mat->ncols);
    mat->rem_nrows = mat->nrows;
    mat->weight = 0;

    int bmin = mat->nrows, bmax = 0;

    /* heavy columns already have wt < 0 */
    tooheavy = (char *) malloc (mat->ncols * sizeof(char));
    memset (tooheavy, 0, mat->ncols * sizeof(char));

    /* find the skip heaviest columns */
    wskip = (int*) malloc ((skip + 1) * sizeof(int));
    memset (wskip, 0, skip * sizeof(int));
    for(j = 0; j < mat->ncols; j++)
      {
        wc = -mat->wt[j];
        k = skip;
        while ((k > 0) && (wc > wskip[k - 1]))
          wskip[k] = wskip[k - 1], k--;
        wskip[k] = wc;
      }
    for (j = 0; j < mat->ncols; j++)
      {
        wc = -mat->wt[j];
        if ((skip > 0) && (wc >= wskip[skip - 1]))
          {
#if DEBUG >= 1
            fprintf(stderr, "Burying j=%d (wt=%d)\n", j, wc);
#endif
            tooheavy[j] = 1;
            mat->nburied += 1;
            if (wc > bmax) bmax = wc;
            if (wc < bmin) bmin = wc;
          }
      }
    free (wskip);
    fprintf(stderr, "# Number of buried columns is %d", mat->nburied);
    fprintf(stderr, " (min weight=%d, max weigth=%d)\n", bmin, bmax);

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
            if(ps->nc >= lbuf){
                lbuf <<= 1;
                fprintf(stderr, "Warning: doubling lbuf in readmat;");
                fprintf(stderr, " new value is %d\n", lbuf);
                free(buf);
                buf = (typerow_t *)malloc(lbuf * sizeof(typerow_t));
            }
            ibuf = 0;
#ifdef FOR_FFS
            int32_t previousj = -1;
#endif
            for(int k = 0 ; k < ps->nc; k++) {
                int32_t j = ps->cols[k];
                ASSERT_ALWAYS (0 <= j && j < mat->ncols);
#if !defined(USE_MPI)
                if(mat->wt[GETJ(mat, j)] < 0)
                    nb_heavy_j++;
#endif
                // always store j in the right interval...!
#ifdef FOR_FFS
                if((j >= jmin) && (j < jmax) && j != previousj)
#else
                if((j >= jmin) && (j < jmax))
#endif
                {
                    // this will be the weight in the current slice
                    if ((tooheavy == NULL) || (tooheavy[j] == 0)) {
                        mat->weight++;
#ifdef FOR_FFS
                        buf[ibuf++] = (typerow_t) { .id = j, .e = 1};
#else
                        buf[ibuf++] = j;
#endif
                        if(mat->wt[GETJ(mat, j)] > 0){ // redundant test?
                            mat->R[GETJ(mat, j)][0]++;
                            mat->R[GETJ(mat, j)][mat->R[GETJ(mat, j)][0]] = i;
                        }
                    }
#ifdef FOR_FFS
                    previousj = j ;
#endif
                }
#ifdef FOR_FFS
                else if  (previousj == j && ibuf != 0 &&
                                ((tooheavy == NULL) || (tooheavy[j] == 0)))
                    buf[ibuf-1].e++;
#endif
            }
#ifdef FOR_FFS
          int32_t j = mat->ncols-1;
          if((j >= jmin) && (j < jmax))
            {
              if ((tooheavy == NULL) || (tooheavy[j] == 0))
                {
                  mat->weight++;
                  buf[ibuf++] = (typerow_t) { .id = j, .e = 1};
                  if(mat->wt[GETJ(mat, j)] > 0) // redundant test?
                    {
                      mat->R[GETJ(mat, j)][0]++;
                      mat->R[GETJ(mat, j)][mat->R[GETJ(mat, j)][0]] = i;
                    }
                }
            }
#endif

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
            mat->rows[i] = (typerow_t*) malloc ((ibuf+1) * sizeof (typerow_t));
            matLengthRow(mat, i) = ibuf;
            memcpy(mat->rows[i]+1, buf, ibuf * sizeof(typerow_t));
            // sort indices in val to ease row merges
            qsort(mat->rows[i]+1, ibuf, sizeof(typerow_t), cmp);
#if DEBUG >= 1
            /* check all indices are distinct, otherwise this is a bug of
               purge */
            for (int k = 1; k < ibuf; k++)
                if (matCell(mat, i, k) == matCell(mat, i, k+1))
                {
                    fprintf (stderr, "Error, duplicate ideal %x in row %i\n",
                            matCell(mat, i, k), i);
                    exit (1);
                }
#endif
        }
    }
    fprintf(stderr, "read %d rows in %.1f s (%.1f MB/s, %.1f rows/s)\n", 
                    ps->rrows, ps->dt, ps->mb_s, ps->rows_s);
    fprintf (stderr, "\n"); /* to keep last output on screen */
#if !defined(USE_MPI)
    fprintf(stderr, "Number of heavy rows: %d\n", nh);
#endif
    // we need to keep informed of what really happens; this will be an upper
    // bound on the number of active columns, I guess
    mat->rem_ncols = mat->ncols;
    free (buf); 
    if (tooheavy != NULL)
        free (tooheavy);
    return 1;
}

void
print_row(filter_matrix_t *mat, int i)
{
    fprintRow (stdout, mat->rows[i]);
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
    free(mat->R[GETJ(mat, j)]);
    mat->R[GETJ(mat, j)] = NULL;
}

// Don't touch to R[j][0]!!!!
void
remove_i_from_Rj(filter_matrix_t *mat, int i, int j)
{
    // be dumb for a while
    int k;

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
    for(k = 1; k <= matLengthRow(mat, i); k++)
	if(matCell(mat, i, k) == j)
	    break;
    ASSERT(k <= matLengthRow(mat, i));
    // crunch
    for(++k; k <= matLengthRow(mat, i); k++)
        matCell(mat, i, k-1) = matCell(mat, i, k);
    matLengthRow(mat, i) -= 1;
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
weightSum(filter_matrix_t *mat, int i1, int i2, MAYBE_UNUSED int32_t j)
{
    int k1, k2, w, len1, len2;

    len1 = (isRowNull(mat, i1) ? 0 : matLengthRow(mat, i1));
    len2 = (isRowNull(mat, i2) ? 0 : matLengthRow(mat, i2));
    if((len1 == 0) || (len2 == 0))
        fprintf(stderr, "i1=%d i2=%d len1=%d len2=%d\n", i1, i2, len1, len2);
#ifdef FOR_FFS /* look for the exponents of j in i1 i2*/
    int e1 = 0, e2 = 0;
    int d;
    int l;
    for (l = 1 ; l <= matLengthRow(mat, i1) ; l++)
        if (matCell(mat, i1, l) == j)
            e1 = mat->rows[i1][l].e;
    for (l = 1 ; l <= matLengthRow(mat, i2) ; l++)
        if (matCell(mat, i2, l) == j)
            e2 = mat->rows[i2][l].e;

    ASSERT (e1 != 0 && e2 != 0);

    d  = (int) gcd_int64 ((int64_t) e1, (int64_t) e2);
    e1 /= -d;
    e2 /= d;
#endif
    k1 = k2 = 1;
    w = 0;
    while((k1 <= len1) && (k2 <= len2))
    {
        if(matCell(mat, i1, k1) < matCell(mat, i2, k2))
        {
#ifdef FOR_FFS
            w += weight_ffs (e2 * mat->rows[i1][k1].e);
#else
            w++;
#endif
            k1++;
        }
        else if(matCell(mat, i1, k1) > matCell(mat, i2, k2))
        {
#ifdef FOR_FFS
            w += weight_ffs (e1 * mat->rows[i2][k2].e);
#else
            w++;
#endif
            k2++;
        }
        else
        {
#ifdef FOR_FFS
            w += weight_ffs (e2*mat->rows[i1][k1].e + e1*mat->rows[i2][k2].e);
#endif
            k1++; 
            k2++;
        }
    }
    // finish with k1
    for( ; k1 <= matLengthRow(mat, i1); k1++)
#ifdef FOR_FFS
            w += weight_ffs (e2 * mat->rows[i1][k1].e);
#else
            w++;
#endif
    // finish with k2
    for( ; k2 <= matLengthRow(mat, i2); k2++)
#ifdef FOR_FFS
            w += weight_ffs (e1 * mat->rows[i2][k2].e);
#else
            w++;
#endif
    return w;
}

/* put in ind[0]..ind[m-1] the indices of the m (active) rows containing j,
   and return the total weight of all those rows */
int
fillTabWithRowsForGivenj(int32_t *ind, filter_matrix_t *mat, int32_t j)
{
  int ni = 0, k, i, w = 0;

  for (k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
    if ((i = mat->R[GETJ(mat, j)][k]) != -1)
      {
        ind[ni++] = i;
        w += matLengthRow(mat, i);
        if (ni == mat->wt[GETJ(mat, j)])
          break;
      }
  return w;
}
