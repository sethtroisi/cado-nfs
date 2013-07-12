#include "cado.h"

#include <string.h>
#include <stdio.h>

#include "portability.h"
#include "utils.h"
#include "filter_utils.h"
#include "merge_opts.h"
#include "filter_matrix.h"
#include "sparse.h"

typedef struct {
  index_t i;
  index_t w;
} perm_t;

//tmp
ideal_merge_t **rel_compact;
weight_t *ideals_weight;


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
initMat (filter_matrix_t *mat)
{
    index_t **R;

    mat->rows = (typerow_t **) malloc (mat->nrows * sizeof (typerow_t *));
    ASSERT_ALWAYS (mat->rows != NULL);
    mat->wt = (int*) malloc (mat->ncols * sizeof (int));
    memset (mat->wt, 0, mat->ncols * sizeof (int));
    R = (index_t **) malloc (mat->ncols * sizeof(index_t *));
    ASSERT_ALWAYS(R != NULL);
    mat->R = R;

    rel_compact = (ideal_merge_t **) malloc (mat->nrows*sizeof(ideal_merge_t*));
    ideals_weight = (weight_t *) malloc (mat->ncols * sizeof (weight_t));
    MEMSETZERO(ideals_weight, mat->ncols);
}

void
clearMat (filter_matrix_t *mat)
{
  index_t i, j;

  for (i = 0; i < mat->nrows; i++)
    free (mat->rows[i]);
  free (mat->rows);
  free (mat->wt);
  for (j = 0; j < mat->ncols - 0; j++)
    free (mat->R[j]);
  free (mat->R);

  free(rel_compact);
  free(ideals_weight);
}

void
checkData(filter_matrix_t *mat)
{
  index_t j, nbj = 0;

    for(j = 0; j < mat->ncols; j++)
	if(mat->wt[j])
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
    int32_t jmin = 0, jmax = mat->ncols;
    for (; purgedfile_stream_get(ps,NULL) >= 0;) 
      {
#ifdef FOR_DL
            int32_t previousj = -1;
#endif
        for (int k = 0; k < ps->nc; k++) 
          {
            int32_t j = ps->cols[k];
            if (LIKELY((j >= jmin) && (j < jmax)))
              {
#ifdef FOR_DL
                /* For DL we do not want to count multiplicity here*/
                if (j != previousj)
                    mat->wt[j]++;
                previousj = j;
#else
                mat->wt[j]++;
#endif
              }
          }
      }
}

/* initialize Rj[j] for light columns, i.e., for those of weight <= cwmax */
void
fillmat (filter_matrix_t *mat)
{
  index_t j, jmin = 0, jmax = mat->ncols;
  int wj;
  index_t *Rj;

  for(j = jmin; j < jmax; j++)
  {
    wj = mat->wt[j];
    if (wj <= mat->cwmax)
    {
      Rj = (index_t *) malloc((wj + 1) * sizeof(index_t));
      Rj[0] = 0; /* last index used */
      mat->R[j] = Rj;
    }
    else /* weight is larger than cwmax */
    {
      mat->wt[j] = -mat->wt[j]; // trick!!!
      mat->R[j] = NULL;
    }
  }
}

void *
thread_insert (buf_arg_t *arg)
{
  unsigned int j;
  unsigned long cpy_cpt_rel_b;
  buf_rel_t *my_rel;

  cpy_cpt_rel_b = cpt_rel_b;
  for ( ; ; )
  {
    while (cpt_rel_a == cpy_cpt_rel_b)
    {
      if (!is_finish())
        NANOSLEEP();
      else if (cpt_rel_a == cpy_cpt_rel_b)
        pthread_exit(NULL);
    }

    j = (unsigned int) (cpy_cpt_rel_b & (SIZE_BUF_REL - 1));
    my_rel = &(arg->rels[j]);

    if (cpt_rel_a == cpy_cpt_rel_b + 1)
      NANOSLEEP();

    //FIXME big bug this does not take into account the exponent....
    arg->info.nprimes += insert_rel_in_table_with_e (my_rel, 0, 0, rel_compact,
                                                      ideals_weight);
    arg->info.W += (double) my_rel->nb;

    test_and_print_progress_now ();
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}

/* Callback function called by prempt_scan_relations */
/* Reads a matrix file.

   mat->wt is already filled:
   * if mat->wt[j] > 0, it is the weight of column j that we consider
   * if mat->wt[j] < 0, we don't consider column j

   FIXME mat->weight is correct on exit.
*/

void
filter_matrix_read (filter_matrix_t *mat, const char *purgedname)
{
    index_t j;
    info_mat_t info;

    fprintf(stderr, "Reading matrix of %d rows and %d columns: excess is %d\n",
            mat->rem_nrows, mat->rem_ncols, mat->rem_nrows - mat->rem_ncols);
    mat->weight = 0;

    char *fic[2] = {(char *) purgedname, NULL};
    info = process_rels (fic, &thread_insert, NULL, 0, NULL, NULL, STEP_MERGE);

    ASSERT_ALWAYS (info.nprimes == (index_t) mat->rem_ncols);
    ASSERT_ALWAYS (info.nrels == (index_t) mat->rem_nrows);

    if (mat->nburied)
    {
      // Buried the 'nburied'-first cols (should be the heaviest by construction)
      /* heavy columns already have wt < 0 */
      int bmin = mat->nrows, bmax = 0;

      for (j = 0; j < mat->nburied; j++)
      {
        int wc = -mat->wt[j];
        // not sure an assert is relevant; maybe just a warning
        //ASSERT_ALWAYS (wc > 0);
        if (wc > bmax)
          bmax = wc;
        if (wc < bmin)
          bmin = wc;
#if DEBUG >= 1
        fprintf(stderr, "Burying j=%d (wt=%d)\n", j, wc);
#endif
      }
      fprintf(stderr, "# Number of buried columns is %"PRid"", mat->nburied);
      fprintf(stderr, " (min weight=%d, max weigth=%d)\n", bmin, bmax);
    }
    else
      fprintf(stderr, "# No columns were buried.\n");

    typerow_t buf[REL_MAX_SIZE];
    index_t i;
    for (i = 0; i < mat->nrows; i++)
    {
      weight_t w = 0;
      weight_t next = 0;
      ideal_merge_t *p;
      for (p = rel_compact[i]; p->id != UMAX(p->id); p++)
        w++;
      ASSERT_ALWAYS (w != 0);

      for(unsigned int k = 0 ; k < w; k++) 
      {
        index_t j = rel_compact[i][k].id;
        ASSERT_ALWAYS (j < mat->ncols);
        if (j >= mat->nburied)
        {
          mat->weight++;
#ifndef FOR_DL
          buf[next] = j;
          ASSERT(rel_compact[i][k].e == 1);
#else
          buf[next].id = rel_compact[i][k].id;
          buf[next].e = rel_compact[i][k].e;
#endif
          next++;
          if(mat->wt[j] > 0)
          {
            mat->R[j][0]++;
            mat->R[j][mat->R[j][0]] = i;
          }
        }
      }
      
      mat->rows[i] = (typerow_t*) malloc ((next+1) * sizeof (typerow_t));
      matLengthRow(mat, i) = next;
      memcpy(mat->rows[i]+1, buf, next * sizeof(typerow_t));
      // sort indices in val to ease row merges
      qsort(mat->rows[i]+1, next, sizeof(typerow_t), cmp);
    }
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
    free(mat->R[j]);
    mat->R[j] = NULL;
}

// Don't touch to R[j][0]!!!!
void
remove_i_from_Rj(filter_matrix_t *mat, int i, int j)
{
  // be dumb for a while
  unsigned int k;

  if(mat->R[j] == NULL)
  {
	  fprintf(stderr, "Row %d already empty\n", j);
	  return;
  }
  for(k = 1; k <= mat->R[j][0]; k++)
	  if((int) mat->R[j][k] == i)
    {
#if DEBUG >= 2
      fprintf(stderr, "Removing row %d from R[%d]\n", i, j);
#endif
      mat->R[j][k] = UMAX(index_t);
      break;
	  }
}

// cell M[i, j] is incorporated in the data structure. It is used
// later on in cases where i is not already present in row[j].
void
add_i_to_Rj(filter_matrix_t *mat, int i, int j)
{
  // be semi-dumb for a while
  unsigned int k;
    
#if DEBUG >= 2
  fprintf(stderr, "Adding row %d to R[%d]\n", i, j);
#endif
  for(k = 1; k <= mat->R[j][0]; k++)
    if(mat->R[j][k] == UMAX(index_t))
      break;
  if(k <= mat->R[j][0])
  {
  // we have found a place where it is -1
  mat->R[j][k] = i;
  }
  else
  {
#if DEBUG >= 2
    fprintf(stderr, "WARNING: reallocing things in add_i_to_Rj for R[%d]\n", j);
#endif
    int l = mat->R[j][0]+2;
    mat->R[j] = (index_t *)realloc(mat->R[j], l * sizeof(index_t));
    mat->R[j][l-1] = i;
    mat->R[j][0] = l-1;
  }
}

// Remove j from R[i] and crunch R[i].
void
remove_j_from_row(filter_matrix_t *mat, int i, int j)
{
  unsigned int k;

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
	if((int) matCell(mat, i, k) == j)
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

int
weightSum(filter_matrix_t *mat, int i1, int i2, MAYBE_UNUSED int32_t j)
{
  unsigned int k1, k2, w, len1, len2;

    len1 = (isRowNull(mat, i1) ? 0 : matLengthRow(mat, i1));
    len2 = (isRowNull(mat, i2) ? 0 : matLengthRow(mat, i2));
    if((len1 == 0) || (len2 == 0))
        fprintf(stderr, "i1=%d i2=%d len1=%d len2=%d\n", i1, i2, len1, len2);
#ifdef FOR_DL /* look for the exponents of j in i1 i2*/
    int e1 = 0, e2 = 0;
    int d;
    unsigned int l;
    for (l = 1 ; l <= matLengthRow(mat, i1) ; l++)
        if ((int) matCell(mat, i1, l) == j)
            e1 = mat->rows[i1][l].e;
    for (l = 1 ; l <= matLengthRow(mat, i2) ; l++)
        if ((int) matCell(mat, i2, l) == j)
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
#ifdef FOR_DL
            w += weight_ffs (e2 * mat->rows[i1][k1].e);
#else
            w++;
#endif
            k1++;
        }
        else if(matCell(mat, i1, k1) > matCell(mat, i2, k2))
        {
#ifdef FOR_DL
            w += weight_ffs (e1 * mat->rows[i2][k2].e);
#else
            w++;
#endif
            k2++;
        }
        else
        {
#ifdef FOR_DL
            w += weight_ffs (e2*mat->rows[i1][k1].e + e1*mat->rows[i2][k2].e);
#endif
            k1++; 
            k2++;
        }
    }
    // finish with k1
    for( ; k1 <= matLengthRow(mat, i1); k1++)
#ifdef FOR_DL
            w += weight_ffs (e2 * mat->rows[i1][k1].e);
#else
            w++;
#endif
    // finish with k2
    for( ; k2 <= matLengthRow(mat, i2); k2++)
#ifdef FOR_DL
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
  int ni = 0, w = 0;
  unsigned int k;
  index_t i;

  for (k = 1; k <= mat->R[j][0]; k++)
    if ((i = mat->R[j][k]) != UMAX(index_t))
      {
        ind[ni++] = i;
        w += matLengthRow(mat, i);
        if (ni == mat->wt[j])
          break;
      }
  return w;
}
