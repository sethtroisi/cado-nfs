#include "cado.h"

#include <string.h>
#include <stdio.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "portability.h"
#include "utils_with_io.h"
#include "filter_config.h"
#include "merge_replay_matrix.h"
#include "sparse.h"
#include "markowitz.h"

/***************** memory allocation on R[j] *********************************/

static void
mallocRj (filter_matrix_t *mat, int j, int32_t w)
{
  if (w == 0)
    {
      /* This can happen because we do not renumber columns in purge,
	 thus the range of the j-values is larger than the number of
	 columns. We simply ignore those columns. */
      mat->R[j] = NULL;
      return;
    }
  mat->R[j] = (index_t *) malloc((w + 1) * sizeof(index_t));
  FATAL_ERROR_CHECK(mat->R[j] == NULL, "Cannot allocate memory");
  mat->R[j][0] = 0; /* last index used */
}

static void
reallocRj (filter_matrix_t *mat, int j, int32_t w)
{
  mat->R[j] = (index_t *) realloc (mat->R[j], (w + 1) * sizeof(index_t));
  FATAL_ERROR_CHECK(mat->R[j] == NULL, "Cannot reallocate memory");
  mat->R[j][0] = w;
}

void
freeRj (filter_matrix_t *mat, int j)
{
    free (mat->R[j]);
    mat->R[j] = NULL;
}

/****************************** heap management ******************************/

/* we precompute the weights for w < 256 to avoid the cost of multiple calls
   to comp_weight_function */
static float comp_weight[256];

static inline float
comp_weight_function (int32_t w MAYBE_UNUSED)
{
#define USE_WEIGHT_LAMBDA 0 /* LAMBDA = 0 seems to be better */
#if USE_WEIGHT_LAMBDA == 0
    return 1.0;
#elif USE_WEIGHT_LAMBDA == 1
    return powf (2.0 / 3.0, (float) (w - 2));
#elif USE_WEIGHT_LAMBDA == 2
    return powf (0.5, (float) (w - 2));
#elif USE_WEIGHT_LAMBDA == 3
    return powf (0.8, (float) (w - 2));
#elif USE_WEIGHT_LAMBDA == 4
    return 1.0 / log2f ((float) w);
#elif USE_WEIGHT_LAMBDA == 5
    return 2.0 / (float) w;
#elif USE_WEIGHT_LAMBDA == 6
    return 4.0 / (float) (w * w);
#else
#error "Invalid value of USE_WEIGHT_LAMBDA"
#endif
}

static void
heap_init (heap H, uint32_t nrows)
{
  H->list = malloc (nrows * sizeof (uint32_t));
  H->size = 0;
  H->alloc = nrows;
  H->index = malloc (nrows * sizeof (uint32_t));
  for (unsigned long i = 0; i < nrows; i++)
    H->index[i] = UINT32_MAX;
  for (int32_t w = 0; w < 256; w++)
    comp_weight[w] = comp_weight_function (w);
  printf ("Using weight function lambda=%d for clique removal\n",
          USE_WEIGHT_LAMBDA);
}

static void
heap_clear (heap H)
{
  free (H->list);
  free (H->index);
}

static float
heapRowWeight (uint32_t i, filter_matrix_t *mat)
{
  float W = 0.0;
  int32_t w;
  uint32_t j, k;

  ASSERT(mat->rows[i] != NULL);
  for (k = 1; k <= matLengthRow (mat, i); k++)
    {
      j = matCell (mat, i, k);
      w = mat->wt[j];
      ASSERT(w != 0);
      if (w < 0)
        w = -w;
      if (w > 255)
        w = 255; /* saturate to 255 */
      W += comp_weight[w];
    }
  return W;
}

#if 0
static void
check_heap (heap H, filter_matrix_t *mat)
{
  for (unsigned int i = 0; i < H->size; i++)
    {
      uint32_t j = H->list[i];
      ASSERT_ALWAYS(H->index[j] == i);
      ASSERT_ALWAYS(mat->rows[j] != NULL);
    }
}
#endif

/* move down entry n of heap */
static void
moveDown (heap H, filter_matrix_t *mat, uint32_t n)
{
  uint32_t i = H->list[n], j;
  float w;

  w = heapRowWeight (i, mat);
  while (2 * n + 1 < H->size)
    {
      uint32_t left = 2 * n + 1, right = 2 * n + 2, son;
      if (right >= H->size ||
          heapRowWeight (H->list[left], mat) > heapRowWeight (H->list[right], mat))
        son = left; /* compare with left son */
      else
        son = right; /* compare with right son */
      if (heapRowWeight (H->list[son], mat) > w)
        {
          j = H->list[son];
          H->list[n] = j;
          H->index[j] = n;
          n = son;
        }
      else
        break; /* i has larger weight */
    }
  H->list[n] = i;
  H->index[i] = n;
}

/* remove relation i */
void
heap_delete (heap H, filter_matrix_t *mat, uint32_t i)
{
  uint32_t n = H->index[i];

  ASSERT(H->size > 0);

  /* here, the row of cell n of the heap is invalid */

  /* put last entry in position n (if not already last) */
  H->size --;
  if (n < H->size)
    {
      H->list[n] = H->list[H->size];
      /* adjust heap */
      moveDown (H, mat, n);
    }
}

/* 1,2 -> 0, 3,4 -> 1, 5,6 -> 2, ... */
#define PARENT(i) (((i)+1)/2-1)

/* add relation i */
void
heap_push (heap H, filter_matrix_t *mat, uint32_t i)
{
  uint32_t n, j;
  float w;

  if (H->index[i] == UINT32_MAX)
    {
      ASSERT(H->size < H->alloc);

      w = heapRowWeight (i, mat);
      n = H->size;

      /* move parents down the tree as long as they have smaller weight */
      while (n > 0 && heapRowWeight (H->list[PARENT(n)], mat) < w)
        {
          j = H->list[PARENT(n)];
          H->list[n] = j;
          H->index[j] = n;
          n = PARENT(n);
        }
      /* insert new element */
      H->list[n] = i;
      H->index[i] = n;
      H->size ++;
    }
  else /* relation was already in heap */
    {
      uint32_t j;
      n = H->index[i];
      w = heapRowWeight (i, mat);
      /* if new weight is smaller than old one, move down the heap */
      if ((2 * n + 1 < H->size && w < heapRowWeight (H->list[2 * n + 1], mat)) ||
          (2 * n + 2 < H->size && w < heapRowWeight (H->list[2 * n + 2], mat)))
        moveDown (H, mat, n);
      else /* move up the heap */
        {
          while (n > 0 && w > heapRowWeight (H->list[PARENT(n)], mat))
            {
              j = H->list[PARENT(n)];
              H->list[n] = j;
              H->index[j] = n;
              n = PARENT(n);
            }
          /* insert updated element */
          H->list[n] = i;
          H->index[i] = n;
        }
    }
}

/* return index i of relation with larger weight in heap H */
uint32_t
heap_pop (heap H, filter_matrix_t *mat MAYBE_UNUSED)
{
  ASSERT(H->size > 0);
  return H->list[0];
}

/*****************************************************************************/

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

/* initialize the sparse matrix mat */
void
initMat (filter_matrix_t *mat, int maxlevel, uint32_t keep,
         uint32_t nburied)
{
  mat->keep  = keep;
  mat->mergelevelmax = maxlevel;
  mat->cwmax = maxlevel;
  ASSERT_ALWAYS (mat->cwmax < 255); /* 255 is reserved for saturated values */
  mat->nburied = nburied;

  mat->weight = 0;
  mat->tot_weight = 0;
  mat->rem_ncols = 0;

  mat->rows = (typerow_t **) malloc (mat->nrows * sizeof (typerow_t *));
  ASSERT_ALWAYS (mat->rows != NULL);
  mat->wt = (int32_t *) malloc (mat->ncols * sizeof (int32_t));
  memset (mat->wt, 0, mat->ncols * sizeof (int32_t));
  mat->R = (index_t **) malloc (mat->ncols * sizeof(index_t *));
  ASSERT_ALWAYS(mat->R != NULL);

  /* initialize the heap for heavy rows */
  heap_init (mat->Heavy, mat->nrows);
}

void
clearMat (filter_matrix_t *mat)
{
  uint64_t i, j;

  for (i = 0; i < mat->nrows; i++)
    free (mat->rows[i]);
  free (mat->rows);
  free (mat->wt);
  for (j = 0; j < mat->ncols; j++)
    freeRj (mat, j);
  free (mat->R);

  /* free the memory for the heap of heavy rows */
  heap_clear (mat->Heavy);
}

#ifndef FOR_DL
/* Renumber the non-zero columns to contiguous values [0, 1, 2, ...]
 * We can use this function only for factorization because in this case we do
 * not need the indexes of the columns (contrary to DL where the indexes of the
 * column are printed in the history file). */
static void
renumber_columns (filter_matrix_t *mat)
{
  index_t *p = NULL;

  p = (index_t *) malloc (mat->ncols * sizeof (index_t));
  ASSERT_ALWAYS (p != NULL);

  /* first compute the mapping of column indices */
  index_t h = 0;
  for (uint64_t j = 0; j < mat->ncols; j++)
  {
    /* at this point wt[j] = 0 if ideal j never appears in relations,
     * or wt[j] > 0 if ideal j appears in relations */
    if (mat->wt[j] != 0)
    {
      /* index j is mapped to h, with h <= j */
      p[j] = h;
#if TRACE_COL >= 0
      if (TRACE_COL == h || TRACE_COL == j)
        printf ("TRACE_COL: column %lu is renumbered into %u\n", j, h);
#endif
      mat->wt[h] = mat->wt[j];
      h++;
    }
  }
  /* h should be equal to rem_ncols, which is the number of columns with
   * non-zero weight */
  ASSERT_ALWAYS(h == mat->rem_ncols);

  /* Realloc mat->wt */
  mat->wt = realloc (mat->wt, h * sizeof (int32_t));
  /* Reset mat->ncols to behave as if they were no non-zero columns. */
  mat->ncols = h;

  /* apply mapping to the rows. As p is a non decreasing function, the rows are
   * still sorted after this operation. */
  for (uint64_t i = 0; i < mat->nrows; i++)
    for (index_t j = 1; j <= mat->rows[i][0]; j++)
      mat->rows[i][j] = p[mat->rows[i][j]];

  free (p);
}
#endif

/* initialize Rj[j] for light columns, i.e., for those of weight <= cwmax */
static void
InitMatR (filter_matrix_t *mat)
{
  index_t h;
  int32_t w;
  int32_t wmax = mat->cwmax;

#ifndef FOR_DL
  renumber_columns (mat);
#endif

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (h = 0; h < mat->ncols; h++)
    {
      w = mat->wt[h];
      if (w <= wmax)
        mallocRj (mat, h, w);
      else /* weight is larger than cwmax */
        {
          mat->wt[h] = -mat->wt[h]; // trick!!!
          mat->R[h] = NULL;
        }
    }
}

static void
recompute_weights (filter_matrix_t *mat)
{
  uint64_t i, h;

  /* reset the weights to zero */
  memset (mat->wt, 0, mat->ncols * sizeof (int32_t));

  /* recompute the weights */
  for (i = 0; i < mat->nrows; i++)
    if (mat->rows[i] != NULL) /* row is still active */
      for (unsigned int k = 1 ; k <= matLengthRow(mat, i); k++)
        {
          h = matCell(mat, i, k);
          mat->wt[h] ++;
#if TRACE_COL >= 0
          if (h == TRACE_COL)
            printf ("TRACE_COL: found ideal %lu in row %lu\n", h, i);
#endif
        }
}

static void
reinitMatR (filter_matrix_t *mat)
{
  index_t h;
  int32_t w;
  int32_t wmax = mat->cwmax;

  for (h = 0; h < mat->ncols; h++)
    {
      w = mat->wt[h];
      if (w <= wmax)
        {
          reallocRj (mat, h, w);
          mat->R[h][0] = 0; /* we reset the weights to 0 */
        }
      else /* weight is larger than cwmax */
        {
          mat->wt[h] = -mat->wt[h]; // trick!!!
          mat->R[h] = NULL;
        }
    }
}

static void
fillR (filter_matrix_t *mat)
{
  uint64_t i, h;

  for (i = 0; i < mat->nrows; i++)
    if (mat->rows[i] != NULL) /* row i is still active */
      for (unsigned int k = 1 ; k <= matLengthRow(mat, i); k++)
        {
          h = matCell(mat, i, k);
          if (mat->wt[h] > 0)
            {
              ASSERT_ALWAYS(mat->R[h] != NULL);
              mat->R[h][0]++;
              mat->R[h][mat->R[h][0]] = i;
            }
        }

#if TRACE_COL >= 0
  h = TRACE_COL;
  printf ("TRACE_COL: weight of ideal %lu is %d\n", h, mat->wt[h]);
  ASSERT_ALWAYS(mat->wt[h] <= 0 || (uint32_t) mat->wt[h] == mat->R[h][0]);
#endif
}

void
recomputeR (filter_matrix_t *mat)
{
  /* recompute the column weights */
  recompute_weights (mat);

  /* re-allocate the matrix R */
  reinitMatR (mat);

  /* fill the matrix R */
  fillR (mat);

  /* recompute the Markowitz structure */
  MkzClear (mat, 0);
  MkzInit (mat, 0);
}

void
matR_disable_cols (filter_matrix_t *mat, const char *infilename)
{
  FILE *file = NULL;
  char buf[256];

  file = fopen_maybe_compressed (infilename, "r");
  ASSERT_ALWAYS (file != NULL);

  int stop = 0;
  while (!stop)
  {
    if (fgets(buf, 256, file) == NULL)
      stop = 1;
    else if (buf[0] != '#')
    {
      size_t n = strnlen(buf, 256);
      ASSERT_ALWAYS(n != 256);

      index_t h;
      int ret = sscanf (buf, "%" SCNid "\n", &h);
      ASSERT_ALWAYS (ret == 1);
      if (h < mat->ncols)
      {
        if (mat->R[h] != NULL)
          freeRj (mat, h);
        if (mat->wt[h] > 0)
          mat->wt[h] = -mat->wt[h]; // trick!!!
      }
    }
  }

  fclose_maybe_compressed (file, infilename);
}

/* callback function called by filter_rels */
void * insert_rel_into_table (void *context_data, earlyparsed_relation_ptr rel)
{
  filter_matrix_t *mat = (filter_matrix_t *) context_data;

  /* XXX For now can't use my_malloc, because rows are realloc later */
  mat->rows[rel->num] = (typerow_t*) malloc ((rel->nb + 1) * sizeof (typerow_t));
  FATAL_ERROR_CHECK(mat->rows[rel->num] == NULL, "Cannot allocate memory");
  matLengthRow(mat, rel->num) = rel->nb;
  mat->tot_weight += rel->nb;

  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    exponent_t e = rel->primes[i].e;
    /* For factorization, they should not be any multiplicity here.
       For DL we do not want to count multiplicity in mat->wt */
#ifndef FOR_DL
    ASSERT_ALWAYS (e == 1);
#endif
    if (mat->wt[h] == 0)
    {
      mat->wt[h] = 1;
      mat->rem_ncols++;
    }
    else if (mat->wt[h] != SMAX(int32_t))
      mat->wt[h]++;

    setCell(mat->rows[rel->num][i+1], h, e);
  }

  /* sort indices to ease row merges */
  qsort(&(mat->rows[rel->num][1]), rel->nb, sizeof(typerow_t), cmp_typerow_t);

  return NULL;
}


void
filter_matrix_read (filter_matrix_t *mat, const char *purgedname)
{
  uint64_t nread, i, h;
  char *fic[2] = {(char *) purgedname, NULL};

  /* read all rels */
  nread = filter_rels(fic, (filter_rels_callback_t) &insert_rel_into_table, mat,
                      EARLYPARSE_NEED_INDEX, NULL, NULL);
  ASSERT_ALWAYS(nread == mat->nrows);
  mat->rem_nrows = nread;

  /* print weight count */
  uint64_t *nbm;
  nbm = (uint64_t *) malloc ((mat->mergelevelmax + 1) * sizeof (uint64_t));
  memset (nbm, 0, (mat->mergelevelmax + 1) * sizeof (uint64_t));
  for (h = 0; h < mat->ncols; h++)
    if (mat->wt[h] <= mat->mergelevelmax)
      nbm[mat->wt[h]]++;
  for (h = 1; h <= (uint64_t) mat->mergelevelmax; h++)
    printf ("There are %" PRIu64 " column(s) of weight %" PRIu64 "\n", nbm[h], h);
  ASSERT_ALWAYS(mat->rem_ncols == mat->ncols - nbm[0]);
  free (nbm);

  /* Bury heavy coloumns. The 'nburied' heaviest column are buried. */
  /* Buried columns are not took into account by merge. */
  uint64_t weight_buried = 0;
  if (mat->nburied)
  {
    uint64_t *heaviest = NULL;
    heaviest = (uint64_t *) malloc (mat->nburied * sizeof (uint64_t));
    for (i = 0; i < mat->nburied; i++)
      heaviest[i] = UMAX(uint64_t);

    /* Find the 'nburied' heaviest columns */
    for (i = 0; i < mat->ncols; i++)
    {
      int32_t w_cur = mat->wt[i];
      uint64_t j = mat->nburied;

      while (j > 0 && (heaviest[j-1] == UMAX(uint64_t) ||
                       w_cur > mat->wt[heaviest[j-1]]))
      {
        if (j < mat->nburied)
          heaviest[j] = heaviest[j-1];
        j--;
      }
      if (j < mat->nburied)
        heaviest[j] = i;
    }

    int32_t buried_max = mat->wt[heaviest[0]];
    int32_t buried_min = mat->wt[heaviest[mat->nburied-1]];

    /* Compute weight of buried part of the matrix. */
    for (i = 0; i < mat->nburied; i++)
    {
      int32_t w = mat->wt[heaviest[i]];
      weight_buried += w;
#if DEBUG >= 1
      fprintf(stderr, "# Burying j=%" PRIu64 " (wt = %" PRId32 ")\n",
                      heaviest[i], w);
#endif
    }
    printf("# Number of buried columns is %" PRIu64 " (min_weight=%" PRId32 ", "
           "max_weight=%" PRId32 ")\n", mat->nburied, buried_min, buried_max);
    printf("# Weight of the buried part of the matrix: %" PRIu64 "\n",
           weight_buried);

    /* Remove buried columns from rows in mat structure */
    printf("# Start to remove buried columns from rels...\n");
    fflush (stdout);
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < mat->nrows; i++)
    {
      unsigned int k = 1;
      for (unsigned int j = 1; j <= matLengthRow(mat, i) ; j++)
      {
        index_t cur_id = matCell(mat, i, j);
        int32_t w = mat->wt[cur_id];
        unsigned int is_buried;
        if (w < buried_min)
          is_buried = 0;
        else if (w > buried_min)
          is_buried = 1;
        else /* w == buried_min */
        {
          unsigned int l = mat->nburied;
          is_buried = 0;
          while (l > 0 && w == mat->wt[heaviest[l-1]])
          {
            if (heaviest[l-1] == cur_id)
              is_buried = 1;
            l--;
          }
        }
        if (!is_buried) /* not a buried column */
          mat->rows[i][k++] = mat->rows[i][j];
      }
      matLengthRow(mat, i) = k-1;
      mat->rows[i] = (typerow_t*) realloc (mat->rows[i], k * sizeof (typerow_t));
    }
    printf ("# Done\n");

    /* compute the matrix weight */
    mat->weight = 0;
    for (i = 0; i < mat->nrows; i++)
      mat->weight += matLengthRow(mat, i);

    free (heaviest);
  }
  else
  {
    printf("# No columns were buried.\n");
    mat->weight = mat->tot_weight;
  }
  printf("# Weight of the active part of the matrix: %" PRIu64 "\n# Total "
         "weight of the matrix: %" PRIu64 "\n", mat->weight, mat->tot_weight);
  ASSERT_ALWAYS (mat->weight + weight_buried == mat->tot_weight);

  /* Allocate mat->R[h] if necessary */
  InitMatR (mat);

  /* Re-read all rels (in memory) to fill-in mat->R */
  printf("# Start to fill-in columns of the matrix...\n");
  fflush (stdout);
  fillR (mat);
  printf ("# Done\n");
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
remove_i_from_Rj(filter_matrix_t *mat, index_t i, int j)
{
  unsigned int k, n = mat->R[j][0];

  /* R[j] should not be empty */
  ASSERT(mat->R[j] != NULL);

  for (k = 1; k <= n; k++)
    if (mat->R[j][k] == i)
      {
        mat->R[j][k] = mat->R[j][n];
        mat->R[j][0] = n - 1;
        return;
      }
  ASSERT_ALWAYS(0);
}

// cell M[i, j] is incorporated in the data structure. It is used
// later on in cases where i is not already present in row[j].
void
add_i_to_Rj(filter_matrix_t *mat, int i, int j)
{
  int l = mat->R[j][0] + 1;

  reallocRj (mat, j, l);
  mat->R[j][l] = i;
}

/* return the weight of the relation obtained when adding relations i1 and i2
*/
int
weightSum(filter_matrix_t *mat, int i1, int i2, MAYBE_UNUSED int32_t j)
{
  unsigned int k1, k2, w, len1, len2;

    len1 = (isRowNull(mat, i1) ? 0 : matLengthRow(mat, i1));
    len2 = (isRowNull(mat, i2) ? 0 : matLengthRow(mat, i2));
#if DEBUG >= 1
    if((len1 == 0) || (len2 == 0))
        fprintf(stderr, "i1=%d i2=%d len1=%d len2=%d\n", i1, i2, len1, len2);
#endif
#ifdef FOR_DL /* look for the exponents of j in i1 and i2*/
    int e1 = 0, e2 = 0;
    int d;
    unsigned int l;
    for (l = 1 ; l <= len1 ; l++)
        if ((int) matCell(mat, i1, l) == j)
            e1 = mat->rows[i1][l].e;
    for (l = 1 ; l <= len2 ; l++)
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
            w += (e2 * mat->rows[i1][k1].e == 0) ? 0 : 1;
#else
            w++;
#endif
            k1++;
        }
        else if(matCell(mat, i1, k1) > matCell(mat, i2, k2))
        {
#ifdef FOR_DL
            w += (e1 * mat->rows[i2][k2].e == 0) ? 0 : 1;
#else
            w++;
#endif
            k2++;
        }
        else
        {
#ifdef FOR_DL
            w += (e2*mat->rows[i1][k1].e + e1*mat->rows[i2][k2].e == 0) ? 0 : 1;
#endif
            k1++;
            k2++;
        }
    }
#ifdef FOR_DL
    // finish with k1
    for( ; k1 <= matLengthRow(mat, i1); k1++)
      w += (e2 * mat->rows[i1][k1].e == 0) ? 0 : 1;
    // finish with k2
    for( ; k2 <= matLengthRow(mat, i2); k2++)
      w += (e1 * mat->rows[i2][k2].e == 0) ? 0 : 1;
#else
    w += matLengthRow(mat, i1) + 1 - k1;
    w += matLengthRow(mat, i2) + 1 - k2;
#endif
    return w;
}

/* put in ind[0]..ind[m-1] the indices of the m (active) rows containing j */
void
fillTabWithRowsForGivenj(int32_t *ind, filter_matrix_t *mat, int32_t j)
{
  int ni = 0;
  unsigned int k;

  for (k = 1; k <= mat->R[j][0]; k++)
    ind[ni++] = mat->R[j][k];
}
