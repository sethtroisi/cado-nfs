/* purge --- perform singleton removal and clique removal
 * 
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015
 * Cyril Bouvier, Alain Filbois, Francois Morain, Paul Zimmermann
 * 
 * This file is part of CADO-NFS.
 * 
 * CADO-NFS is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * CADO-NFS is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with CADO-NFS; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
*/

/* References:
 * On the Number Field Sieve Integer Factorisation Algorithm,
 * Stefania Cavallar, PhD Thesis, University of Leiden, 2002.
 *
 * The filtering step of discrete logarithm and integer factorization
 * algorithms.
 * Cyril Bouvier, preprint, 22 pages. 2013.
 */

/*
 * Important remark:
 *   A relation corresponds to a row of the matrix and a ideal corresponds to
 *   a column of the matrix.
 *
 * This program works in two passes over the relation files:
 * - the first pass loads in memory only indexes of columns >= col_min_index
 *   and keeps a count of the weight of each column in cols_weight.
 *   Then, a first step of singleton removal is performed followed by 'nsteps'
 *   steps of singleton removal and clique removal, in order to obtained the
 *   final excess 'keep'.
 * - the second pass goes through the relations again, and dumps the remaining
 *   ones in the format needed by 'merge'.

 * This program uses the following data structures:
 * row_used[i]    - non-zero iff row i is kept (so far)
 * row_compact[i] - list of indexes of columns (greater than or equal to
 *                  col_min_index) of the row i (terminated by a -1 sentinel)
 * cols_weight [h] - weight of the column h in current rows (saturates at 256)
 */

/*
 * Exit value:
 * - 0 if enough relations
 * - 1 if an error occurred (then we should abort the factorization)
 * - 2 if not enough relations
 */

#include "cado.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>  /* for _O_BINARY */
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <pthread.h>
#include <errno.h>
#include <pthread.h>
#ifdef HAVE_LIBGEN_H
#include <libgen.h>
#endif

#include "portability.h"

#include "filter_common.h"
#include "purge_matrix.h"
#include "singleton_removal.h"

//#define USE_CAVALLAR_WEIGHT_FUNCTION

/********************* comp_t struct (clique) ********************************/

/* A clique is a connected component of the graph where the nodes are the rows
   and the edges are the columns of weight 2.
   /!\ It is not a clique is the sense of graph theory.
   We will try to use the "connected component" terminology instead.
*/

typedef struct {
  float w;   /* Weight of the connected component */
  uint64_t i; /* smallest row appearing in the connected component */
} comp_t;

int
comp_cmp_weight (const void *p, const void *q)
{
  float x = ((comp_t *)p)->w;
  float y = ((comp_t *)q)->w;
  return (x <= y ? 1 : -1);
}

/* Contribution of each column of weight w to the weight of the connected
 * component. */
static inline float
comp_weight_function (weight_t w)
{
#ifdef USE_CAVALLAR_WEIGHT_FUNCTION
  if (w >= 3)
    return ldexpf(1, -(w - 1));
  else if (w == 2)
    return 0.25;
  else
    return 0.0;
#else
  if (w >= 3)
    return powf(0.8, (float) (w - 2));
  else if (w == 2)
    return 0.125;
  else
    return 0.0;
#endif
}

/* print info on the weight function that is used */
static inline void
comp_print_info_weight_function ()
{
  fprintf(stdout, "Weight function used during clique removal:\n"
                  "  0     1     2     3     4     5     6     7\n"
                  "0.000 0.000 ");
  for (unsigned int k = 2; k < 8; k++)
    fprintf(stdout, "%0.3f ", comp_weight_function((weight_t) k));
  fprintf(stdout, "\n");
}

/******************** uint64_buffer struct ***********************************/

/* Double buffer:
 * |---done----|---todo----|---free----|
 * | | | | | | | | | | | | | | | | | | |
 *  ^           ^           ^           ^
 *  begin       next_todo   next_free   end
 */
struct uint64_buffer_s {
  uint64_t *begin, *next_todo, *next_free, *end;
};
typedef struct uint64_buffer_s uint64_buffer_t[1];
typedef struct uint64_buffer_s * uint64_buffer_ptr;
typedef const struct uint64_buffer_s * uint64_buffer_srcptr;

#define UINT64_BUFFER_MIN_SIZE 32

/* Init function for uint64_buffer_t */
static inline void
uint64_buffer_init (uint64_buffer_ptr buf, size_t size)
{
  buf->begin = (uint64_t *) malloc_check (sizeof(uint64_t) * size);
  buf->next_todo = buf->begin;
  buf->next_free = buf->begin;
  buf->end = buf->begin + size;
}

/* Reset the buffer and add element to the todo list. */
/* Assume buffer is already initialized with an array of length at least 1. */
static inline void
uint64_buffer_reset_with_one_element (uint64_buffer_ptr buf, uint64_t element)
{
  buf->begin[0] = element;
  buf->next_todo = buf->begin;
  buf->next_free = buf->begin + 1;
}

/* Clear function for uint64_buffer_t */
static inline void
uint64_buffer_clear (uint64_buffer_ptr buf)
{
  free(buf->begin);
}

/* return non-zero if target is in buf
   return 0 otherwise */
static inline int
uint64_buffer_is_in (uint64_buffer_srcptr buf, uint64_t target)
{
  for (uint64_t *p = buf->begin; p < buf->next_free; p++)
    if (UNLIKELY(*p == target))
      return 1;
  return 0;
}

/* This function grows a possible too small uint64_buffer_t. */
static inline void
uint64_buffer_resize (uint64_buffer_ptr buf)
{
  if (UNLIKELY(buf->next_free >= buf->end))
  {
    size_t save_next_todo = buf->next_todo - buf->begin;
    size_t save_next_free = buf->next_free - buf->begin;
    size_t new_size = (buf->end - buf->begin) << 1;
    buf->begin = (uint64_t*) realloc (buf->begin, new_size * sizeof(uint64_t));
    ASSERT_ALWAYS (buf->begin != NULL);
    buf->next_todo = buf->begin + save_next_todo;
    buf->next_free = buf->begin + save_next_free;
    buf->end = buf->begin + new_size;
  }
}

/* push fonction for uint64_buffer_t */
static inline void
uint64_buffer_push_todo (uint64_buffer_ptr buf, uint64_t new_element)
{
  uint64_buffer_resize (buf); /* check that space is enough */
  *(buf->next_free)++ = new_element;
}

/* push fonction for uint64_buffer_t
   return UMAX(uint64_t) if buffer is empty */
static inline uint64_t
uint64_buffer_pop_todo (uint64_buffer_ptr buf)
{
  uint64_t pop_element;
  if (buf->next_todo == buf->next_free)
    pop_element = UMAX(pop_element);
  else
    pop_element = *((buf->next_todo)++);
  return pop_element;
}

/***************** Sorted binary tree of comp_t struct ***********************/

/* Binary tree of comp_t sorted by weight of the connected component.
 * The two sons of a node of index n are the nodes of indexes 2n+1 and 2n+1.
 * The property is: the weight of the comp_t father is less than its 2 sons.
 * so root of the tree is the connected component with the smallest weight.
 * When the tree is full, a new node is inserted in the tree iff its weight is
 * greater than the weight of the root (and the root is removed from the tree).
 * At the end the tree contained the 'size' (= 'alloc') heaviest comp_t.
 */
struct comp_sorted_bin_tree_s {
  size_t alloc;
  size_t size;
  comp_t * tree;
};
typedef struct comp_sorted_bin_tree_s comp_sorted_bin_tree_t[1];
typedef struct comp_sorted_bin_tree_s * comp_sorted_bin_tree_ptr;
typedef const struct comp_sorted_bin_tree_s * comp_sorted_bin_tree_srcptr;

static inline void
comp_sorted_bin_tree_init (comp_sorted_bin_tree_ptr T, size_t max_size)
{
  /* If max_size is even, T->alloc = max_size + 1, in order to avoid a
   * node-father with only one node-son.a We will have one more connected
   * component than needed but we do not care (the computing cost in
   * negligible) */
  if (max_size & 1)
    T->alloc = max_size;
  else
    T->alloc = max_size + 1;

  T->tree = (comp_t *) malloc_check (sizeof(comp_t) * T->alloc);
  T->size = 0;
}

static inline void
comp_sorted_bin_tree_clear (comp_sorted_bin_tree_ptr T)
{
  free(T->tree);
  T->size = T->alloc = 0;
}

static inline void
comp_sorted_bin_tree_qsort (comp_sorted_bin_tree_ptr T,
                            int (*cmp_fct)(const void *, const void *))
{
  qsort (T->tree, T->size, sizeof(comp_t), cmp_fct);
}

static void
comp_sorted_bin_tree_insert (comp_sorted_bin_tree_ptr T, comp_t new_node)
{
  /* If the tree is not full, add new_node at the end and go up. */
  if (UNLIKELY (T->size < T->alloc))
  {
    size_t cur = T->size;
    while (cur)
    {
      size_t father = (cur - 1) >> 1;
      if (UNLIKELY (new_node.w >= T->tree[father].w))
        break;
      T->tree[cur] = T->tree[father];
      cur = father;
    }
    T->size++;
    T->tree[cur] = new_node;
  }
  /* If the tree is full, insert new_node iff the weight of new_node if greater
   * than the weight of the root. In this case, replace the root by new node
   * and go down */
  else if (UNLIKELY(new_node.w > T->tree[0].w))
  {
    size_t cur = 0;
    for (;;) /* In this loop, always, T->tree[cur].w < new_node.w */
    {
      float lightest_weight, other_weight;
      size_t lightest_son = cur * 2 + 1; /* son is son1 */
      /* If there is no son: cur is a leaf. Found! */
      if (UNLIKELY (lightest_son >= T->alloc))
        break;
      lightest_weight = T->tree[lightest_son].w; /* weight of son 1 */
      other_weight = T->tree[lightest_son + 1].w; /* weight of son 2 */
      if (UNLIKELY(other_weight < lightest_weight))
      {
        lightest_son++; /* lightest_son is now son 2 */
        lightest_weight = other_weight;
      }
      /* Now lightest_son and lightest_weight are correct */
      if (UNLIKELY (new_node.w <= lightest_weight)) /* Found! */
        break;
      /* The lightest_son has a weight lighter than new_node, so we have to go
       * down the tree, but before the son replaces its father. */
      T->tree[cur] = T->tree[lightest_son];
      cur = lightest_son;
    }
    T->tree[cur] = new_node;
  }
  /* else do nothing */
}

/************* Code for clique (= connected component) removal ***************/

/* Compute connected component beginning at row clique->i
 * The weight of the connected component is written in clique->w
 *
 * Multithread version.
 *
 * Return the number of rows in the connected component or 0 if it exists a
 * row i1 with i1 < clique->i in the connected component (it implies that this
 * component has already been found, by this thread or another, when we have
 * treated the row i1.
 */
static uint64_t
compute_connected_component_mt (comp_t *clique, purge_matrix_srcptr mat,
                                uint64_buffer_ptr row_buffer)
{
  uint64_t cur_row;
  index_t *h;

  /* Reset buffer and add clique->i to the list of row to explore */
  uint64_buffer_reset_with_one_element (row_buffer, clique->i);
  clique->w = 0.; /* Set initial weight of the connected component to 0. */

  /* Loop on all connected rows */
  while ((cur_row = uint64_buffer_pop_todo (row_buffer)) != UMAX(cur_row))
  {
    /* Loop on all columns of the current row */
    for (h = mat->row_compact[cur_row]; *h != UMAX(*h); h++)
    {
      index_t cur_h = *h;
      weight_t cur_h_weight = mat->cols_weight[cur_h];
      clique->w += comp_weight_function (cur_h_weight);
      if (UNLIKELY(cur_h_weight == 2))
      {
        uint64_t the_other_row = mat->sum2_row[cur_h] - cur_row;
        /* First, if the_other_row < clique.i, the connected component was
         * already found (by this thread or another). return 0 */
        if (the_other_row < clique->i)
          return 0;
        /* If the_other_row is not already in the buffer, add it as a todo. */
        if (!uint64_buffer_is_in (row_buffer, the_other_row))
          uint64_buffer_push_todo (row_buffer, the_other_row);
      }
    }
  }

  /* Return the nb of rows in the connected component */
  return (row_buffer->next_todo - row_buffer->begin);
}

/* Delete connected component of row current_row
 *
 * WARNING: this code itself is multithread compatible, but
 * it calls delete_row, which is NOT compatible!
 */
static void
delete_connected_component (purge_matrix_ptr mat, uint64_t cur_row,
                            uint64_buffer_ptr row_buffer)
//, uint64_t *ncols,
{
  index_t *h;

  /* Reset buffer and add clique->i to the list of row to explore */
  uint64_buffer_reset_with_one_element (row_buffer, cur_row);

  /* Loop on all connected rows */
  while ((cur_row = uint64_buffer_pop_todo (row_buffer)) != UMAX(cur_row))
  {
    /* Loop on all columns of the current row */
    for (h = mat->row_compact[cur_row]; *h != UMAX(*h); h++)
    {
      index_t cur_h = *h;
      weight_t cur_h_weight = mat->cols_weight[cur_h];
      /* We might have some H->hashcount[h] = 3, which is decreased to 2, but
       * we don't want to treat that case. Thus we check in addition that
       * mat->sum2_row[h] <> 0, which only occurs when H->hashcount[h] = 2
       * initially.*/
      if (UNLIKELY(cur_h_weight == 2 && mat->sum2_row[cur_h]))
      {
        uint64_t the_other_row = mat->sum2_row[cur_h] - cur_row;
        /* If the_other_row is not already in the buffer, add it as a todo. */
        if (!uint64_buffer_is_in (row_buffer, the_other_row))
          uint64_buffer_push_todo (row_buffer, the_other_row);
      }
    }
  }

  /* Now, we deleted all rows explored */
  for (uint64_t *pt = row_buffer->begin; pt < row_buffer->next_todo; pt++)
    purge_matrix_delete_row (mat, *pt);
    /* mat->nrows and mat->ncols are updated by purge_matrix_delete_row. */
}

/***** Functions to compute mat->sum2_row array (mono and multi thread) *******/

#ifdef HAVE_SYNC_FETCH /* Multithread code */
/* The main structure for the working pthreads pool which compute mat->sum2_row
 */
typedef struct sum2_mt_data_s {
  // Read only part
  unsigned int pthread_number;
  unsigned int nthreads;
  // Read-write part
  purge_matrix_ptr mat;
  pthread_t pthread;
} sum2_mt_data_t;

static void *
compute_sum2_row_mt (void *pt)
{
  sum2_mt_data_t *data = (sum2_mt_data_t *) pt;
  purge_matrix_ptr mat = data->mat;
  uint64_t j = mat->nrows_init / data->nthreads;
  uint64_t i = (j * data->pthread_number) & ~((size_t) (BV_BITS - 1));
  uint64_t end = (data->pthread_number == data->nthreads - 1) ? mat->nrows_init
                : (j * (data->pthread_number + 1)) & ~((size_t) (BV_BITS - 1));
  bv_t bv, *pbv;
  index_t h, *myrowcompact;

  pbv = mat->row_used->p + (i >> LN2_BV_BITS);
  for (; i < end; i += BV_BITS)
  {
    for (j = i, bv = *pbv++; bv; ++j, bv >>= 1)
    {
      if (LIKELY(bv & 1))
      {
        myrowcompact = mat->row_compact[j];
        for (;;)
        {
          h = *myrowcompact++;
          if (UNLIKELY (h == UMAX(h)))
            break;
          if (LIKELY (mat->cols_weight[h] == 2))
            __sync_add_and_fetch (mat->sum2_row + h, j);
        }
      }
    }
  }
  pthread_exit (NULL);
  return NULL;
}
#endif

static inline void
compute_sum2_row (purge_matrix_ptr mat, unsigned int nthreads)
{
  memset(mat->sum2_row, 0, mat->col_max_index * sizeof(uint64_t));
#ifndef HAVE_SYNC_FETCH /* monothread */
  if (nthreads > 1)
    fprintf (stdout, "# INFO: Cannot use multithread code for compute_sum2_row:"
                     " HAVE_SYNC_FETCH is not defined\n");
  for (uint64_t i = 0; i < mat->nrows_init; i++)
  {
    if (bit_vector_getbit(mat->row_used, (size_t) i))
    {
      index_t h, *myrowcompact;
      for (myrowcompact = mat->row_compact[i]; (h = *myrowcompact++) != UMAX(h);)
        if (mat->cols_weight[h] == 2)
          mat->sum2_row[h] += i;
    }
  }
#else
  sum2_mt_data_t *th_data = malloc_check (nthreads * sizeof (sum2_mt_data_t));
  for (size_t i = nthreads; i--; )
  {
    th_data[i].pthread_number = i;
    th_data[i].nthreads = nthreads;
    th_data[i].mat = mat;
    if (pthread_create (&(th_data[i].pthread), NULL, compute_sum2_row_mt,
                                                     (void *) (th_data + i)))
    {
      perror ("compute_sum2_row pthread creation failed\n");
      exit (1);
    }
  }
  for (size_t i = nthreads; i--; )
    pthread_join (th_data[i].pthread, NULL);
  free (th_data);
#endif
}

/*************** Multithread code for clique removal *************************/

/* The main structure for the working pthreads pool which compute the
   heaviest cliques */
typedef struct comp_mt_thread_data_s {
  unsigned int th_id; /* read only */
  unsigned int nthreads; /* read only */
  pthread_t pthread;
  purge_matrix_srcptr mat; /* Read only */
  comp_sorted_bin_tree_t comp_tree; /* sorted tree computed by the thread */
} comp_mt_thread_data_t;
/* Size of the block of rows treated by a thread during the computation of the
 * connected component (has to be of the form: BV_BITS << x)*/
#define COMP_MT_ROWS_BLOCK (BV_BITS<<4)

/* This MT function fill the tree pth->comp_tree with the pth->comp_tree->alloc
 * heaviest connected components whose smallest row belongs in
 *   [ (data->th_id + j*nthreads) * COMP_MT_ROWS_BLOCK,
                        (data->th_id + j*nthreads + 1) * COMP_MT_ROWS_BLOCK],
 * with j=0,1,2.. until the interval does not intersect [0..mat->nrows_init-1].
 */
static void *
compute_sorted_list_of_connected_components_mt (void *pt)
{
  comp_mt_thread_data_t *data = (comp_mt_thread_data_t *) pt;
  uint64_t end_step_rows;
  bv_t bv, *pbv;
  uint64_buffer_t buf;
  comp_t clique;

  /* Init of the structure & malloc. */
  uint64_buffer_init (buf, UINT64_BUFFER_MIN_SIZE);

  // Now the first begin of the search
  clique.i = data->th_id * COMP_MT_ROWS_BLOCK;
  if (UNLIKELY (clique.i >= data->mat->nrows_init)) pthread_exit (NULL);
  // And the first end
  end_step_rows = clique.i + COMP_MT_ROWS_BLOCK;

  pbv = data->mat->row_used->p + (clique.i >> LN2_BV_BITS);
  while (clique.i < data->mat->nrows_init)
  {
    uint64_t save_i = clique.i;
    for (bv = *pbv++; bv; ++(clique.i), bv >>= 1)
    {
      if (LIKELY(bv & 1))
      {
        unsigned int nb_rows = compute_connected_component_mt (&(clique),
                                              data->mat, buf);
        if (UNLIKELY (!nb_rows))
          continue;
        comp_sorted_bin_tree_insert (data->comp_tree, clique);
      }
    }
    clique.i = save_i + BV_BITS;
    // Have we reach the actuel end ?
    if (UNLIKELY (clique.i == end_step_rows))
    {
      // The new begin & end.
      end_step_rows += data->nthreads * COMP_MT_ROWS_BLOCK;
      clique.i = end_step_rows - COMP_MT_ROWS_BLOCK;
      pbv = data->mat->row_used->p + (clique.i >> LN2_BV_BITS);
    }
  }
  uint64_buffer_clear (buf);
  /* Re-order the connected component by decreasing weight. */
  comp_sorted_bin_tree_qsort (data->comp_tree, comp_cmp_weight);
  pthread_exit (NULL);
  return NULL;
}

/******************* Functions to print stats ********************************/

void
print_stats_uint64 (FILE *out, uint64_t *w, uint64_t len, char name[],
                    char unit[], int verbose)
{
  uint64_t av = 0, min = UMAX(uint64_t), max = 0, std = 0, nb_nzero = 0;
  for (uint64_t i = 0; i < len; i++)
  {
    if (w[i] > 0)
    {
      nb_nzero++;
      if (w[i] < min)
        min = w[i];
      if (w[i] > max)
        max = w[i];
      av += w[i];
      std += w[i]*w[i];
    }
  }

  double av_f = ((double) av) / ((double) nb_nzero);
  double std_f = sqrt(((double) std) / ((double) nb_nzero) - av_f*av_f);

  fprintf (out, "# STATS on %s: #%s = %" PRIu64 "\n", name, name, nb_nzero);
  fprintf (out, "# STATS on %s: min %s = %" PRIu64 "\n", name, unit, min);
  fprintf (out, "# STATS on %s: max %s = %" PRIu64 "\n", name, unit, max);
  fprintf (out, "# STATS on %s: av %s = %.2f\n", name, unit, av_f);
  fprintf (out, "# STATS on %s: std %s = %.2f\n", name, unit, std_f);

  if (verbose > 1)
  {
    uint64_t *nb_w = NULL;
    nb_w = (uint64_t *) malloc ((max-min+1) * sizeof (uint64_t));
    ASSERT_ALWAYS (nb_w != NULL);
    memset (nb_w, 0, (max-min+1) * sizeof (uint64_t));
    for (uint64_t i = 0; i < len; i++)
      if (w[i] > 0)
        nb_w[w[i]-min]++;

    for (uint64_t i = 0; i < max-min+1; i++)
    {
      if (nb_w[i] > 0)
        fprintf (out, "# STATS on %s: #%s of %s %" PRIu64 " : %" PRIu64
                      "\n", name, name, unit, min+i, nb_w[i]);
    }
    free (nb_w);
  }
  fflush (out);
}

void
print_stats_columns_weight (FILE *out, purge_matrix_srcptr mat, int verbose)
{
  uint64_t *w = NULL;
  index_t *h = 0;
  w = (uint64_t *) malloc (mat->col_max_index * sizeof (uint64_t));
  ASSERT_ALWAYS (w != NULL);
  memset (w, 0, mat->col_max_index * sizeof (uint64_t));

  for (uint64_t i = 0; i < mat->nrows_init; i++)
    if (bit_vector_getbit(mat->row_used, (size_t) i))
      for (h = mat->row_compact[i]; *h != UMAX(*h); h++)
        w[*h]++;

  print_stats_uint64 (out, w, mat->col_max_index, "cols", "weight", verbose);
  free (w);
}

void
print_stats_rows_weight (FILE *out, purge_matrix_srcptr mat, int verbose)
{
  uint64_t *w = NULL;
  index_t *h = 0;
  w = (uint64_t *) malloc (mat->nrows_init * sizeof (uint64_t));
  ASSERT_ALWAYS (w != NULL);
  memset (w, 0, mat->nrows_init * sizeof (uint64_t));

  for (uint64_t i = 0; i < mat->nrows_init; i++)
    if (bit_vector_getbit(mat->row_used, (size_t) i))
      for (h = mat->row_compact[i]; *h != UMAX(*h); h++)
        w[i]++;

  print_stats_uint64 (out, w, mat->nrows_init, "rows", "weight", verbose);
  free (w);
}

void print_stats_on_cliques (FILE *out, purge_matrix_srcptr mat, int verbose)
{
  uint64_t *len = NULL;

  len = (uint64_t *) malloc (mat->nrows_init * sizeof (uint64_t));
  ASSERT_ALWAYS (len != NULL);
  memset (len, 0, mat->nrows_init * sizeof (uint64_t));

  uint64_buffer_t buf;
  uint64_buffer_init (buf, UINT64_BUFFER_MIN_SIZE);
  for (uint64_t i = 0; i < mat->nrows_init; i++)
  {
    if (bit_vector_getbit(mat->row_used, (size_t) i))
    {
      comp_t c = {.i = i, .w = 0.0};
      uint64_t nrows = compute_connected_component_mt (&c, mat, buf);
      len[i] = nrows;
    }
  }

  print_stats_uint64 (out, len, mat->nrows_init, "cliques", "length", verbose);

  uint64_buffer_clear (buf);
  free (len);
}

/***************************** Clique removal stage *************************/

static void
cliques_removal (purge_matrix_ptr mat, int64_t target_excess,
                 unsigned int nthreads, int verbose)
{
  int64_t excess = purge_matrix_compute_excess (mat);
  uint64_t nb_clique_deleted = 0;
  size_t i, max_nb_comp_per_thread;

  /* If the excess is smaller than target_excess, then we have nothing to do. */
  if (excess <= target_excess)
    return;

  /* if we are in monothread, we increase max_nb_comp_per_thread by 25% to take
   * into account the fact that something a little bit more connected components
   * are needed to achieve the targeted excess.
   */
  max_nb_comp_per_thread = (size_t) (excess - target_excess);
  if (nthreads == 1)
    max_nb_comp_per_thread+= max_nb_comp_per_thread/4;

  /* First, collect sums for columns with weight 2.
   * If HAVE_SYNC_FETCH is defined, this part is done with the previous
   * multithread function; if not, it's done in sequential, immediatly.
   */
  compute_sum2_row (mat, nthreads);
  fprintf(stdout, "Cliq. rem.: computed mat->sum2_row at %2.2lf\n", seconds());
  fflush (stdout);

  /* Then, each thread searches for its "max_nb_comp_per_thread" heaviest
   * connected compoents and sorts them by decreasing weigh.
   */
  comp_mt_thread_data_t *th_data = (comp_mt_thread_data_t *)
                      malloc_check (nthreads * sizeof(comp_mt_thread_data_t));
  memset (th_data, 0, nthreads * sizeof (comp_mt_thread_data_t));
  for (i = nthreads; i--; )
  {
    th_data[i].th_id = i;
    th_data[i].nthreads = nthreads;
    th_data[i].mat = mat;
    comp_sorted_bin_tree_init (th_data[i].comp_tree, max_nb_comp_per_thread);
    if (pthread_create (&(th_data[i].pthread), NULL,
                        compute_sorted_list_of_connected_components_mt,
                        (void *) &(th_data[i])))
    {
      perror ("compute_sorted_list_of_connected_components_mt\n");
      exit (1);
    }
  }
  for (i = nthreads; i--; )
    pthread_join (th_data[i].pthread, NULL);
  fprintf(stdout, "Cliq. rem.: computed heaviest connected components at "
                  "%2.2lf\n", seconds());
  fflush (stdout);

  if (verbose > 0)
    print_stats_on_cliques (stdout, mat, verbose);

  /* At this point, in each pth[i].comp_tree we have pth[i].comp_tree->size
     connected components order by decreasing weight. */
  size_t *next_clique = NULL;
  uint64_buffer_t buf;
  next_clique = (size_t *) malloc (nthreads * sizeof (next_clique));
  ASSERT_ALWAYS (next_clique != NULL);
  memset (next_clique, 0, nthreads * sizeof (next_clique));
  uint64_buffer_init (buf, UINT64_BUFFER_MIN_SIZE);

  while (mat->nrows > target_excess + mat->ncols)
  {
    comp_t max_comp = {.i = 0, .w = -1.0};
    size_t max_thread = 0;

    for (i = 0; i < nthreads; ++i)
    {
      if (next_clique[i] < th_data[i].comp_tree->size)
      {
        comp_t cur_comp = th_data[i].comp_tree->tree[next_clique[i]];
        if (cur_comp.w > max_comp.w)
        {
          max_comp = cur_comp;
          max_thread = i;
        }
      }
    }
    if (max_comp.w < 0.0)
    {
      fprintf (stderr, "Cliq. rem.: Warning, all lists of connected components"
                       " are empty\n");
      break;
    }

    delete_connected_component (mat, max_comp.i, buf);
    next_clique[max_thread]++;
    nb_clique_deleted++;
  }

  free (next_clique);
  uint64_buffer_clear (buf);
  
  fprintf(stdout, "Cliq. rem.: deleted %" PRIu64 " heaviest connected "
                  "components at %2.2lf\n", nb_clique_deleted, seconds());
  if (verbose > 0)
    fprintf(stdout, "# INFO: max_nb_comp_per_thread=%zu target_excess="
                    "%" PRId64 "\n", max_nb_comp_per_thread, target_excess);
  fflush (stdout);

  /* We can suppress pth[i] and pth itself */
  for (i = 0; i < nthreads; ++i)
    comp_sorted_bin_tree_clear (th_data[i].comp_tree);
  free (th_data);
}

/*****************************************************************************/

/* write final values to stdout */
/* This output, incl. "Final values:", is required by the script */
static void
print_final_values (purge_matrix_ptr mat, double weight)
{
  int64_t excess = purge_matrix_compute_excess (mat);
  fprintf (stdout, "Final values:\nnrows=%" PRIu64 " ncols=%" PRIu64 " "
                   "excess=%" PRId64 "\nweight=%1.0f weight*nrows=%1.2e\n",
                   mat->nrows, mat->ncols, excess, weight,
                   weight * (double) mat->nrows);
  fflush (stdout);
}

/* If nsteps is negative, then the value is chosen by the function. */
static void singletons_and_cliques_removal(purge_matrix_ptr mat, int nsteps,
                                           int64_t final_excess,
                                           double required_excess,
                                           unsigned int nthreads, int verbose)
{
  uint64_t oldnrows = 0;
  int64_t oldexcess, excess, target_excess;
  int count;

  /* First step of singletons removal */
  fprintf(stdout, "\nStep 0: only singleton removal\n");
  excess = singleton_removal (mat, nthreads, verbose);

  if (excess <= 0) /* covers case nrows = ncols = 0 */
  {
    fprintf (stdout, "number of rows <= number of columns\n");
    print_final_values (mat, 0);
    exit(2);
  }

  if ((double) excess < required_excess * ((double) mat->ncols))
  {
    fprintf (stdout, "(excess / ncols) = %.2f < %.2f. See -required_excess "
                     "argument.\n", ((double) excess / (double) mat->ncols),
                     required_excess);
    print_final_values (mat, 0);
    exit(2);
  }

  /* If nsteps was not given in the command line, adjust nsteps in
     [1..DEFAULT_PURGE_NSTEPS] so that each step removes at least about 1% wrt
     the number of columns */
  if (nsteps < 0)
  {
    if ((uint64_t) excess / DEFAULT_PURGE_NSTEPS < mat->ncols / 100)
      nsteps = 1 + (100 * excess) / mat->ncols;
    else
      nsteps = DEFAULT_PURGE_NSTEPS;
  }

  int64_t chunk = excess / nsteps;

  fprintf(stdout, "# INFO: number of clique removal steps: %d\n", nsteps);
  fprintf(stdout, "# INFO: At each step, excess will be decreased by "
                  "%" PRId64 "\n", chunk);
  fflush (stdout);

  /* nsteps steps of clique removal + singletons removal */
  for (count = 0; count < nsteps && excess > 0; count++)
  {
    oldnrows = mat->nrows;
    oldexcess = excess;
    target_excess = excess - chunk;
    if (target_excess < final_excess)
      target_excess = final_excess;
    fprintf(stdout, "\nStep %u on %u: target excess is %" PRId64 "\n",
                    count + 1, nsteps, target_excess);
    fflush(stdout);

    /* prints some stats on columns and rows weight if verbose > 0. */
    if (verbose > 0)
    {
      print_stats_columns_weight (stdout, mat, verbose);
      print_stats_rows_weight (stdout, mat, verbose);
    }

    cliques_removal (mat, target_excess, nthreads, verbose);
    excess = singleton_removal (mat, nthreads, verbose);

    fprintf(stdout, "This step removed %" PRId64 " rows and decreased excess "
                    "by %" PRId64 "\nEach excess row deleted %2.2lf rows\n",
                    (int64_t) (oldnrows-mat->nrows), (oldexcess-excess),
                    (double) (oldnrows-mat->nrows) / (double) (oldexcess-excess));
  }


  /* May need an extra step of clique removal + singletons removal if excess is
     still larger than keep. It may happen due to the fact that each clique does
     not make the excess go down by one but can (rarely) left the excess
     unchanged. */
  if (excess > final_excess && nsteps > 0)
  {
    oldnrows = mat->nrows;
    oldexcess = excess;
    target_excess = final_excess;

    fprintf(stdout, "\nStep extra: target excess is %" PRId64 "\n",
                    target_excess);
    fflush(stdout);

    /* prints some stats on columns and rows weight if verbose > 0. */
    if (verbose > 0)
    {
      print_stats_columns_weight (stdout, mat, verbose);
      print_stats_rows_weight (stdout, mat, verbose);
    }

    cliques_removal (mat, target_excess, nthreads, verbose);
    excess = singleton_removal (mat, nthreads, verbose);

    fprintf(stdout, "This step removed %" PRId64 " rows and decreased excess "
                    "by %" PRId64 "\nEach excess row deleted %2.2lf rows\n",
                    (int64_t) (oldnrows-mat->nrows), (oldexcess-excess),
                    (double) (oldnrows-mat->nrows) / (double) (oldexcess-excess));
  }
}

/*************** Callback functions called by filter_rels ********************/

/* Callback function called by filter_rels on pass 1 */
void *
insert_rel_into_table (purge_matrix_ptr mat, earlyparsed_relation_ptr rel)
{
  ASSERT_ALWAYS(rel->num < mat->nrows_init);
  // TODO: following fct should have only 2 args: mat and rel
  mat->ncols += insert_rel_in_table_no_e(rel, mat->col_min_index, mat->row_compact,
                                         mat->cols_weight);

  return NULL;
}

/* Data struct for thread_print (called by filter_rels on pass 2) */
struct data_second_pass_s
{
  double W; /* Total weight of the matrix (counting only remaining rows) */
  bit_vector_ptr remaining_rows;
  /* fd[0]: for printing kept relations */
  /* fd[1]: for printing deleted relations */
  FILE * fd[2];
};
typedef struct data_second_pass_s data_second_pass_t[1];
typedef struct data_second_pass_s * data_second_pass_ptr;

/* Callback function called by filter_rels on pass 2 */
void *
thread_print(data_second_pass_ptr arg, earlyparsed_relation_ptr rel)
{
  if (bit_vector_getbit(arg->remaining_rows, rel->num))
  {
    arg->W += rel->nb;
    fputs(rel->line, arg->fd[0]);
  }
  else if (arg->fd[1] != NULL)
    fputs(rel->line, arg->fd[1]);
  return NULL;
}


/*********** utility functions for purge binary ****************/

  /* Build the file list (ugly). It is the concatenation of all
   *  b s p
   * where:
   *    b is the basepath (empty if not given)
   *    s ranges over all subdirs listed in the subdirlist (empty if no
   *    such list)
   *    p ranges over all paths listed in the filelist.
   */
char **filelist_from_file_with_subdirlist(const char *basepath,
					  const char *filelist,
					  const char *subdirlist)
{
    /* count the number of files in the filelist */
    int nfiles = 0;
    int nsubdirs = 0;
    char **fl = filelist_from_file(NULL, filelist, 0);
    for (char **p = fl; *p; p++, nfiles++);

    char **sl = filelist_from_file(basepath, subdirlist, 1);
    for (char **p = sl; *p; p++, nsubdirs++);

    char ** fic = (char **) malloc((nsubdirs * nfiles + 1) * sizeof(char *));
    ASSERT_ALWAYS(fic != NULL);

    char **full = fic;
    for (char **f = fl; *f; f++) {
	for (char **s = sl; *s; s++, full++) {
	    int ret = asprintf(full, "%s/%s", *s, *f);
	    ASSERT_ALWAYS(ret >= 0);
	}
    }
    *full = NULL;
    filelist_clear(fl);
    filelist_clear(sl);
    return fic;
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "filelist", "file containing a list of input files");
  param_list_decl_usage(pl, "subdirlist",
                               "file containing a list of subdirectories");
  param_list_decl_usage(pl, "basepath", "path added to all file in filelist");
  param_list_decl_usage(pl, "out", "outfile for remaining relations");
  param_list_decl_usage(pl, "nrels", "number of initial relations");
  param_list_decl_usage(pl, "col-max-index", "upper bound on the number of "
                                  "columns (must be at least the number\n"
                  "                   of prime ideals in renumber table)");
  param_list_decl_usage(pl, "col-min-index", "do not take into account columns"
                                             " with indexes <= col-min-index");
  param_list_decl_usage(pl, "keep", "wanted excess at the end of purge "
                                    "(default " STR(DEFAULT_FILTER_EXCESS) ")");
  param_list_decl_usage(pl, "nsteps", "maximal number of steps of clique "
                                      "removal (default: chosen in [1.."
                                             STR(DEFAULT_PURGE_NSTEPS) "])");
  param_list_decl_usage(pl, "required_excess", "\% of excess required at the "
                            "end of the 1st singleton removal step (default "
                            STR(DEFAULT_PURGE_REQUIRED_EXCESS) ")");
  param_list_decl_usage(pl, "outdel", "outfile for deleted relations (for DL)");
  param_list_decl_usage(pl, "t", "number of threads (default "
                                             STR(DEFAULT_PURGE_NTHREADS) ")");
  param_list_decl_usage(pl, "v", "(switch) verbose mode");
  param_list_decl_usage(pl, "force-posix-threads", "(switch)");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}


/*************************** main ********************************************/

int main(int argc, char **argv)
{
    char * argv0 = argv[0];
    param_list pl;
    uint64_t col_min_index_arg = UMAX(uint64_t);
    char ** input_files;
    uint64_t col_max_index_arg = 0;
    uint64_t nrows_init_arg = 0;
    purge_matrix_t mat; /* All info regarding the matrix is in this struct */
    int64_t keep = DEFAULT_FILTER_EXCESS; /* maximun final excess */
    int nsteps = -1; /* negative value means chosen by purge */
    double required_excess = DEFAULT_PURGE_REQUIRED_EXCESS;
    unsigned int nthreads = DEFAULT_PURGE_NTHREADS;
    int verbose = 0;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;		/* Binary open for all files */
#endif

    double wct0 = wct_seconds();

    param_list_init(pl);
    declare_usage(pl);
    argv++,argc--;

    param_list_configure_switch (pl, "-v", &verbose);
    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

    if (argc == 0)
      usage (pl, argv0);

    /* read all command-line parameters */
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    /* read command-line parameters */
    param_list_parse_uint64(pl, "nrels", &nrows_init_arg);
    param_list_parse_uint64(pl, "col-max-index", &col_max_index_arg);
    param_list_parse_int64(pl, "keep", &keep);

    /* Only look at columns of index >= col-min-index */
    param_list_parse_uint64(pl, "col-min-index", &col_min_index_arg);


    param_list_parse_uint(pl, "t", &nthreads);
    param_list_parse_int(pl, "nsteps", &nsteps);
    param_list_parse_double(pl, "required_excess", &required_excess);

    /* These three parameters specify the set of input files, of the form
     * <base path>/<one of the possible subdirs>/<one of the possible
     * file names>
     *
     * possible subdirs are lister in the file passed as subdirlist.
     * Ditto for possible file names.
     *
     * file names need not be basenames, i.e. they may contain directory
     * components. subdirlist and basepath are optional.
     */
    const char *basepath = param_list_lookup_string(pl, "basepath");
    const char *subdirlist = param_list_lookup_string(pl, "subdirlist");
    const char *filelist = param_list_lookup_string(pl, "filelist");
    const char *purgedname = param_list_lookup_string(pl, "out");
    const char *deletedname = param_list_lookup_string(pl, "outdel");

    if (param_list_warn_unused(pl))
    {
      fprintf(stderr, "Error, unused parameters are given\n");
      usage(pl, argv0);
    }

    if (!purgedname) {
        fprintf(stderr, "Error, option -out is mandatory\n");
        exit(EXIT_FAILURE);
    }
    /* }}} */


    /*{{{ argument checking, and some statistics for things related to
     * the hash table. It needs several static parameters on the command
     * line. This is cumbersome, but while it can probably be avoided, it
     * also hard to do so efficiently */
    if ((basepath || subdirlist) && !filelist)
    {
      fprintf(stderr, "Error, -basepath / -subdirlist only valid with -filelist\n");
      usage(pl, argv0);
    }
    if ((filelist != NULL) + (argc != 0) != 1) {
      fprintf(stderr, "Error, provide either -filelist or freeform file names\n");
      usage(pl, argv0);
    }
    if (nrows_init_arg == 0)
    {
      fprintf(stderr, "Error, missing -nrels command line argument "
                      "(or nrels = 0)\n");
      usage(pl, argv0);
    }
    if (col_max_index_arg == 0)
    {
      fprintf(stderr, "Error, missing -col-max-index command line argument "
                      "(or col-max-index = 0)\n");
      usage(pl, argv0);
    }
    if (col_min_index_arg == UMAX(uint64_t))
    {
      fprintf(stderr, "Error, missing -col-min-index command line argument\n");
      usage(pl, argv0);
    }
    if (col_min_index_arg >= col_max_index_arg)
    {
      fprintf(stderr, "Error, col-min-index >= col-max-index\n");
      usage(pl, argv0);
    }
    /* If col_max_index_arg > 2^32, then we need index_t to be 64-bit */
    if (((col_max_index_arg >> 32) != 0) && sizeof(index_t) < 8)
    {
      fprintf(stderr, "Error, -col-max-index is too large for a 32-bit "
                      "program\nSee #define __SIZEOF_INDEX__ in typedefs.h\n");
      exit(EXIT_FAILURE);
    }
    if (nthreads == 0)
    {
      fprintf(stderr, "Error, -t should be non-zero\n");
      exit(EXIT_FAILURE);
    }

    /* Printing relevant information */
    comp_print_info_weight_function ();

    fprintf(stdout, "# INFO: number of rows: %" PRIu64 "\n", nrows_init_arg);
    fprintf(stdout, "# INFO: maximum possible index of a column: %" PRIu64
                    "\n", col_max_index_arg);
    fprintf(stdout, "# INFO: number of threads: %u\n", nthreads);
    fprintf(stdout, "# INFO: number of clique removal steps: ");
    if (nsteps < 0)
      fprintf(stdout, "will be chosen by the program\n");
    else
      fprintf(stdout, "%d\n", nsteps);
    fprintf(stdout, "# INFO: target excess: %" PRId64 "\n", keep);
    fflush (stdout);
    /*}}}*/

    /* }}} */

    purge_matrix_init (mat, nrows_init_arg, col_min_index_arg,
                       col_max_index_arg);


    set_antebuffer_path(argv0, param_list_lookup_string(pl, "path_antebuffer"));
    /* }}} */

    /*{{{ build the list of input files from the given args
     * If no filelist is given, files are on the command-line.
     * If no subdirlist is given, files are easily construct from basepath and
     * filelist.
     * If subdirlist is given, it is a little bit trickier, see
     * filelist_from_file_with_subdirlist for more details. */
    if (!filelist)
      input_files = argv;
    else if (!subdirlist)
      input_files = filelist_from_file(basepath, filelist, 0);
    else
      input_files = filelist_from_file_with_subdirlist(basepath, filelist, 
                                                       subdirlist);
    /*}}}*/

    /****************** Begin interesting stuff *************************/
    fprintf(stdout, "\nPass 1, reading and storing columns with index h >= "
                    "%" PRIu64 "\n", mat->col_min_index);

    /* first pass over relations in files */
    /* Note: Now that we no longer take a bitmap on input, all
     * relations are considered active at this point, so that we
     * do not need to pass a bitmap to filter_rels */
    mat->nrows = filter_rels(input_files,
                        (filter_rels_callback_t) &insert_rel_into_table, mat,
                        EARLYPARSE_NEED_INDEX, NULL, NULL);

    if (mat->nrows != mat->nrows_init)
    {
      fprintf(stderr, "Error, -nrels value should match the number of scanned "
                      "relations\nexpected %" PRIu64 " relations, found "
                      "%" PRIu64 "\n", mat->nrows_init, mat->nrows);
      abort();
    }

    /* Take into account the memory allocated for all mat->row_compact[i] */
    purge_matrix_clear_row_compact_update_mem_usage (mat);

    /* prints some stats on columns and rows weight if verbose > 0. */
    if (verbose > 0)
    {
      print_stats_columns_weight (stdout, mat, verbose);
      print_stats_rows_weight (stdout, mat, verbose);
    }

    /* MAIN FUNCTIONS: do singletons and cliques removal. */
    singletons_and_cliques_removal (mat, nsteps, keep, required_excess,
                                    nthreads, verbose);

    /* prints some stats on columns and rows weight if verbose > 0. */
    if (verbose > 0)
    {
      print_stats_columns_weight (stdout, mat, verbose);
      print_stats_rows_weight (stdout, mat, verbose);
    }

    /* XXX: Are these tests useful ? Already checked after first call to
     * removal_all_singletons in singletons_and_cliques_removal. */
    if (mat->nrows < mat->ncols)
    {
      fprintf (stdout, "number of rows <= number of columns\n");
      print_final_values (mat, 0);
      exit(2);
    }
    if (mat->nrows == 0 || mat->ncols == 0)
    {
      fprintf(stdout, "number of rows or number of columns is 0\n");
      print_final_values (mat, 0);
      exit(2);
    }

    /* free mat->row_compact[i] and mat->row_compact. We no longer need them */
    purge_matrix_clear_row_compact (mat);

    /****** Pass 2: reread the relation files and write output file(s) ******/
    fprintf(stdout, "\nPass 2, reading and writing output file%s...\n",
                    deletedname == NULL ? "" : "s");
    data_second_pass_t data2;
    memset(data2, 0, sizeof(data_second_pass_t));
    data2->remaining_rows = mat->row_used;

    if (!(data2->fd[0] = fopen_maybe_compressed(purgedname, "w")))
    {
      fprintf(stderr, "Error, cannot open file %s for writing.\n", purgedname);
      exit(1);
    }

    if (deletedname != NULL)
    {
      if (!(data2->fd[1] = fopen_maybe_compressed(deletedname, "w")))
      {
        fprintf(stderr, "Error, cannot open file %s for writing.\n",
                        deletedname);
        exit(1);
      }
      /* Write the header line for the file of deleted relations. */
      fprintf(data2->fd[1], "# %" PRIu64 "\n", mat->nrows_init - mat->nrows);
    }

    /* Write the header line for the file of remaining relations:
     * compute last index i such that cols_weight[i] != 0
     */
    {
      uint64_t last_used = mat->col_max_index - 1;
      while (mat->cols_weight[last_used] == 0)
        last_used--;

      fprintf(data2->fd[0], "# %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
                            mat->nrows, last_used + 1, mat->ncols);
    }

    /* second pass over relations in files */
    filter_rels(input_files, (filter_rels_callback_t) &thread_print, data2,
                EARLYPARSE_NEED_LINE, NULL, NULL);

    /* write final values to stdout */
    /* This output, incl. "Final values:", is required by the script */
    print_final_values (mat, data2->W);

    /* Free allocated stuff */
    if (filelist)
      filelist_clear(input_files);

    fclose_maybe_compressed(data2->fd[0], purgedname);
    if (data2->fd[1])
      fclose_maybe_compressed(data2->fd[1], deletedname);

    purge_matrix_clear (mat);
    /* print usage of time and memory */
    print_timing_and_memory(wct0);

    param_list_clear(pl);

    return 0;
}
