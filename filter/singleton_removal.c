#include "cado.h"
#include <stdio.h>
#include <stdlib.h>

#include "portability.h"

#include "filter_common.h"
#include "purge_matrix.h"
#include "singleton_removal.h"

/* If HAVE_SYNC_FETCH is not defined, we will use mutex for multithreaded
 * version of the code. May be too slow. */
#ifndef HAVE_SYNC_FETCH
pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER;
#endif

/* Code for singletons removal --- monothread version */

void singleton_removal_oneiter_mono (purge_matrix_ptr mat)
{
  for (uint64_t i = 0; i < mat->nrows_init; i++)
  {
    if (bit_vector_getbit(mat->row_used, (size_t) i))
    {
      index_t *row;
      for (row = mat->row_compact[i]; *row != UMAX(*row); row++)
      {
        if (mat->cols_weight[*row] == 1)
        {
          /* mat->nrows and mat->ncols are updated by purge_matrix_delete_row.*/
          purge_matrix_delete_row (mat, i);
          break;
        }
      }
    }
  }
}


/* Code for singletons removal --- multithread version */

typedef struct sing_rem_mt_data_s
{
  /* Read only */
  uint64_t begin, end;
  /* Read / Write */
  uint64_t sup_nrow, sup_ncol;
  purge_matrix_ptr mat;
} sing_rem_mt_data_t;

/* Hightest criticality for performance. I inline all myself. */
static void *
singleton_removal_mt_thread (void *arg)
{
  sing_rem_mt_data_t *data = (sing_rem_mt_data_t *) arg;
  purge_matrix_ptr mat = data->mat;
  index_t *row;
  weight_t *w;
  bv_t j;

  data->sup_nrow = data->sup_ncol = 0;
  for (uint64_t i = data->begin; i < data->end; i++)
  {
    j = (((bv_t) 1) << (i & (BV_BITS - 1)));

    if (mat->row_used->p[i >> LN2_BV_BITS] & j) /* Still active ? */
    {
      for (row = mat->row_compact[i]; *row != UMAX(*row); row++)
      {
        if (UNLIKELY(mat->cols_weight[*row] == 1)) /* We found a singleton */
        {
          for (row = mat->row_compact[i]; *row != UMAX(*row); row++)
          {
            w = &(mat->cols_weight[*row]);
            ASSERT(*w);
            /* Decrease only if not equal to the maximum value */
            /* If weight becomes 0, we just remove a column */
#ifdef HAVE_SYNC_FETCH
            if (*w < UMAX(*w) && !__sync_sub_and_fetch(w, 1))
              (data->sup_ncol)++;
#else /* else we use mutex to protect the substraction on w */
            pthread_mutex_lock (&lock);
            if (*w < UMAX(*w) && !(--(*w)))
              (data->sup_ncol)++;
            pthread_mutex_unlock (&lock);
#endif
          }
          /* We do not free mat->row_compact[i] as it is freed later with
             my_malloc_free_all */
          /* This is thread-safe because we know that all rows in
             data->row_used->p[i >> LN2_BV_BITS] are handled by the same thread
           */
          mat->row_used->p[i >> LN2_BV_BITS] &= ~j;
          (data->sup_nrow)++;
          break;
        }
      }
    }
  }
  return NULL;
}

void
singleton_removal_oneiter_mt (purge_matrix_ptr mat, unsigned int nthreads)
{
  pthread_attr_t attr;
  pthread_t *threads;
  sing_rem_mt_data_t *th_data;
  uint64_t nrows_per_thread, k;
  int err;

  th_data = (sing_rem_mt_data_t *) malloc (nthreads*sizeof(sing_rem_mt_data_t));
  ASSERT_ALWAYS(th_data != NULL);
  threads = (pthread_t *) malloc (nthreads * sizeof(pthread_t));
  ASSERT_ALWAYS(threads != NULL);

  /* nrows_per_thread MUST be a multiple of BV_BITS to ensure that each address
   * of the bit_vector array is accessed by only one thread. */
  nrows_per_thread = (mat->nrows_init / nthreads) & ((uint64_t) ~(BV_BITS - 1));
  k = 0;
  for (unsigned int i = 0; i < nthreads; i++)
  {
    th_data[i].mat = mat;
    th_data[i].begin = k;
    k += nrows_per_thread;
    th_data[i].end = k;
  }
  th_data[nthreads-1].end = mat->nrows_init;

  pthread_attr_init (&attr);
  pthread_attr_setstacksize (&attr, 1 << 16);
  pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
  for (unsigned int i = 0; i < nthreads; i++)
  {
    if ((err = pthread_create(&threads[i], &attr, &singleton_removal_mt_thread,
                              (void *) &(th_data[i]))))
    {
      fprintf(stderr, "Error, pthread_create failed in singleton_removal_"
                      "oneiter_mt: %d. %s\n", err, strerror(errno));
      abort();
    }
  }
  for (unsigned i = 0; i < nthreads; i++)
  {
    pthread_join (threads[i], NULL);
    mat->nrows -= th_data[i].sup_nrow;
    mat->ncols -= th_data[i].sup_ncol;
  }
  pthread_attr_destroy(&attr);
  free(threads);
  free(th_data);
}

/* Perform a complete singleton removal step:  call singleton_removal_oneiter_*
 * until there is no more singleton.
 * Return the excess at the end *
 */
int64_t
singleton_removal (purge_matrix_ptr mat, unsigned int nthreads, int verbose)
{
  int64_t excess;
  uint64_t oldnrows;
  unsigned int iter = 0;
  do
  {
    oldnrows = mat->nrows;
    excess = purge_matrix_compute_excess (mat);
    if (verbose >= 0) /* if no quiet */
    {
      if (iter == 0)
        fprintf(stdout, "Sing. rem.: begin with: ");
      else
        fprintf(stdout, "Sing. rem.:   iter %03u: ", iter);
      fprintf(stdout, "nrows=%" PRIu64 " ncols=%" PRIu64 " excess=%" PRId64 " "
                      "at %2.2lf\n", mat->nrows, mat->ncols, excess, seconds());
      fflush(stdout);
    }

  if (nthreads > 1)
    singleton_removal_oneiter_mt (mat, nthreads);
  else
    singleton_removal_oneiter_mono (mat);

    iter++;
  } while (oldnrows != mat->nrows);

  if (verbose >= 0) /* if no quiet */
    fprintf(stdout, "Sing. rem.:   iter %03u: No more singletons, finished at "
                    "%2.2lf\n", iter, seconds());

  return excess;
}
