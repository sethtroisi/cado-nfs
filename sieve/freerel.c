/*
 * Program: free relations
 * Original author : F. Morain
 * Purpose: creating free relations in a suitable format
 * Modified / rewritten by C. Bouvier (and others)
 * Multi-thread code by A. Filbois

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* Model of this code : One producer -> Many consumers/producers -> one consumer.
   1. First thread produces SIZE_BUF_PRIMES primes by buffer, in NB_BUFS buffers
      (named primes) for each working thread.

      When all the primes are produced, this thread produces one buffer with
      only one "false" prime = (p_r_values_t) (-1) for each pthread consumer.
      It's the end marker.

   2. nb_pthread threads load a primes buffer by a classical one producer/many
      consumer models (2 pseudo semaphores BY pair of producer/consumer).
      Each pthread produces 2 buffers: one buffer of free_relations, each on the
      form of prime + a set of (deg(polynom1)+deg(polynom2)+1 free_rels_buf_t
      type (unsigned long or unsigned int); and one buffer of roots, on the form
      of a ASCII array.

      When a thread loads a primes buffer which begins and contains only
      (p_r_values_t) (-1), the thread produces a roots buffer which contains
      only (p_r_values_t) (-1) and an empty free relations buffer, and exits.

   3. With a many producers/one consumer model, the principal programs loads the
      roots & free_buffers, writes sequentially the roots, computes the real
      renumber index (sum of all previous index) for the free relations, and
      writes these relations.

      When a roots buffer contains exactly one (p_r_values_t) (-1), the job is
      done.
*/


#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>   /* for _O_BINARY */
#include <string.h>
#include <math.h>
#include <pthread.h>

#include <gmp.h>
#include "mod_ul.c"
#include "portability.h"
#include "utils.h"
#include "typedefs.h"

/* 2 buffers are sufficient; in fact it works with only one! */
#define NB_BUFFERS_PER_THREAD (1<<1)
#define NB_PRIMES_PER_BUFFER (1<<10)
#define CACHELINESIZE 64
#define INIT_SIZE_BUF_CHAR (1<<15) /* Init size for the char buffer. */
/* Max number of char necessary to write all the roots in hexa for ONE prime p.
 * Must be greater than (max_nb_char_per_roots)*MAXDEGREE*NB_POLYS_MAX.
 * max_nb_char_per_roots is 8 (a root is always smaller than a prime and a
 * prime is smaller than 2^64 (this bound is not tight at all, in practice a
 * prime is smaller 2^32~2^34 at most))
 */
#define MAX_SIZE_PER_PRIME (1<<10)
#if MAX_SIZE_PER_PRIME < 8*MAXDEGREE*NB_POLYS_MAX
  #error "MAX_SIZE_PER_PRIME is too small."
#endif
#if INIT_SIZE_BUF_CHAR < MAX_SIZE_PER_PRIME
  #error "INIT_SIZE_BUF_CHAR must be greater than MAX_SIZE_PER_PRIME."
#endif

/************************** Unnamed semaphores ********************************/
/* Unnamed Posix semaphores don't exist in MacOS, named Posix semaphores are
 * persistent and need a interrupt handler to clean them, so unnamed semaphores
 * are implemented here with pthread_cond_t and pthread_mutex_t, and of course a
 * value (unsigned int).
 */
struct sema_s
{
  unsigned int value;
  pthread_mutex_t mutex;
  pthread_cond_t conditional;
};

typedef struct sema_s sema_t[1];
typedef struct sema_s * sema_ptr;

/* These functions are sem_init, sem_post, sem_wait and sem_destroy with mutex &
 * pthread conditional variables.
 * Careful:
 *  1. it's not possible to share these pseudo semaphores between processus, but
 *     only between threads.
 *  2. I use pthread_cond_signal here, not pthread_cond_broadcast! So sema_post
 *     wakes only one sleeper!
 */

void sema_init (sema_ptr sema, unsigned int value_init)
{
  sema->value = value_init;
  if (pthread_mutex_init(&(sema->mutex), NULL))
  {
    perror ("Error in pthread_mutex_init\n");
    exit (1);
  }
  if (pthread_cond_init(&(sema->conditional), NULL))
  {
    perror ("Error in pthread_cond_init\n");
    exit (1);
  }
}

void sema_post (sema_ptr sema)
{
  pthread_mutex_lock (&(sema->mutex));
  if (!(sema->value)++)
    pthread_cond_signal (&(sema->conditional)); /* CAREFUL: no broadcast here! */
  pthread_mutex_unlock (&(sema->mutex));
}

void sema_wait (sema_ptr sema)
{
  pthread_mutex_lock (&(sema->mutex));
  /* The loop is need for the "spurious wakeup", but these is only at most 1
   * sleeper by sema... so, not really useful. */
  while (!(sema->value))
    pthread_cond_wait (&(sema->conditional), &(sema->mutex));
  sema->value--;
  pthread_mutex_unlock (&(sema->mutex));
}

void sema_destroy (sema_ptr sema)
{
  pthread_mutex_destroy (&(sema->mutex));
  pthread_cond_destroy (&(sema->conditional));
}

/************************** Thread buffer *************************************/
struct th_buf_s
{
  /* protected by primes_full/primes_empty */
  unsigned long primes[NB_PRIMES_PER_BUFFER]; /* array of prime number */
  uint64_t nprimes_in; /* nb of primes in the array on input */
  /* protected by roots_full/roots_empty */
  uint64_t nprimes_out; /* nb of treated primes */
  uint64_t nfreerels;
  uint64_t size_local_renum_tab; /* nb of entry in the local table */
  unsigned long freerels_p[NB_PRIMES_PER_BUFFER];
  uint64_t freerels_first_index[NB_PRIMES_PER_BUFFER];
  char *local_renum_tab; /* The local renumbering table is kept as a char*. */
  char *local_renum_tab_cur;
  uint64_t local_renum_tab_alloc;
};

typedef struct th_buf_s th_buf_t[1];
typedef struct th_buf_s * th_buf_ptr;
typedef const struct th_buf_s * th_buf_srcptr;

static inline void
th_buf_init (th_buf_ptr buf, size_t size)
{
  buf->local_renum_tab = malloc_aligned (size * sizeof (char), CACHELINESIZE);
  buf->local_renum_tab_alloc = size;
}

static inline void
th_buf_reset (th_buf_ptr buf)
{
  buf->nprimes_out = 0;
  buf->nfreerels = 0;
  buf->size_local_renum_tab = 0;
  buf->local_renum_tab_cur = buf->local_renum_tab;
}

/* This function grows the char * if the remaining space is lower than
 * min_needed.
 */
static inline void
th_buf_char_resize (th_buf_ptr buf, size_t min_needed)
{
  char *end = buf->local_renum_tab + buf->local_renum_tab_alloc;
  if (UNLIKELY(buf->local_renum_tab_cur + min_needed >= end))
  {
    size_t save_current = buf->local_renum_tab_cur - buf->local_renum_tab,
    new_size = buf->local_renum_tab_alloc << 1;
    buf->local_renum_tab = (char *) realloc (buf->local_renum_tab,
                                                      new_size * sizeof (char));
    ASSERT_ALWAYS (buf->local_renum_tab != NULL);
    buf->local_renum_tab_cur = buf->local_renum_tab + save_current;
    buf->local_renum_tab_alloc = new_size;
  }
}

static inline void
th_buf_clear (th_buf_ptr buf)
{
  free (buf->local_renum_tab);
}

/******************** Data struct for sharing among threads *******************/
/* The main structure for the working pthreads pool */
typedef struct freerel_th_data_s {
  /* 4 pseudo semaphores. roots_* are also used for free_rel buffer. */
  sema_t primes_full, primes_empty, roots_full, roots_empty;
  /* Read only part */
  unsigned int nthreads;
  unsigned int th_id;
  unsigned long pmin, pmax, lpb[NB_POLYS_MAX], lpbmax;
  mpz_poly_ptr pols[NB_POLYS_MAX];
  int nb_polys;
  /* Read-write part */
  th_buf_t bufs[NB_BUFFERS_PER_THREAD];
} freerel_th_data_t;

/*************************** Main functions ***********************************/

/* This function produces all the primes from 2 to pth->lpbmax,
   for the workers pthreads pool (see next functions). The primes are in
   buffers; each buffer contains SIZE_BUF_PRIMES primes (optimal: 1024),
   except the last one of course.
   After all buffers have been produced, this function creates again
   nb_pthreads buffers with only the end marker ((p_r_values_t) (-1)),
   in order to stop each worker.
*/
static void *
pthread_primes_producer (void *arg)
{
  /* This thread get the array of all the thread's data. */
  freerel_th_data_t *data = (freerel_th_data_t *) arg;
  unsigned long lpbmax = data->lpbmax, p = 2;
  size_t nthreads = data->nthreads, cur_buf = 0, cur_th = 0;

  prime_info pi;
  prime_info_init (pi);
  while (p <= lpbmax)
  {
    /* Here, we are sure that we can write at least one prime. */
    sema_wait (data[cur_th].primes_empty); /* Need an empty primes buffer */
    th_buf_ptr local_buf = data[cur_th].bufs[cur_buf];
    local_buf->nprimes_in = 0;

    do /* Main loop that produces primes */
    {
      local_buf->primes[local_buf->nprimes_in++] = p;
      p = getprime_mt (pi); /* get next prime */
    } while (local_buf->nprimes_in < NB_PRIMES_PER_BUFFER && p <= lpbmax);

    sema_post (data[cur_th].primes_full); /* Drop the primes buffer */

    if (UNLIKELY (++cur_th == nthreads))
    {
      cur_th = 0;
      if (UNLIKELY (++cur_buf == NB_BUFFERS_PER_THREAD))
        cur_buf = 0;
    }
  }
  prime_info_clear (pi);

  /* We have to produce a special end buffer in each pth[].primes */
  size_t bak_th = cur_th;
  do
  {
    sema_wait (data[cur_th].primes_empty);
    th_buf_ptr local_buf = data[cur_th].bufs[cur_buf];
    local_buf->primes[0] = (unsigned long) (-1);
    local_buf->nprimes_in = 1;
    sema_post (data[cur_th].primes_full);

    if (UNLIKELY (++cur_th == nthreads))
    {
      cur_th = 0;
      if (UNLIKELY (++cur_buf == NB_BUFFERS_PER_THREAD))
        cur_buf = 0;
    }
  } while (cur_th != bak_th);

  return NULL;
}

/* This function takes a primes buffer and generates the corresponding
   roots and free relations.
   The roots are directly on the form of a ASCII buffer.
   The free relations must be recomputed by the main program, so a free relation
   contains only its prime and ((degree (polynom1) + degree (polynom2)) incremental
   index.
   This functions stops when the primes buffer contains only ((p_r_values_t) (-1)),
   the end marker.
*/
static void *
pthread_roots_and_free_rels_producer (void *arg)
{
  freerel_th_data_t *my_data = arg;
  unsigned long p;
  int nb_roots[NB_POLYS_MAX];
  size_t cur_buf = 0;
  unsigned long computed_roots[NB_POLYS_MAX][32]; // With malloc, the computed_roots of 2 pthreads
  // may be in the same cacheline.
  // A computed_root with a fixed size is faster.
  for(int k = 0; k < my_data->nb_polys; k++)
      ASSERT_ALWAYS (my_data->pols[k]->deg < 32);

  for (;;)
  {
    sema_wait (my_data->roots_empty); // Need an empty roots & free relations buffer
    sema_wait (my_data->primes_full); // Need a produced primes buffer

    th_buf_ptr local_buf = my_data->bufs[cur_buf];
    unsigned long *my_p = local_buf->primes;

    // Have we received the special end buffer in primes buffer ?
    if (UNLIKELY(local_buf->nprimes_in == 1 && my_p[0] == (unsigned long) (-1)))
    {
      // Yes: we create a special end buffer in roots buffer, and exit
      local_buf->freerels_p[0] = (unsigned long) (-1);
      local_buf->nfreerels = 1;
      sema_post (my_data->primes_empty); // Drop my now consumed primes buffer (not useful)
      sema_post (my_data->roots_full);   // Give my produced roots & free_rels buffer
      return NULL; /* The job is done: exit the endless loop and the pthread. */
    }

    th_buf_reset (local_buf);

    // I prefer to separate the 2 cases (rat & alg) and (alg & alg).
    // The goal is to suppress the first case in few months
    // TODO: make this MNFS-compliant!
    if (my_data->pols[0]->deg == 1)
    {
      // rational & algebraic polynomials, special case
      while (local_buf->nprimes_out < local_buf->nprimes_in)
      {
        p = *my_p++;  // p is the current prime
        // First, we compute the roots of alg
        if (UNLIKELY(p > my_data->lpb[1]))
        {
          nb_roots[0] = 1;
          nb_roots[1] = 0;
        }
        else
        {
          nb_roots[0] = (p < my_data->lpb[0]) ? 1 : 0;
          nb_roots[1] = mpz_poly_roots_ulong (computed_roots[1], my_data->pols[1], p);
          if (UNLIKELY(nb_roots[1] != my_data->pols[1]->deg &&
              mpz_divisible_ui_p (my_data->pols[1]->coeff[my_data->pols[1]->deg], p)))
          {
            computed_roots[1][nb_roots[1]++] = computed_roots[1][0];
            computed_roots[1][0] = p; // p is inserted in first place (for the next sort) as a root
          }
        }

        /* Second, we fill my_roots->current buffer by the computed roots.
         * We are sure that there is enough space to write all the roots for one
         * prime p, but once this is finished, we must check if the size of the
         * buffer need to be increased.
         */
        local_buf->local_renum_tab_cur +=
          renumber_write_p_buffer_rat_alg (local_buf->local_renum_tab_cur, p,
                                           nb_roots[0], computed_roots[1],
                                           nb_roots[1]);
        th_buf_char_resize (local_buf, MAX_SIZE_PER_PRIME);

    
        // Third, we fill my_free_rels->current buffer by the possible free rels, on the form
        // [1] p
        // [(degree of alg + 1)] ++(my_free_rels->renumber)
        if (UNLIKELY(nb_roots[0] && nb_roots[1] == my_data->pols[1]->deg &&
             p >= my_data->pmin && p <= my_data->pmax))
        {
          local_buf->freerels_p[local_buf->nfreerels] = p;
          local_buf->freerels_first_index[local_buf->nfreerels] =
                                              local_buf->size_local_renum_tab;
          local_buf->nfreerels++;
        }
        local_buf->size_local_renum_tab += nb_roots[0] + nb_roots[1];
        local_buf->nprimes_out++;
      } /* Next p in current primes buffer */
    }
    else
    {
      // algebraic polynomials poly, "normal" case
      while (local_buf->nprimes_out < local_buf->nprimes_in)
      {
        p = *my_p++;  // p is the current prime
        // First, we compute the roots of all polys
        for (size_t my_alg = 0; my_alg < (size_t)my_data->nb_polys; ++my_alg)
        {
          if (LIKELY (p < my_data->lpb[my_alg]))
          {
            nb_roots[my_alg] = mpz_poly_roots_ulong (computed_roots[my_alg],
                                                     my_data->pols[my_alg], p);
            if (UNLIKELY(nb_roots[my_alg] != my_data->pols[my_alg]->deg &&
                mpz_divisible_ui_p (my_data->pols[my_alg]->coeff[my_data->pols[my_alg]->deg], p)))
            {
              computed_roots[my_alg][nb_roots[my_alg]++] = computed_roots[my_alg][0];
              computed_roots[my_alg][0] = p; // p is inserted in first place (for the next sort) as a root
            }
          }
          else
            nb_roots[my_alg] = 0;
        }
    
        /* Second, we fill my_roots->current buffer by the computed roots.
         * We are sure that there is enough space to write all the roots for one
         * prime p, but once this is finished, we must check if the size of the
         * buffer need to be increased.
         */
        local_buf->local_renum_tab_cur +=
          renumber_write_p_buffer_2algs (local_buf->local_renum_tab_cur, p,
                       computed_roots[0], nb_roots[0],
                       computed_roots[1], nb_roots[1]);
        th_buf_char_resize (local_buf, MAX_SIZE_PER_PRIME);

        // Third, we fill my_free_rels->current buffer by the possible free rels, on the form
        // [1] : p
        // [(degree of alg + 1)] : ++(my_free_rels->renumber)
        if (UNLIKELY(nb_roots[0] == my_data->pols[0]->deg &&
                     nb_roots[1] == my_data->pols[1]->deg &&
                     p >= my_data->pmin && p <= my_data->pmax))
        {
          local_buf->freerels_p[local_buf->nfreerels] = p;
          local_buf->freerels_first_index[local_buf->nfreerels] =
                                              local_buf->size_local_renum_tab;
          local_buf->nfreerels++;
        }
        local_buf->size_local_renum_tab += nb_roots[0] + nb_roots[1];
        local_buf->nprimes_out++;
      } /* Next p in current primes buffer */
    }

    sema_post(my_data->primes_empty); // Drop my now consumed primes buffer
    sema_post(my_data->roots_full);   // Give my produced roots & free_rels buffer
    if (UNLIKELY (++(cur_buf) == NB_BUFFERS_PER_THREAD)) /* Go to next buffer */
      cur_buf = 0;
  } // Endless loop
}

/* generate all free relations up to the large prime bound */
/* generate the renumbering table */
static unsigned long
allFreeRelations (cado_poly pol, unsigned long pmin, unsigned long pmax,
                  unsigned long lpb[NB_POLYS_MAX], renumber_t renumber_table,
                  size_t nthreads, const char *outfilename)
{
  ASSERT_ALWAYS(pol->nb_polys == 2); // TMP!!!!!

  /* open outfile */
  FILE *outfile = NULL;
  outfile = fopen_maybe_compressed (outfilename, "w");
  ASSERT_ALWAYS (outfile != NULL);

  /* precompute sum(pols[k]->deg for k in [1..nb_polys]) */
  unsigned int sum_degs = 0; /* sum of the degrees of all polynomials */
  for(int k = 0; k < pol->nb_polys; k++)
      sum_degs += pol->pols[k]->deg;

  /* Set lbpmax lbpmin and handle pmax and pmin.
   * We generate all free relations from pmin and up to the *minimum*
   * of the two large prime bounds, since larger primes will never
   * occur on both sides.
   * We generate the renumbering table from the first prime (2) and
   * up to the *maximum* of the two large prime bounds.
   */
  unsigned long lpbmax, lpbmin; /* MAX(lpb[0],lpb[1],...), MIN(..) */
  lpbmin = 1UL << lpb[0]; lpbmax = 0;
  for(int k = 0; k < pol->nb_polys; k++)
  {
    ASSERT_ALWAYS (lpb[k] < sizeof(unsigned long) * CHAR_BIT);
    lpb[k] = 1UL << lpb[k];
    lpbmax = MAX(lpbmax, lpb[k]);
    lpbmin = MIN(lpbmin, lpb[k]);
  }

  if (pmax && pmax > lpbmin)
  {
    fprintf (stderr, "Error: pmax is greater than MIN(lpb[])\n");
    exit (1);
  }
  else if (!pmax)
    pmax = lpbmin;

  printf ("Generating freerels for %lu <= p <= %lu\n", pmin, pmax);
  printf ("Generating renumber table for 2 <= p <= %lu\n", lpbmax);
  fflush (stdout);

  /* Init statistics */
  uint64_t nprimes_total = 0; /* Total number of the primes */
  uint64_t nfreerels_total = 0; /* Total number of the free relations */
  stats_data_t stats; /* struct for printing progress */
  /* will print report at 2^10, 2^11, ... 2^23 computed primes and every
   * 2^23 primes after that */
  stats_init (stats, stdout, &nprimes_total, 23, "Treated", "primes", "", "p");

  /* Init threads and threads data */
  freerel_th_data_t *th_data = NULL; /* data of the threads */
  pthread_t primes_producer_th; /* pthread for primes producer */
  pthread_t *primes_consumer_th = NULL;

  th_data = (freerel_th_data_t *) malloc (nthreads * sizeof (freerel_th_data_t));
  ASSERT_ALWAYS (th_data != NULL);
  primes_consumer_th = (pthread_t *) malloc (nthreads * sizeof (pthread_t));
  ASSERT_ALWAYS (primes_consumer_th != NULL);

  for (size_t i = 0; i < nthreads; i++)
  {
    sema_init (th_data[i].roots_full, 0);
    sema_init (th_data[i].roots_empty, NB_BUFFERS_PER_THREAD);
    sema_init (th_data[i].primes_full, 0);
    sema_init (th_data[i].primes_empty, NB_BUFFERS_PER_THREAD);
    th_data[i].th_id = i;
    th_data[i].nthreads = nthreads;
    th_data[i].pmin = pmin;
    th_data[i].pmax = pmax;
    th_data[i].lpbmax = lpbmax;
    th_data[i].nb_polys = pol->nb_polys;
    for(int k = 0; k < pol->nb_polys; k++)
    {
      th_data[i].lpb[k] = lpb[k];
      th_data[i].pols[k] = pol->pols[k];
    }
    for (size_t j = 0; j < NB_BUFFERS_PER_THREAD; j++)
      th_buf_init (th_data[i].bufs[j], INIT_SIZE_BUF_CHAR);
  }

  if (pthread_create (&(primes_producer_th), NULL, &pthread_primes_producer,
                      (void *) th_data))
  {
    perror ("pthread_primes_producer creation failed\n");
    exit (1);
  }
  for (size_t i = 0; i < nthreads; i++)
  {
    if (pthread_create (&(primes_consumer_th[i]), NULL,
                        &pthread_roots_and_free_rels_producer,
                        (void *) &(th_data[i])))
    {
      perror ("pthread_roots_and_free_rels_producer pthread creation failed\n");
      exit (1);
    }
  }

  size_t cur_th = 0;
  size_t cur_buf = 0;
  // Main loop: load the roots & free rels buffers, and print them
  for (;;)
  {
    sema_wait (th_data[cur_th].roots_full); // Need a full root & free rels buffer
    th_buf_ptr local_buf = th_data[cur_th].bufs[cur_buf];

    // Is it the special end buffer (with only the end marker ?)
    if (local_buf->nfreerels == 1 &&
        local_buf->freerels_p[0] == (unsigned long) (-1))
    {
      /* All is done: when a pthread which produces roots & free rels gives this
	 special end buffer, all the next pthreads have finished and give also another
	 special end buffer. So it's not useful to read all; one is sufficient.
	 But it's dirty, so I prefer "read" (consume the semaphores) them all.
	 Morever, some ASSERTs after the loop verify the values of the pseudo semaphores;
	 if all is OK, these values must be equal to their initializations.
      */
      sema_post (th_data[cur_th].roots_empty);
      for (size_t th_bak = cur_th;;)
      {
        if (UNLIKELY(++cur_th == nthreads))
          cur_th = 0;
        if (UNLIKELY (cur_th == th_bak))
          break;
        sema_wait (th_data[cur_th].roots_full);
        /* Nothing to do after this, we drop directly this special end buffer. */
        sema_post (th_data[cur_th].roots_empty);
      }
      break; // End of the main loop
    }

    /* We write the roots in one ASCII block */
    size_t n = local_buf->local_renum_tab_cur - local_buf->local_renum_tab;
    fwrite ((void *) local_buf->local_renum_tab, n * sizeof (char), 1,
            renumber_table->file);

    // We have to recompute the real index of the renumber table for the free rels
    for (uint64_t k = 0; k < local_buf->nfreerels; k++)
    {
      uint64_t index = local_buf->freerels_first_index[k] + renumber_table->size;
      uint64_t last_index = index + sum_degs;
      fprintf (outfile, "%lx,0:%" PRIx64, local_buf->freerels_p[k], index);
      index++;   
      for (; index < last_index; index++)
        fprintf (outfile, ",%" PRIx64, index);
      fputc ('\n', outfile);
    }
    renumber_table->size += local_buf->size_local_renum_tab;
    nfreerels_total += local_buf->nfreerels;
    nprimes_total += local_buf->nprimes_out;

    /* Drop the consumed root & free rels buffer */
    sema_post (th_data[cur_th].roots_empty);

    if (stats_test_progress(stats))
      stats_print_progress (stats, nprimes_total, 0, 0, 0);

    if (UNLIKELY (++cur_th == nthreads))
    {
      cur_th = 0;
      if (UNLIKELY (++cur_buf == NB_BUFFERS_PER_THREAD))
        cur_buf = 0;
    }
  }

  // All is done!
  stats_print_progress (stats, nprimes_total, 0, 0, 1);
  if (pthread_join (primes_producer_th, NULL))
  {
    perror ("Error, pthread primes producer stops abnormally\n");
    exit (1);
  }
  for (size_t i = 0; i < nthreads; i++)
  {
    if (pthread_join (primes_consumer_th[i], NULL))
    {
      perror ("Error, one of the threads pthread_roots_and_free_rels_producer "
              "stops abnormally\n");
      exit (1);
    }
    ASSERT(th_data[i].roots_full->value == 0);
    ASSERT(th_data[i].roots_empty->value == NB_BUFFERS_PER_THREAD);
    ASSERT(th_data[i].primes_full->value == 0);
    ASSERT(th_data[i].primes_empty->value == NB_BUFFERS_PER_THREAD);
    sema_destroy (th_data[i].roots_full);
    sema_destroy (th_data[i].roots_empty);
    sema_destroy (th_data[i].primes_full);
    sema_destroy (th_data[i].primes_empty);
    for (size_t j = 0; j < NB_BUFFERS_PER_THREAD; j++)
      th_buf_clear (th_data[i].bufs[j]);
  }
  free (th_data);
  fclose_maybe_compressed (outfile, outfilename);
  return nfreerels_total;
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "renumber", "output file for renumbering table");
  param_list_decl_usage(pl, "out", "output file for free relations");
  for (unsigned int i = 0; i < NB_POLYS_MAX; i++)
  {
    char desc[64], name[8];
    snprintf (desc, 64, "large prime bound on side %u", i);
    snprintf (name, 8, "lpb%u", i);
    param_list_decl_usage(pl, name, desc);
  }
  param_list_decl_usage(pl, "pmin", "do not create freerel below this bound");
  param_list_decl_usage(pl, "pmax", "do not create freerel beyond this bound");
  param_list_decl_usage(pl, "badideals", "file describing bad ideals (for DL)");
  param_list_decl_usage(pl, "addfullcol", "(switch) add a column of 1 in the matrix (for DL)");
  param_list_decl_usage(pl, "t", "number of threads");
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
  char *argv0 = argv[0];
  cado_poly cpoly;
  unsigned long pmin = 2, pmax = 0, nfree;
  renumber_t renumber_table;
  int add_full_col = 0;
  unsigned long lpb[NB_POLYS_MAX] = { 0 };
  unsigned long nb_pthreads = 1;

  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  param_list_configure_switch(pl, "-addfullcol", &add_full_col);

#ifdef HAVE_MINGW
  _fmode = _O_BINARY;     /* Binary open for all files */
#endif

  argv++, argc--;
  if (argc == 0)
    usage (pl, argv0);

  for( ; argc ; )
  {
    if (param_list_update_cmdline(pl, &argc, &argv))
      continue;
    FILE *f;
    if ((f = fopen(argv[0], "r")) != NULL)
    {
      param_list_read_stream(pl, f, 0);
      fclose(f);
      argv++,argc--;
      continue;
    }
    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    usage (pl, argv0);
  }
  /* print command-line arguments */
  verbose_interpret_parameters (pl);
  param_list_print_command_line (stdout, pl);
  fflush (stdout);

  const char * polyfilename = param_list_lookup_string(pl, "poly");
  const char * outfilename = param_list_lookup_string(pl, "out");
  const char * badidealsfilename = param_list_lookup_string(pl, "badideals");
  const char * renumberfilename = param_list_lookup_string(pl, "renumber");
  for (unsigned int i = 0; i < NB_POLYS_MAX; i++)
  {
    char name[8];
    snprintf (name, 8, "lpb%u", i);
    param_list_parse_ulong(pl, name, &lpb[i]);
  }
  param_list_parse_ulong(pl, "pmin", &pmin);
  param_list_parse_ulong(pl, "pmax", &pmax);
  param_list_parse_ulong(pl, "t"   , &nb_pthreads);

  if (polyfilename == NULL)
  {
    fprintf (stderr, "Error, missing -poly command line argument\n");
    usage (pl, argv0);
  }
  if (renumberfilename == NULL)
  {
    fprintf (stderr, "Error, missing -renumber command line argument\n");
    usage (pl, argv0);
  }
  if (outfilename == NULL)
  {
    fprintf (stderr, "Error, missing -out command line argument\n");
    usage (pl, argv0);
  }

  cado_poly_init(cpoly);
  if (!cado_poly_read (cpoly, polyfilename))
  {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

  for (int i = 0; i < cpoly->nb_polys; i++)
  {
    if (lpb[i] == 0)
    {
      fprintf (stderr, "Error, missing -lpb%u command line argument\n", i);
      usage (pl, argv0);
    }
  }

  if (nb_pthreads == 0 || nb_pthreads > 512)
  {
    fprintf (stderr, "Error, the number of threads is incorrect, it must be "
                      "between 1 and 512.\n");
  }

  if (param_list_warn_unused(pl))
    usage (pl, argv0);

  int ratside = cado_poly_get_ratside (cpoly);
  renumber_init_for_writing (renumber_table, cpoly->nb_polys, ratside,
                                                            add_full_col, lpb);
  renumber_write_open (renumber_table, renumberfilename, badidealsfilename,
                       cpoly);

  nfree = allFreeRelations (cpoly, pmin, pmax, lpb, renumber_table,
                            (size_t) nb_pthreads, outfilename);

  /* /!\ Needed by the Python script. /!\ */
  fprintf (stderr, "# Free relations: %lu\n", nfree);
  fprintf (stderr, "Renumbering struct: nprimes=%" PRIu64 "\n",
                   renumber_table->size);

  renumber_write_close (renumber_table, renumberfilename);
  renumber_clear (renumber_table);
  cado_poly_clear (cpoly);
  param_list_clear(pl);

  return 0;
}
