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

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>   /* for _O_BINARY */
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>

#include <gmp.h>
#include "mod_ul.c"
#include "portability.h"
#include "utils.h"
#include "typedefs.h"

// All next define doesn't need to be a power of two... but it's smarter!
#define NB_BUFS (1<<1)             // 2 buffers are sufficient; in fact it works with only one
#define SIZE_BUF_PRIMES (1<<12)    // Size of a bloc of generated primes for one pthread
#define SIZE_BUF_FREE_RELS (1<<10) // Nb of free rels for the primes bloc
#define MIN_BUF_FREE_RELS (1<<8)   // Minimal size in BYTES of free rels buffer before it grows
#define SIZE_BUF_ROOTS (1<<17)     // Nb of roots for the primes bloc; NB: in ASCII!
#define MIN_BUF_ROOTS (1<<10)      // Max size in ASCII of all the roots for ONE p (1024 bytes)

/* Model of this code : One producer -> Many consumers/producers -> one consumer.
   1. One thread produces SIZE_BUF_PRIMES primes by buffer, in NB_BUFS buffers
      (named primes) for each working thread.

      When all the primes are produced, this thread produces one buffer with
      only one "false" prime = (p_r_values_t) (-1). It's the mark of this end.

   2. Each nb_pthread threads load a primes buffers for this thread by a classical
      one produced/many consumer models (2 semaphores).
      Each pthread produces 2 buffers: one of free_relations, each on the form of
      a set of (deg(polynom1)+deg(polynom2)+1 free_rels_buf_t type (unsigned long or
      unsigned int); and one of roots, on the form of a ASCII array.

      When a thread loads a buffer which begins with (p_r_values_t) (-1), the thread
      produces a roots buffer which contains only (p_r_values_t) (-1) and an empty free
      relations buffer, and exit.

   3. With a many producer/one consumer model, the principal programs load the
      roots & free_buffers, writes sequentially the roots, computed the real renumber
      index (sum of all previous index) for the free relations, and write these relations.

      When a roots buffers contains at the beginning (p_r_values_t) (-1), the program
      increments a counter; when this counter == nb_pthreads, the job is done.
*/

// Must be EXACTLY the same than begin/current/end
// in the 2 buffers array of each pthread structure. See below.
typedef struct buf_s {
  void *begin, *current, *end;
} buf_t;

/* Structure for a buffer which keeps the free relations.
   Problem: we cannot compute directly the true renumber_table->size value:
   it's the sum of all previous values in the others pthreads.
   So, the buffer cannot be ASCII, but must stock directly
   a free relation on the form of one type "p" (so, p_r_values_t)
   and [sum degrees of the 2 poly] type renumber_t.
*/
#if __SIZEOF_P_R_VALUES__ == 4 && __SIZEOF_INDEX__ == 4
typedef uint32_t free_rels_buf_t;
#else
typedef uint64_t free_rels_buf_t;
#endif
typedef struct free_rels_s {
  free_rels_buf_t *begin, *current, *end; // Same than buf_t
  index_t renumber;
  unsigned int nb_p, nb_free;             // max SIZE_BUF_PRIMES
} free_rels_t;

typedef struct roots_s {
  char *begin, *current, *end;            // Same than buf_t
} roots_t;

typedef struct primes_s {
  p_r_values_t *begin, *current, *end;
} primes_t;

/* The main structure for the consumer producer pthread pool */
typedef struct producer_s {
  struct {
    sem_t full, empty;  // Classical 2 semaphores consumer/producer model
  } primes_sem;
  struct {
    sem_t full, empty;  // Classical 2 semaphores consumer/producer model
  } roots_sem;          // roots_sem is also free_relations_sem
// Now the read only part
  unsigned int pthread_number, number_of_pthreads;
  unsigned long pmin, pmax, lpb[2], lpbmax;
  size_t deg[2];
  mpz_t *coeff[2];
  mpz_poly_ptr pols[2];
// Now the read-write part
  pthread_t pthread;
  size_t current_buf; // [0, NB_NBUFS[
  free_rels_t free_rels [NB_BUFS];
  roots_t roots [NB_BUFS];
  primes_t primes [NB_BUFS];
} producer_t;

// careful! min_buf is in BYTES, no in nb of the simple type of the buffer
inline void
resize_buf (buf_t *buf, size_t min_buf) {
  if (UNLIKELY(buf->current + min_buf > buf->end)) {
    size_t ind_current = buf->current - buf->begin,
      new_lg = (buf->end - buf->begin) << 1;
    buf->begin = realloc (buf->begin, new_lg);
    if (!buf->begin) {
      perror ("Realloc error\n");
      exit (1);
    }
    buf->current = buf->begin + ind_current;
    buf->end = buf->begin + new_lg;
  }
}

static void *
pthread_primes_producer (void *ptvoid) {
  producer_t *pth = (producer_t *) ptvoid;
  unsigned long lpbmax = pth->lpbmax, p = 2;
  size_t number_of_pthreads = pth->number_of_pthreads, current_buf = 0, i = 0;

  while (p <= lpbmax){
    // here, we are sure we can produce at least one prime
    sem_wait (&(pth[i].primes_sem.empty)); // Need an empty primes buffer
    primes_t my_primes = pth[i].primes[current_buf]; // Careful! it's a local copy!
    my_primes.current = my_primes.begin;

    do { // Main loop to produce primes
      *my_primes.current++ = (p_r_values_t) p;
      p = getprime (p);
    } while (my_primes.current < my_primes.end && p <= lpbmax);

    pth[i].primes[current_buf] = my_primes;
    sem_post (&(pth[i].primes_sem.full)); // Drop the now full primes buffer

    if (UNLIKELY (++i >= number_of_pthreads)) {
      i = 0;
      if (UNLIKELY (++current_buf >= NB_BUFS)) current_buf = 0;
    }

  }
  getprime (0); // Free structures in getprime

  // We have to produce a special end buffer in each pth[].primes
  size_t j = i;
  do {

    sem_wait (&(pth[i].primes_sem.empty));
    primes_t *my_primes = pth[i].primes + current_buf; // Careful! it's a pointer!
    *(my_primes->begin) = (p_r_values_t) (-1);
    my_primes->current = my_primes->begin + 1;
    sem_post (&(pth[i].primes_sem.full));

    if (UNLIKELY (++i >= number_of_pthreads)) {
      i = 0;
      if (UNLIKELY (++current_buf >= NB_BUFS)) current_buf = 0;
    }

  } while (j != i);

  pthread_exit (NULL);
}

static void *
pthread_roots_and_free_rels_producer (void *ptvoid) {
  producer_t *my_pth = ptvoid;
  unsigned long p;
  size_t i, nb_roots[2];
  unsigned long *computed_roots[2];
  computed_roots[0] = malloc_check (my_pth->deg[0] * sizeof(*(computed_roots[0])));
  computed_roots[1] = malloc_check (my_pth->deg[1] * sizeof(*(computed_roots[1])));

  for (;;) {

    sem_wait (&(my_pth->roots_sem.empty)); // Need a empty roots & free relations buffer
    sem_wait (&(my_pth->primes_sem.full)); // Need a produced primes buffer

    roots_t      *my_roots     = my_pth->roots +     my_pth->current_buf;
    free_rels_t  *my_free_rels = my_pth->free_rels + my_pth->current_buf;
    primes_t     *my_primes    = my_pth->primes +    my_pth->current_buf;
    p_r_values_t *my_p         = my_primes->begin;

    // Have we received the special end buffer in primes buffer ?
    if (UNLIKELY(my_primes->current - my_p == 1 && *my_p == (p_r_values_t) -1)) {
      // Yes: we create a special end buffer in roots buffer, and exit
      * (p_r_values_t *) my_roots->begin = (p_r_values_t) -1;
      my_roots->current = (char *) ((void *) my_roots->begin + sizeof (p_r_values_t));
      sem_post (&(my_pth->primes_sem.empty)); // Drop my now consumed primes buffer (not useful)
      sem_post (&(my_pth->roots_sem.full));   // Give my produced roots & free_rels buffer
      // OK, the job is done
      free (computed_roots[0]);
      free (computed_roots[1]);
      pthread_exit (NULL); // The exit of the endless loop and of the pthread
    }

    my_roots->current = my_roots->begin;
    my_free_rels->current = my_free_rels->begin;
    my_free_rels->renumber = 0;
    my_free_rels->nb_p = 0;
    my_free_rels->nb_free = 0;

    // I prefer to separate the 2 cases (rat & alg) and (alg & alg).
    // The goal is to suppress the first case in fews months
    if (my_pth->deg[0] == 1)

      while (my_p < my_primes->current) { // rational & algebraic polynomes, special case
	p = (unsigned long) *my_p++; // p is the current prime
	
	// First, we compute the roots of alg
	if (UNLIKELY(p > my_pth->lpb[1])) {
	  nb_roots[0] = 1;
	  nb_roots[1] = 0;
	} else {
	  nb_roots[0] = (p < my_pth->lpb[0]) ? 1 : 0;
	  nb_roots[1] = mpz_poly_roots_ulong (computed_roots[1], my_pth->pols[1], p);
	  if (UNLIKELY(nb_roots[1] != my_pth->deg[1] &&
		       mpz_divisible_ui_p (my_pth->coeff[1][my_pth->deg[1]], p))) {
	    computed_roots[1][nb_roots[1]++] = computed_roots[1][0];
	    computed_roots[1][0] = p; // p is insered in first place for (the next sort) as a root
	  }
	}

	// Second, we fill my_roots->current buffer by the computed roots
	my_roots->current += renumber_write_p_rat_alg (p, nb_roots[0], computed_roots[1], nb_roots[1], my_roots->current);
	resize_buf ((buf_t *) &(my_roots->begin), MIN_BUF_ROOTS);
	
	// Third, we fill my_free_rels->current buffer by the possible free rels, on the form
	// [1 unsigned long] : p
	// [(degree of alg + 1) unsigned long] : ++(my_free_rels->renumber)
	if (UNLIKELY(nb_roots[0] && nb_roots[1] == my_pth->deg[1] &&
		     p >= my_pth->pmin && p <= my_pth->pmax)) {
	  *(my_free_rels->current)++ = p;
	  for (i = 0; i <= nb_roots[1]; ++i)
	    my_free_rels->current[i] = my_free_rels->renumber + i;
	  my_free_rels->current += i;
	  resize_buf ((buf_t *) &(my_free_rels->begin), MIN_BUF_FREE_RELS);
	  ++(my_free_rels->nb_free);
	}
	my_free_rels->renumber += nb_roots[0] + nb_roots[1];

	++(my_free_rels->nb_p);
      } // Next p in current primes buffer

    else

      while (my_p < my_primes->current) { // 2 algebraics polynomes, normal case
	p = *my_p++; // p is the current prime

	// First, we compute the roots of alg1 and alg2 - same treatment for both
	for (size_t my_alg = 0; my_alg < 2; ++my_alg) {
	  if (LIKELY (p < my_pth->lpb[my_alg])) {
	    nb_roots[my_alg] = mpz_poly_roots_ulong (computed_roots[my_alg], my_pth->pols[my_alg], p);
	    if (UNLIKELY(nb_roots[my_alg] != my_pth->deg[my_alg] &&
			 mpz_divisible_ui_p (my_pth->coeff[my_alg][my_pth->deg[my_alg]], p))) {
	      computed_roots[my_alg][nb_roots[my_alg]++] = computed_roots[my_alg][0];
	      computed_roots[my_alg][0] = p; // p is insered in first place (for the next sort) as a root
	    }
	  }
	  else
	    nb_roots[my_alg] = 0;
	}
	
	// Second, we fill my_roots->current buffer by the computed roots
	my_roots->current += renumber_write_p_2algs (p, computed_roots[0], nb_roots[0], computed_roots[1], nb_roots[1], my_roots->current);
	resize_buf ((buf_t *) &(my_roots->begin), MIN_BUF_ROOTS);

	// Third, we fill my_free_rels->current buffer by the possible free rels, on the form
	// [1 unsigned long] : p
	// [(degree of alg + 1) unsigned long] : ++(my_free_rels->renumber)
	if (UNLIKELY(nb_roots[0] == my_pth->deg[0] && nb_roots[1] == my_pth->deg[1] &&
		     p >= my_pth->pmin && p <= my_pth->pmax)) {
	  *(my_free_rels->current)++ = p;
	  for (i = 0; i < nb_roots[0] + nb_roots[1]; ++i)
	    my_free_rels->current[i] = my_free_rels->renumber + i;
	  my_free_rels->current += i;
	  resize_buf ((buf_t *) &(my_free_rels->begin), MIN_BUF_FREE_RELS);
	  ++(my_free_rels->nb_free);
	}
	my_free_rels->renumber += nb_roots[0] + nb_roots[1];

	++(my_free_rels->nb_p);
      } // Next p in current primes buffer

    sem_post(&(my_pth->primes_sem.empty)); // Drop my now consumed primes buffer
    sem_post(&(my_pth->roots_sem.full));   // Give my produced roots & free_rels buffer

    if (UNLIKELY (++(my_pth->current_buf) >= NB_BUFS)) my_pth->current_buf = 0;

  } // Endless loop
}

/* generate all free relations up to the large prime bound */
/* generate the renumbering table */
static unsigned long MAYBE_UNUSED
allFreeRelations (cado_poly pol, unsigned long pmin, unsigned long pmax,
                  unsigned long lpb[2], renumber_t renumber_table, size_t nb_pthreads,
                  const char *outfilename)
{
  size_t i, j;
  FILE *fpout = fopen_maybe_compressed (outfilename, "w");
  unsigned int sum_degs_add_one = pol->rat->deg + pol->alg->deg + 1;
  unsigned long lpbmax;            // MAX(lpb[0],lpb[1])
  uint64_t total_p = 0;               // Total of the primes
  p_r_values_t total_free_rels = 0;// Total of the free relations
  pthread_t primes_producer_pth;   // pthread for primes producer
  producer_t *pth;                 // primes consumers/root+free_rels producers main array
  size_t current_buf = 0;          // current index buffer(s) for all the pthreads = [0, NB_NBUFS[

  ASSERT_ALWAYS(lpb[0] < sizeof(unsigned long) * CHAR_BIT);
  lpb[0] = 1UL << lpb[0];
  ASSERT_ALWAYS(lpb[1] < sizeof(unsigned long) * CHAR_BIT);
  lpb[1] = 1UL << lpb[1];
  lpbmax = MAX(lpb[0], lpb[1]);

  if (pmax) {
    if (UNLIKELY (pmax > MIN(lpb[0], lpb[1]))) {
      fprintf (stderr, "Error: pmax is greater than MIN(lpb[])\n");
      exit (1);
    }
  }
  else
    pmax = MIN(lpb[0], lpb[1]);

  printf ("Generating freerels for %lu <= p <= %lu\n", pmin, pmax);
  printf ("Generating renumber table for 2 <= p <= %lu\n", lpbmax);
  fflush (stdout);

  stats_data_t stats; /* struct for printing progress */
  /* will print report at 2^10, 2^11, ... 2^23 computed primes and every
   * 2^23 primes after that */
  stats_init (stats, stdout, &total_p, 23, "Looked into", "primes", "", "p");

  /* We generate all free relations from pmin and up to the *minimum*
     of the two large prime bounds, since larger primes will never
     occur on both sides.
     We generate the renumbering table from the first prime (2) and
     up to the *maximum* of the two large prime bounds.
  */

  // First, all the roots & free rels producers, which immediatly wait for primes.
  // Boring but needed initializations.
  pth = malloc_check (nb_pthreads * sizeof (*pth));
  for (i = 0; i < nb_pthreads; ++i) {
    if (sem_init (&(pth[i].roots_sem.full),   0, 0      ) ||
	sem_init (&(pth[i].roots_sem.empty),  0, NB_BUFS) ||
	sem_init (&(pth[i].primes_sem.full),  0, 0      ) ||
	sem_init (&(pth[i].primes_sem.empty), 0, NB_BUFS)) {
      perror ("sem_init error\n");
      exit (1);
    }
    pth[i].current_buf        = 0;
    pth[i].pthread_number     = i;
    pth[i].number_of_pthreads = nb_pthreads;
    pth[i].pmin               = pmin;
    pth[i].pmax               = pmax;
    pth[i].lpb[0]             = lpb[0];            pth[i].lpb[1]   = lpb[1];
    pth[i].lpbmax             = lpbmax;
    pth[i].deg[0]             = pol->rat->deg;     pth[i].deg[1]   = pol->alg->deg;
    pth[i].coeff[0]           = pol->rat->coeff;   pth[i].coeff[1] = pol->alg->coeff;
    pth[i].pols[0]            = pol->pols[0];      pth[i].pols[1]  = pol->pols[1];
    for (j = NB_BUFS; j--;) {
      pth[i].free_rels[j].begin = pth[i].free_rels[j].current =
	malloc_check (SIZE_BUF_FREE_RELS * sum_degs_add_one * sizeof (*(pth[i].free_rels[j].begin)));
      pth[i].roots[j].begin =     pth[i].roots[j].current =
	malloc_check (SIZE_BUF_ROOTS     * sizeof (*(pth[i].roots[j].begin)));
      pth[i].primes[j].begin =    pth[i].primes[j].current =
	malloc_check (SIZE_BUF_PRIMES    * sizeof (*(pth[i].primes[j].begin)));
      pth[i].free_rels[j].end = pth[i].free_rels[j].begin + SIZE_BUF_FREE_RELS;
      pth[i].roots[j].end     = pth[i].roots[j].begin     + SIZE_BUF_ROOTS;
      pth[i].primes[j].end    = pth[i].primes[j].begin    + SIZE_BUF_PRIMES;
    }
    if (pthread_create (&(pth[i].pthread), NULL, pthread_roots_and_free_rels_producer, (void *) (pth + i))) {
      perror ("pthread_roots_and_free_rels_producer pthread creation failed\n");
      exit (1);
    }
  }
  // Second, run the primes producer
  if (pthread_create (&(primes_producer_pth), NULL, pthread_primes_producer, (void *) pth)) {
    perror ("pthread_primes_producer creation failed\n");
    exit (1);
  }

  // Load the roots & free rels buffers, and print them
  for (i = 0;;) {
    sem_wait (&(pth[i].roots_sem.full)); // Need a full root & free rels buffer

    roots_t *roots = pth[i].roots + current_buf;

    // Is it the special end buffer (with only the end mark ?)
    if ((void *) roots->current - (void *) roots->begin == sizeof (p_r_values_t) &&
	*(p_r_values_t *) roots->begin == (p_r_values_t) (-1)) {
      /* OK, all is done: when a pthread which produces roots & free rels gives this
	 special end buffer, all the next pthreads have finished and give also another
	 special end buffer. So it's not useful to read all; one is sufficient and
	 we stop immediatly. */
      break;
    }

    // We write the roots in one ascii bloc
    fwrite (roots->begin, (void *) roots->current - (void *) roots->begin, 1, renumber_table->file);

    // We have to recomputed the real index of the renumber table for the free rels
    free_rels_t free_rels = pth[i].free_rels [current_buf]; // Careful: local copy!
    ASSERT (!((free_rels.current - free_rels.begin) % sum_degs_add_one));
    for (free_rels_buf_t *pt = free_rels.begin; pt < free_rels.current; pt += sum_degs_add_one) {
      fprintf (fpout, "%" PRIx64 ",0:%" PRIx64, (uint64_t) pt[0], (uint64_t) pt[1] + renumber_table->size);
      for (j = 2; j < sum_degs_add_one; ++j)
	fprintf (fpout, ",%" PRIx64, (uint64_t) pt[j] + renumber_table->size);
      fputc ('\n', fpout);
    }
    renumber_table->size += free_rels.renumber;
    total_free_rels += free_rels.nb_free;
    total_p += free_rels.nb_p;

    sem_post (&(pth[i].roots_sem.empty)); // Drop the consumed root & free rels buffer

    if (stats_test_progress(stats)) stats_print_progress (stats, total_p, 0, 0, 0);

    if (UNLIKELY (++i >= nb_pthreads)) {
      i = 0;
      if (UNLIKELY (++current_buf >= NB_BUFS)) current_buf = 0;
    }
  }

  // All is done!
  stats_print_progress (stats, total_p, 0, 0, 1);

  if (pthread_join (primes_producer_pth, NULL)) {
    perror ("Error on the return of the primes producer pthread\n");
    exit (1);
  }
  for (i = nb_pthreads; i--;) {
    if (pthread_join (pth[i].pthread, NULL)) {
      perror ("Error on the return of one of the roots & free relations producer pthread\n");
      exit (1);
    }
    sem_destroy (&pth[i].roots_sem.full);
    sem_destroy (&pth[i].roots_sem.empty);
    sem_destroy (&pth[i].primes_sem.full);
    sem_destroy (&pth[i].primes_sem.empty);
    for (j = NB_BUFS; j--;) {
      free(pth[i].free_rels[j].begin);
      free(pth[i].roots[j].begin);
      free(pth[i].primes[j].begin);
    }
  }
  free (pth);

  fclose_maybe_compressed (fpout, outfilename);

  return total_free_rels;
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "renumber", "output file for renumbering table");
  param_list_decl_usage(pl, "out", "output file for free relations");
  param_list_decl_usage(pl, "lpbr", "rational large prime bound");
  param_list_decl_usage(pl, "lpba", "algebraic large prime bound");
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
    const char *renumberfilename = NULL;
    const char *badidealsfilename = NULL;
    const char *polyfilename = NULL;
    const char *outfilename = NULL;
    char *argv0 = argv[0];
    cado_poly cpoly;
    unsigned long pmin = 2, pmax = 0, nfree;
    renumber_t renumber_table;
    int add_full_col = 0;
    unsigned long lpb[2] = {0, 0};
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

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv))
            continue;
        FILE *f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f);
            fclose(f);
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage (pl, argv0);
    }
    /* print command-line arguments */
    param_list_print_command_line (stdout, pl);
    printf ("\n");
    fflush(stdout);

    polyfilename = param_list_lookup_string(pl, "poly");
    outfilename = param_list_lookup_string(pl, "out");
    badidealsfilename = param_list_lookup_string(pl, "badideals");
    renumberfilename = param_list_lookup_string(pl, "renumber");


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

    param_list_parse_ulong(pl, "lpbr", &lpb[RATIONAL_SIDE]);
    param_list_parse_ulong(pl, "lpba", &lpb[ALGEBRAIC_SIDE]);
    param_list_parse_ulong(pl, "pmin", &pmin);
    param_list_parse_ulong(pl, "pmax", &pmax);
    param_list_parse_ulong(pl, "t"   , &nb_pthreads);

    if (lpb[0] == 0 || lpb[1] == 0)
    {
      fprintf (stderr, "Error, missing -lpbr or -lpba command line argument\n");
      usage (pl, argv0);
    }

    if (param_list_warn_unused(pl))
    {
      usage (pl, argv0);
    }

    int ratside = cado_poly_get_ratside (cpoly);
    renumber_init_for_writing (renumber_table, ratside, add_full_col, lpb);
    renumber_write_open (renumber_table, renumberfilename, badidealsfilename,
                         cpoly);

    nfree = allFreeRelations (cpoly, pmin, pmax, lpb, renumber_table, (size_t) nb_pthreads,
                              outfilename);

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
