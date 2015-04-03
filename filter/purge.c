/* purge --- remove singletons
 * 
 * Copyright 2008, 2009, 2010, 2011, 2012 Alain Filbois, Francois Morain,
 *                                        Paul Zimmermann
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
 */

/*
 * This program works in two passes over the relation files:
 * - the first pass loads in memory only rational primes >= minpr and algebraic
 *   ideals >= minpa, but stores all ideals in the hash table, and keeps a count
 *   of the number of non-stored ideals for each relation.
 *   By default minpr and minpa are taken as rlim and alim respectively.
 *   Then simultaneously singleton removal is performed, and heavy relations
 *   are discarded, until the final excess is 'keep'.
 * - the second pass goes through the relations again, and dumps the remaining
 *   ones in the format needed by 'merge'.

 * This program uses the following data structures:
 * rel_used[i]    - non-zero iff relation i is kept (so far)
 * rel_compact[i] - list of 'h' indices in H table of considered (p,r) for row i
 *                  (terminated by a -1 sentinel)
 * ideals_weight [h] - number of occurrences of h in current relations

 * Exit value:
 * - 0 if enough relations
 * - 1 if an error occurred (then we should abort the factorization)
 * - 2 if not enough relations
 */

#include "cado.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>		/* for _O_BINARY */
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>
#include <pthread.h>
#ifdef HAVE_LIBGEN_H
#include <libgen.h>
#endif

#include "portability.h"

#include "filter_common.h"

/* Some define, no so interesting. */
#define TO_EXPLORE_MIN_BUFFER 16 // Minimal size of the "stack" to explore cliques
#define EXPLORED_MIN_BUFFER 16   // Id, for the relations already explored */
#define RELS_BLOCK (BV_BITS<<4)    // The size of the relations minimal block : BV_BITS << x

//#define STAT
//#define STAT_VALUES_COEFF //STAT must be defined. Interesting only DL

//#define USE_CAVALLAR_WEIGHT_FUNCTION

/* Main variables */

/* This one is passed to all functions, so it's morally a global */
struct purge_data_s {
    /* for minimal changes, don't put these here _yet_ */
    // uint64_t min_index;
    // index_t **rel_compact;	/* see main documentation */
    // weight_t *ideals_weight;
    info_mat_t info;
    /* fd[0]: printed relations */
    /* fd[1]: deleted relations */
    FILE * fd[2];
};
typedef struct purge_data_s purge_data[1];
typedef struct purge_data_s * purge_data_ptr;
typedef const struct purge_data_s * purge_data_srcptr;

/* A clique is a connected components of the relation R, where R(i1,i2) iff
 * i1 and i2 share a prime of weight 2. */
typedef struct {
  float w;   /* Weight of the clique */
  index_t i; /* smallest relation of the clique (index in rel_compact) */
} comp_t;

uint64_t min_index;
index_t **rel_compact;	/* see main documentation */
weight_t *ideals_weight;

static bit_vector rel_used;
static index_t *sum2_index = NULL;	/*sum of rows index for primes of weight 2 */

static uint64_t nrelmax = 0, nprimemax = 0;
static int64_t keep = DEFAULT_FILTER_EXCESS; /* maximun final excess */
static unsigned int npass = DEFAULT_PURGE_NPASS;
static double required_excess = DEFAULT_PURGE_REQUIRED_EXCESS;
static unsigned int npt = DEFAULT_PURGE_NPT;

#ifdef STAT
uint64_t __stat_weight;
#ifdef STAT_VALUES_COEFF
#define STAT_VALUES_COEFF_LEN 10
uint64_t __stat_nbcoeffofvalue[STAT_VALUES_COEFF_LEN + 1];
#endif
#endif

typedef struct index_buffer_s { // Classical buffer (here, a stack in fact)
  index_t *begin, *current, *end;
} index_buffer_t;

/* The main structure for the working pthreads pool which
   compute the heaviest cliques */
typedef struct pth_s {
// Read only part
  unsigned int pthread_number;
  size_t chunk;
// Read-write part
  pthread_t pthread;
  comp_t *clique_graph;     // Array/binary tree sorted of the heaviest cliques
  size_t size_clique_graph; // max & optimal size of the tree = chunk
  comp_t clique;            // current clique
  index_buffer_t to_explore, explored; // stacks to explore and already explored
} pth_t;

/*****************************************************************************/


/* Delete a relation: set rel_used[i] to 0, update the count of primes
 * in that relation.
 * Warning: we only update the count of primes that we consider, i.e.,
 * primes with index >= min_index.
 * CAREFUL: no multithread compatible with " !(--(*o))) " and
 * bit_vector_clearbit.
 */
static unsigned int delete_relation (uint64_t i)
{
    index_t *tab;
    unsigned int nremoveprimes = 0;
    weight_t *o;

    for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++) {
	o = &(ideals_weight[*tab]);
	ASSERT(*o);
	if (*o < UMAX(*o) && !(--(*o)))
	    nremoveprimes++;
    }
    /* We do not free rel_compact[i] as it is freed with my_malloc_free_all */
    bit_vector_clearbit(rel_used, (size_t) i);
    return nremoveprimes;
}

/*****************************************************************************/
/* Code for clique removal.
 * A clique is a connected components of the relation R, where R(i1,i2) iff
 * i1 and i2 share a prime of weight 2.
 * We remove the heaviest cliques.
 * Each ideal h contributes to the weight of the cliques depending on its
 * weight (see weight_function_clique).
 */

static inline
float weight_function_clique(weight_t w)
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

// This function grows a possible too small index_t buffer.
static inline
void resize_buf_index (index_buffer_t *buf, size_t min_buf) {
  if (UNLIKELY(buf->current + min_buf >= buf->end)) {
    size_t ind_current = buf->current - buf->begin,
      new_lg = (buf->end - buf->begin) << 1;
    buf->begin = realloc (buf->begin, new_lg * sizeof (index_t));
    if (!buf->begin) {
      perror ("Realloc error\n");
      exit (1);
    }
    buf->current = buf->begin + ind_current;
    buf->end = buf->begin + new_lg;
  }
}

/* This function insert a clique in a binary tree on the classical array form:
   2 sons at the index 2n+1 and 2n+2 for a father at the index n.
   The property is: the weight of the clique father is less than its 2 sons.
*/
static void insert_clique_pth (pth_t *pth) {
  size_t son, father, place;
  son = pth->size_clique_graph;
  while (son) {
    father = (son - 1) >> 1;
    if (UNLIKELY (pth->clique.w >= pth->clique_graph[father].w)) break;
    son = father;
  }
  place = son;
  son = (pth->size_clique_graph)++;
  while (place != son) {
    father = (son - 1) >> 1;
    pth->clique_graph[son] = pth->clique_graph[father];
    son = father;
  }
  pth->clique_graph[son] = pth->clique;
}
  
/* This function takes a clique where the clique weight is greater than the weight
   of the clique root of a binary tree build by insert_clique_pth.
   The function searchs the right place in this tree for the clique in respect with
   the property described in insert_clique_pth.
*/
static void replace_clique_pth (pth_t *pth) {
  ASSERT_EXPENSIVE (pth->clique_graph->w < pth->clique.w);
  size_t father = 0;
  for (;;) { // In this loop, always, clique_graph[father].w < weight_clique
    float weight1, weight2;
    size_t son = father * 2 + 1; // son is son1
    if (UNLIKELY (son >= pth->chunk)) break; // No son: father is a leaf. Found!
    weight1 = pth->clique_graph[son].w;
    weight2 = pth->clique_graph[son + 1].w;
    if (UNLIKELY(weight1 > weight2)) {
      ++son; // son is here the son 2
      weight1 = weight2;
    }
    // here, weight1 is the lightest weight and son is its clique index
    if (UNLIKELY (weight1 >= pth->clique.w)) break; // Found!
    // The lightest son has a weight lighter than clique, so
    // we have to go down the tree, but before the son replaces its father.
    pth->clique_graph[father] = pth->clique_graph[son];
    father = son;
  }
  pth->clique_graph[father] = pth->clique;
}

/* Compute connected component of row clique->i for the relation R(i1,i2)
 * if rows i1 and i2 share a prime of weight 2.
 * The total weight of the clique is written in pth->clique.w
 *
 * NB: Multithread version! if it exists a row i1 with i1 < clique_pth->clique.i in
 * the connected component, this component has been already found when we
 * have treated the row i1.
 * In this case, the function returns 0.
 * If no, the function returns the number of connected relations.
 */
static int compute_connected_component_pth (pth_t *pth) {
  index_t *primes_of_current_row, current_prime, current_row, the_other_row, nb_rels = 1;
  weight_t current_ideal_weight;

  pth->to_explore.current = pth->to_explore.begin; // Empty buffer
  pth->explored.current   = pth->explored.begin;   // Empty buffer
  current_row = pth->clique.i;
  pth->clique.w = 0.;
  for (;;) { // Loop on all connected rows
    primes_of_current_row = rel_compact[current_row];
    for (;;) { // Loop on all primes of the current row
      next_row:
      current_prime = *primes_of_current_row++;
      if (UNLIKELY (current_prime == UMAX(current_prime))) break; // exit of the current loop
      current_ideal_weight = ideals_weight[current_prime];
      pth->clique.w += weight_function_clique (current_ideal_weight);
      if (LIKELY(current_ideal_weight == 2)) {
	the_other_row = sum2_index[current_prime] - current_row;
	// First, if the_other_row < original row, it's an already found clique. Bye.
	if (UNLIKELY (the_other_row < pth->clique.i)) return 0;
	// Second, is the_other_row already explored ?
	for (index_t *pt = pth->explored.begin; pt < pth->explored.current; pt++)
	  if (UNLIKELY(*pt == the_other_row)) goto next_row;
	// Third, is the_other_row is already in the to_explore buffer ?
	for (index_t *pt = pth->to_explore.begin; pt < pth->to_explore.current; pt++)
	  if (UNLIKELY(*pt == the_other_row)) goto next_row;
	// No: store the new clique row in order to explore it later
	++nb_rels;
	resize_buf_index (&(pth->to_explore), 0);
	*(pth->to_explore.current)++ = the_other_row;
      }
    }
    // current row is now explored
    resize_buf_index (&(pth->explored), 0);
    *(pth->explored.current)++ = current_row;
    // We need another row to explore, or it's the end
    if (UNLIKELY(pth->to_explore.current == pth->to_explore.begin)) break;
    current_row = *(--(pth->to_explore.current));
  }
  return nb_rels;
}

/* Delete connected component of row current_row
 * CAREFUL: this code itself is multithread compatible, but
 * it calls delete_relation, which is NOT compatible!
 * Warning: we might have some H->hashcount[h] = 3, which is decreased
 * to 2, but we don't want to treat that case. Thus we check in addition
 * that sum2_index[h] <> 0, which only occurs when H->hashcount[h] = 2
 * initially.
 */
static index_t delete_connected_component_nopth (pth_t *pth, index_t current_row, uint64_t *nprimes)
{
  index_t *primes_of_current_row, current_prime, the_other_row;
  weight_t current_ideal_weight;

  pth->to_explore.current = pth->to_explore.begin;   // Empty buffer
  pth->explored.current   = pth->explored.begin;     // Empty buffer
  for (;;) { // Loop on all connected rows
    primes_of_current_row = rel_compact[current_row];
    for (;;) { // Loop on all primes of the current row
      next_row:
      current_prime = *primes_of_current_row++;
      if (UNLIKELY (current_prime == UMAX(current_prime))) break; // exit of the current loop
      current_ideal_weight = ideals_weight[current_prime];
      if (LIKELY(current_ideal_weight == 2 && sum2_index[current_prime])) {
	the_other_row = sum2_index[current_prime] - current_row;
	// is the_other_row already explored ?
	for (index_t *pt = pth->explored.begin; pt < pth->explored.current; pt++)
	  if (UNLIKELY(*pt == the_other_row)) goto next_row;
	// Is the_other_row already in the to explore buffer ?
	for (index_t *pt = pth->to_explore.begin; pt < pth->to_explore.current; pt++)
	  if (UNLIKELY(*pt == the_other_row)) goto next_row;
	// No: store the new clique row in order to explore it later
	resize_buf_index (&(pth->to_explore), 0);
	*(pth->to_explore.current)++ = the_other_row;
      }
    }
    // current row is now explored
    resize_buf_index (&(pth->explored), 0);
    *(pth->explored.current)++ = current_row;
    // We need another row to explore, or it's the end
    if (pth->to_explore.current == pth->to_explore.begin) break;
    current_row = *(--(pth->to_explore.current));
    }
  // Now, we deleted all rows explored
  for (index_t *pt = pth->explored.begin; pt < pth->explored.current; pt++)
    *nprimes -= delete_relation (*pt);
  return (pth->explored.current - pth->explored.begin);
}

/* This MT function search the heaviest pth->chunk cliques in
   interlaced parts of rel_used cliques.
*/
static void *search_chunk_max_cliques (void *pt) {
  pth_t *pth = (pth_t *) pt;
  index_t end_step_rels;
  bv_t bv, *pbv;

  // Init of the structure & malloc.
  // If chunk is even, an artificial clique must exists at [chunk]
  // to avoid a node-father with only one node-son
  if (!(pth->chunk & 1)) { // Artificial clique to avoid a father with only one son
    pth->clique_graph = (comp_t *) malloc_check (sizeof(comp_t) * (pth->chunk + 1));
    pth->clique_graph[pth->chunk].i = UMAX(pth->clique_graph[pth->chunk].i);
    pth->clique_graph[pth->chunk].w = 1E18; // Avoid MAXFLOAT for overflow
  }
  else
    pth->clique_graph = (comp_t *) malloc_check (sizeof(comp_t) * pth->chunk);
  pth->size_clique_graph = 0;
  pth->to_explore.begin = pth->to_explore.current
    = (index_t *) malloc_check (sizeof(index_t) * TO_EXPLORE_MIN_BUFFER);
  pth->to_explore.end = pth->to_explore.begin + TO_EXPLORE_MIN_BUFFER;
  pth->explored.begin = pth->explored.current
    = (index_t *) malloc_check (sizeof(index_t) * EXPLORED_MIN_BUFFER);
  pth->explored.end = pth->explored.begin + EXPLORED_MIN_BUFFER;

  // Now the first begin of the search
  pth->clique.i = pth->pthread_number * RELS_BLOCK;
  if (UNLIKELY (pth->clique.i >= nrelmax)) pthread_exit (NULL);
  // And the first end
  end_step_rels = pth->clique.i + RELS_BLOCK;

  pbv = rel_used->p + (pth->clique.i >> LN2_BV_BITS);
  while (pth->clique.i < nrelmax) {
    index_t save_i = pth->clique.i;
    for (bv = *pbv++; bv; ++(pth->clique.i), bv >>= 1) {
      if (LIKELY(bv & 1)) {
	unsigned int nb_rels_connected = compute_connected_component_pth (pth);
	if (UNLIKELY (!nb_rels_connected)) continue;
	if (UNLIKELY (pth->size_clique_graph < pth->chunk))
	  insert_clique_pth (pth);
	else 
	  if (UNLIKELY(pth->clique.w > pth->clique_graph->w))
	    replace_clique_pth (pth);
      }
    }
    pth->clique.i = save_i + BV_BITS;
    // Have we reach the actuel end ?
    if (UNLIKELY (pth->clique.i == end_step_rels)) {
      // The new begin & end.
      end_step_rels += npt * RELS_BLOCK;
      pth->clique.i = end_step_rels - RELS_BLOCK;
      pbv = rel_used->p + (pth->clique.i >> LN2_BV_BITS);
    }
  }
  pthread_exit (NULL);
  return NULL;
}

#ifdef HAVE_SYNC_FETCH
/* The main structure for the working pthreads pool which
   compute sum2_index, only if HAVE_SYNC_FETCH is defined.
*/
typedef struct sum2_pth_s {
  // Read only part
  unsigned int pthread_number;
  // Read-write part
  pthread_t pthread;
  index_t rels_found;
} sum2_pth_t;

MAYBE_UNUSED  
static void *compute_sum2_index (void *pt) {
  sum2_pth_t *sum2_pth = (sum2_pth_t *) pt;
  index_t j = nrelmax / npt;
  index_t i = (j * sum2_pth->pthread_number) & ~((size_t) (BV_BITS - 1));
  index_t end = (sum2_pth->pthread_number == npt - 1) ? nrelmax :
    (j * (sum2_pth->pthread_number + 1)) & ~((size_t) (BV_BITS - 1));
  bv_t bv, *pbv;
  index_t rels_found = 0, h, *myrelcompact;

  pbv = rel_used->p + (i >> LN2_BV_BITS);
  for (; i < end; i += BV_BITS) {
    for (j = i, bv = *pbv++; bv; ++j, bv >>= 1) {
      if (LIKELY(bv & 1)) {
	++rels_found;
	myrelcompact = rel_compact[j];
	for (;;) {
	  h = *myrelcompact++;
	  if (UNLIKELY (h == UMAX(h))) break;
	  if (LIKELY (ideals_weight[h] == 2))
	    __sync_add_and_fetch (sum2_index + h, j);
	}
      }
    }
  }
  sum2_pth->rels_found = rels_found;
  pthread_exit (NULL);
  return NULL;
}
#endif

static void
cliques_removal(int64_t target_excess, uint64_t * nrels, uint64_t * nprimes)
{
  int64_t excess = (((int64_t) *nrels) - *nprimes);
  size_t i, chunk;
  index_t N = 0;		/* number of rows */

  ASSERT(npt);

  if (excess <= keep || excess <= target_excess) return;
  chunk = (size_t) (excess - target_excess);

  /* first collect sums for primes with weight 2, and compute total weight.
     If HAVE_SYNC_FETCH is defined, this part is done with the previous
     multithread function; if not, it's done in sequential, immediatly.
  */
  memset(sum2_index, 0, nprimemax * sizeof(index_t));
#ifndef HAVE_SYNC_FETCH
  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, (size_t) i)) {
      index_t h, *myrelcompact;
      ++N;
      for (myrelcompact = rel_compact[i];
	   (h = *myrelcompact++) != UMAX(h);)
	if (ideals_weight[h] == 2)
	  sum2_index[h] += i;
    }
#else
  sum2_pth_t *sum2_pth = malloc_check (npt * sizeof (*sum2_pth));
  for (i = npt; i--; ) {
    sum2_pth[i].pthread_number = i;
    if (pthread_create (&(sum2_pth[i].pthread), NULL,
			compute_sum2_index, (void *) (sum2_pth + i))) {
      perror ("compute_sum2_index pthread creation failed\n");
      exit (1);
    }
  }
  for (i = npt; i--; ) {
    pthread_join (sum2_pth[i].pthread, NULL);
    N += sum2_pth[i].rels_found;
  }
  free (sum2_pth);
#endif

  ASSERT_ALWAYS(N == *nrels);
  
  // Second we search in parallel the chunk heaviest cliques.
  // For this, each thread search its chunk heaviest cliques.
  pth_t *pth = malloc_check (npt * sizeof (*pth));
  memset (pth, 0, npt * sizeof (*pth));
  for (i = npt; i--; ) {
    pth[i].pthread_number = i;
    pth[i].chunk = chunk;
    if (pthread_create (&(pth[i].pthread), NULL,
			search_chunk_max_cliques, (void *) (pth + i))) {
      perror ("search_chunk_max_cliques pthread creation failed\n");
      exit (1);
    }
  }
  for (i = npt; i--; ) pthread_join (pth[i].pthread, NULL);

  
  /* At this point, in each pth[i].graph_cliques we have
     chunk cliques at most.
     Boring problem: we need if possible chunk cliques in pth->clique_graph
     (in pth[0].clique_graph). If not enough, we try to add the
     others pth[i].clique_graph, until chunk.
  */
  if (UNLIKELY (pth->size_clique_graph < chunk)) {
    for (i = 1; i < npt; ++i) {
      if (UNLIKELY(pth[i].size_clique_graph + pth->size_clique_graph >= chunk)) {
	while (pth->size_clique_graph < chunk) {
	  pth->clique = pth[i].clique_graph[--pth[i].size_clique_graph];
	  insert_clique_pth (pth);
	}
	// Right, we have chunk cliques in pth[0].clique_graph now
	break;
      }
      // We could insert all pth[i].clique_graph in pth[0].clique_graph
      while (pth[i].size_clique_graph) {
	pth->clique = pth[i].clique_graph[--pth[i].size_clique_graph];
	insert_clique_pth (pth);
      }
    }
  }
  // we do a fusion with pth->clique_graph and all pth[i != 0].clique_graph
  if (LIKELY(pth->size_clique_graph == chunk)) // if NOT, pth[1...npt[.clique_graph are empty
    for (i = 1; i < npt; ++i)
      while (pth[i].size_clique_graph)
	if (UNLIKELY(pth[i].clique_graph[--pth[i].size_clique_graph].w > pth->clique_graph->w)) {
	  pth->clique = pth[i].clique_graph[pth[i].size_clique_graph];
	  replace_clique_pth (pth);
	}
  // We could suppress pth[1...npt[
  for (i = 1; i < npt; ++i) {
    free (pth[i].clique_graph);
    free (pth[i].to_explore.begin);
    free (pth[i].explored.begin);
  }
  // pth->size_clique_graph has at most the chunk heaviest cliques -> suppress them!
  comp_t *pt = pth->clique_graph + pth->size_clique_graph;
  while (pt > pth->clique_graph)
    *nrels -= delete_connected_component_nopth (pth, (--pt)->i, nprimes);
  
  fprintf(stdout, "    deleted %zu heavier connected components at %2.2lf\n",
	  pth->size_clique_graph, seconds());
#if DEBUG >= 1
  fprintf(stdout, "    DEBUG: nb heaviest cliques=%zu chunk=%u target=%u\n",
	  pth->size_clique_graph, chunk, target_excess);
#endif
  // We could suppress pth[0] and pth itself
  free (pth->clique_graph);
  free (pth->to_explore.begin);
  free (pth->explored.begin);
  free (pth);
}

/*****************************************************************************/
/* Code for singletons removal.
 * Exist in multithread if __sync_sub_and_fetch exists.
 */

#ifndef HAVE_SYNC_FETCH
static void onepass_singleton_removal(uint64_t * nrels, uint64_t * nprimes)
{
    index_t *tab;
    uint64_t i, nremoverels = 0, nremoveprimes = 0;

    for (i = 0; i < nrelmax; i++) {
	if (bit_vector_getbit(rel_used, (size_t) i)) {
	    for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++) {
		if (ideals_weight[*tab] == 1) {
		    nremoveprimes += delete_relation(i);
		    nremoverels++;
		    break;
		}
	    }
	}
    }
    *nrels -= nremoverels;
    *nprimes -= nremoveprimes;
}

#else				/* ifndef HAVE_SYNC_FETCH */

typedef struct {
    unsigned int nb;
    pthread_t mt;
    uint64_t begin, end, sup_nrel, sup_npri;
} ti_t;
static ti_t *ti;

/* Hightest criticality for performance. I inline all myself. */
static void onepass_thread_singleton_removal(ti_t * mti)
{
  index_t *tab;
  uint64_t i;
  weight_t *o;
  bv_t j;
  
  mti->sup_nrel = mti->sup_npri = 0;
  for (i = mti->begin; i < mti->end; i++) {
    j = (((bv_t) 1) << (i & (BV_BITS - 1)));
    
    if (rel_used->p[i >> LN2_BV_BITS] & j)
      for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++)
	if (UNLIKELY(ideals_weight[*tab] == 1)) {
	  for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++) {
	    o = &(ideals_weight[*tab]);
	    ASSERT(*o);
	    if (*o < UMAX(*o) && !__sync_sub_and_fetch(o, 1))
	      (mti->sup_npri)++;
	  }
	  /* rel_compact[i] is not freed , it is freed with my_malloc_free_all */
	  rel_used->p[i >> LN2_BV_BITS] &= ~j;
	  (mti->sup_nrel)++;
	  break;
	}
  }
  pthread_exit(NULL);
}

static void
onepass_singleton_parallel_removal(unsigned int nb_thread, uint64_t * nrels,
				   uint64_t * nprimes)
{
    pthread_attr_t attr;
    uint64_t pas, incr;
    unsigned int i;
    int err;

    ti = (ti_t *) malloc(nb_thread * sizeof(ti_t));
    ASSERT_ALWAYS(ti != NULL);
    ti[0].begin = 0;
    pas = (nrelmax / nb_thread) & ((uint64_t) ~ (BV_BITS - 1));
    incr = 0;
    for (i = 0, incr = 0; i < nb_thread - 1;) {
	incr += pas;
	ti[i].nb = i;
	ti[i++].end = incr;
	ti[i].begin = incr;
    }
    ti[i].nb = i;
    ti[i].end = nrelmax;
    pthread_attr_init(&attr);
    pthread_attr_setstacksize(&attr, 1 << 16);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    for (i = 0; i < nb_thread; i++)
	if ((err =
	     pthread_create(&ti[i].mt, &attr,
			    (void *) onepass_thread_singleton_removal,
			    &ti[i]))) {
	    fprintf(stderr,
		    "onepass_singleton_parallel_removal : pthread_create error 1: %d. %s\n",
		    err, strerror(errno));
	    exit(1);
	}
    for (i = 0; i < nb_thread; i++) {
	pthread_join(ti[i].mt, NULL);
	*nrels -= ti[i].sup_nrel;
	*nprimes -= ti[i].sup_npri;
    }
    pthread_attr_destroy(&attr);
    if (ti != NULL)
	free(ti);
    ti = NULL;
}
#endif				/* ifdef HAVE_SYNC_FETCH */

static void
remove_all_singletons(uint64_t * nrels, uint64_t * nprimes, int64_t * excess)
{
    uint64_t oldnrels;
    *excess = (((int64_t) * nrels) - *nprimes);
    fprintf(stdout,
	    "  nrels=%" PRIu64 " nprimes=%" PRIu64 " excess=%" PRId64 "\n",
	    *nrels, *nprimes, *excess);
    do {
	oldnrels = *nrels;
#ifdef HAVE_SYNC_FETCH
	ASSERT(npt);
	onepass_singleton_parallel_removal(npt, nrels, nprimes);
#else
	onepass_singleton_removal(nrels, nprimes);
#endif
	*excess = (((int64_t) * nrels) - *nprimes);
	fprintf(stdout,
		"  new_nrels=%" PRIu64 " new_nprimes=%" PRIu64 " excess=%" PRId64
		"" " at %2.2lf\n", *nrels, *nprimes, *excess, seconds());
    } while (oldnrels != *nrels);
}

static void singletons_and_cliques_removal(uint64_t * nrels, uint64_t * nprimes)
{
    uint64_t oldnrels = 0;
    int64_t oldexcess, excess, target_excess;
    unsigned int count;

    //First step of singletons removal
    remove_all_singletons(nrels, nprimes, &excess);

    if (excess <= 0) {		/* covers case nrel = nprimes = 0 */
	fprintf(stdout, "number of relations <= number of ideals\n");
	exit(2);
    }

    if ((double) excess < required_excess * ((double) *nprimes)) {
	fprintf(stdout,
		"(excess / nprimes) = %.2f < %.2f. See -required_excess "
		"argument.\n", ((double) excess / (double) *nprimes),
		required_excess);
	exit(2);
    }

    /* adjust npass so that each pass removes at least about 1% wrt the number
       of ideals */
    if ((uint64_t) excess / npass < *nprimes / 100)
      npass = 1 + (100 * excess) / *nprimes;

    int64_t chunk = excess / npass;

    //npass pass of clique removal + singletons removal
    for (count = 0; count < npass && excess > 0; count++) {
	oldnrels = *nrels;
	oldexcess = excess;
	target_excess = excess - chunk;
	if (target_excess < keep)
	    target_excess = keep;
	fprintf(stdout, "Step %u on %u: target excess is %" PRId64 "\n",
		count + 1, npass, target_excess);
	cliques_removal(target_excess, nrels, nprimes);

	remove_all_singletons(nrels, nprimes, &excess);
	fprintf(stdout, "  [each excess row deleted %2.2lf rows]\n",
		(double) (oldnrels - *nrels) / (double) (oldexcess - excess));
    }


    /* May need an extra pass of clique removal + singletons removal if excess is
       still larger than keep. It may happen due to the fact that each clique does
       not make the excess go down by one but can (rarely) left the excess
       unchanged. */
    if (excess > keep) {
	oldnrels = *nrels;
	oldexcess = excess;
	target_excess = excess - chunk;
	target_excess = keep;

	fprintf(stdout, "Step extra: target excess is %" PRId64 "\n",
		target_excess);
	cliques_removal(target_excess, nrels, nprimes);

	remove_all_singletons(nrels, nprimes, &excess);
	fprintf(stdout, "  [each excess row deleted %2.2lf rows]\n",
		(double) (oldnrels - *nrels) / (double) (oldexcess - excess));
    }
}

/*****************************************************************************/
/* I/O functions */

/* Callback functions called by filter_rels */
void *insert_rel_into_table(purge_data_ptr arg, earlyparsed_relation_ptr rel)
{
    ASSERT_ALWAYS(rel->num < nrelmax);

    arg->info.nprimes +=
        insert_rel_in_table_no_e(rel, min_index, 
                rel_compact, ideals_weight);
#ifdef STAT
    /* here we also used to accumulate the number of primes above
     * min_index in arg->info.W */
    arg->info.W += earlyparsed_relation_nb_above_min_index(rel, min_index);
#endif

    return NULL;
}

void *thread_print(purge_data_ptr arg, earlyparsed_relation_ptr rel)
{
    if (bit_vector_getbit(rel_used, rel->num)) {
        arg->info.W += rel->nb;
        fputs(rel->line, arg->fd[0]);
    } else if (arg->fd[1] != NULL) {
        fputs(rel->line, arg->fd[1]);
    }
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
  param_list_decl_usage(pl, "nprimes",
                                  "number of prime ideals in renumber table");
  param_list_decl_usage(pl, "minindex", "index of the first considered prime");
  param_list_decl_usage(pl, "keep", "wanted excess at the end of purge "
                                    "(default " STR(DEFAULT_FILTER_EXCESS) ")");
  param_list_decl_usage(pl, "npass", "maximal number of steps of clique removal "
                                     "(default " STR(DEFAULT_PURGE_NPASS) ")");
  param_list_decl_usage(pl, "required_excess", "\% of excess required at the "
                            "end of the 1st singleton removal step (default "
                            STR(DEFAULT_PURGE_REQUIRED_EXCESS) ")");
  param_list_decl_usage(pl, "outdel", "outfile for deleted relations (for DL)");
  param_list_decl_usage(pl, "npthr", "number of threads (default " STR(DEFAULT_PURGE_NPT) ")");
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
    int k;
    param_list pl;
    min_index = UMAX(uint64_t);
    uint64_t nrels, nprimes;
    size_t tot_alloc_bytes = 0;
    char ** input_files;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;		/* Binary open for all files */
#endif

    double wct0 = wct_seconds();

    param_list_init(pl);
    declare_usage(pl);
    argv++,argc--;

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
    param_list_parse_uint64(pl, "nrels", &nrelmax);
    param_list_parse_uint64(pl, "nprimes", &nprimemax);
    param_list_parse_int64(pl, "keep", &keep);

    /* Only look at relations of index above minindex */
    param_list_parse_uint64(pl, "minindex", &min_index);


    param_list_parse_uint(pl, "npthr", &npt);
    param_list_parse_uint(pl, "npass", &npass);
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
    if (!deletedname) {
        fprintf(stderr, "Error, option -outdel is mandatory\n");
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
    if (nrelmax == 0)
    {
      fprintf(stderr, "Error, missing -nrels command line argument "
                      "(or nrels = 0)\n");
      usage(pl, argv0);
    }
    if (nprimemax == 0)
    {
      fprintf(stderr, "Error, missing -nprimes command line argument "
                      "(or nprimes = 0)\n");
      usage(pl, argv0);
    }
    if (min_index > nprimemax)
    {
      fprintf(stderr, "Error, missing -minindex command line argument "
                      "or (minindex > nprimes)\n");
      usage(pl, argv0);
    }
    /* If nrels or nprimes > 2^32, then we need index_t to be 64-bit */
    if (((nprimemax >> 32) != 0 || (nrelmax >> 32) != 0) && sizeof(index_t) < 8)
    {
      fprintf(stderr, "Error, -nrels or -nprimes is too large for a 32-bit "
                      "program\nSee #define __SIZEOF_INDEX__ in typedefs.h\n");
      exit(EXIT_FAILURE);
    }

    /* Printing relevant information */
    fprintf(stdout, "Weight function used during clique removal:\n"
                    "  0     1     2     3     4     5     6     7\n");
    for (k = 0; k < 8; k++)
      fprintf(stdout, "%0.3f ", weight_function_clique((weight_t) k));
    fprintf(stdout, "\n");

    fprintf(stdout, "Number of relations is %" PRIu64 "\n", nrelmax);
    fprintf(stdout, "Number of prime ideals below the two large prime bounds: "
                    "%" PRIu64 "\n", nprimemax);
    /*}}}*/

    /* {{{ Allocating memory. We are keeping track of the total
     * malloc()'ed amount only for informational purposes */

    /* {{{ Some macros for tracking memory-consuming variables */
#define ALLOC_VERBOSE_MALLOC(type_, variable_, amount_) do {		\
    variable_ = (type_ *) malloc(amount_ * sizeof(type_));		\
    ASSERT_ALWAYS(variable_ != NULL);					\
    size_t cur_alloc = amount_ * sizeof(type_);                         \
    tot_alloc_bytes += cur_alloc;                                       \
    fprintf(stdout,							\
            "Allocated " #variable_ " of %zuMB (total %zuMB so far)\n",	\
	    cur_alloc >> 20, tot_alloc_bytes >> 20);                  	\
} while (0)

#define ALLOC_VERBOSE_CALLOC(type_, variable_, amount_) do {		\
    ALLOC_VERBOSE_MALLOC(type_, variable_, amount_);                    \
    /* Do this now so that we crash early if kernel overcommitted memory\
     */									\
    memset(variable_, 0, amount_);					\
} while (0)

#define ALLOC_VERBOSE_BIT_VECTOR(variable_, amount_) do {		\
    bit_vector_init(variable_, amount_);				\
    size_t cur_alloc = bit_vector_memory_footprint(variable_);		\
    tot_alloc_bytes += cur_alloc;					\
    fprintf(stdout,                                                     \
            "Allocated " #variable_ " of %zuMB (total %zuMB so far)\n",	\
	    cur_alloc >> 20, tot_alloc_bytes >> 20);			\
} while (0)
    /* }}} */

    ALLOC_VERBOSE_MALLOC(index_t, sum2_index, nprimemax);


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
    purge_data pd;

    memset(pd, 0, sizeof(purge_data));

    ALLOC_VERBOSE_MALLOC(index_t*, rel_compact, nrelmax);
    ALLOC_VERBOSE_CALLOC(weight_t, ideals_weight, nprimemax);

    fprintf(stdout, "Pass 1, reading and storing ideals with index h >= "
            "%" PRIu64 "\n", min_index);

    nrels = nrelmax;

    /* first pass over relations in files */
    /* Note: Now that we no longer take a bitmap on input, all
     * relations are considered active at this point, so that we
     * do not need to pass a bitmap to filter_rels */
    pd->info.nrels = filter_rels(
            input_files,
            (filter_rels_callback_t) &insert_rel_into_table, pd,
            EARLYPARSE_NEED_INDEX,
            NULL, NULL);

    if (pd->info.nrels != nrels) {
	fprintf(stderr,
		"Error, -nrels value should match the number of scanned "
		"relations\nexpected %" PRIu64 " relations, found %" PRIu64
		"\n", nrelmax, pd->info.nrels);
        abort();
    }

    nprimes = pd->info.nprimes;

    tot_alloc_bytes += get_my_malloc_bytes();
    fprintf(stdout, "Allocated rel_compact[i] %zuMB (total %zuMB so far)\n",
	    get_my_malloc_bytes() >> 20, tot_alloc_bytes >> 20);
#ifdef STAT
    {
	size_t tmp = ((uint64_t) buf_arg.W + nrelmax) * sizeof(index_t);
	double ratio = 100.0 * (double) (((double) tmp) /
					 ((double) get_my_malloc_bytes()));
	fprintf(stdout,
		"STAT: W_active=%1.0f\nSTAT: Should take %zuMB in memory, "
		"take %zuMB (%.2f %%)\n", buf_arg.W, tmp >> 20,
		get_my_malloc_bytes() >> 20, ratio);
    }
#endif

    ALLOC_VERBOSE_BIT_VECTOR(rel_used, nrels);
    bit_vector_set(rel_used, 1);

    singletons_and_cliques_removal(&nrels, &nprimes);
    if (nrels < nprimes) {
      fprintf(stdout, "number of relations < number of ideals\n");
      exit(2);
    }
    if (nrels == 0 || nprimes == 0) {
      fprintf(stdout, "number of relations or number of ideals is 0\n");
      exit(2);
    }

    /* free rel_compact[i] and rel_compact. We no longer need them */
    tot_alloc_bytes -= get_my_malloc_bytes();
    fprintf(stdout, "Freed rel_compact[i] %zuMB (total %zuMB so far)\n",
            get_my_malloc_bytes() >> 20, tot_alloc_bytes >> 20);
    my_malloc_free_all();

    free(rel_compact);
    size_t tmp = (nrelmax * sizeof(index_t *));
    tot_alloc_bytes -= tmp;
    fprintf(stdout, "Freed rel_compact %zuMB (total %zuMB so far)\n",
            tmp >> 20, tot_alloc_bytes >> 20);

    /* reread the relation files and convert them to the new coding */
    fprintf(stdout, "Storing remaining relations...\n");

    if (!(pd->fd[0] = fopen_maybe_compressed(purgedname, "w"))) {
	fprintf(stderr, "Error, cannot open file %s for writing.\n",
		purgedname);
	exit(1);
    }

    if (deletedname != NULL) {
	if (!(pd->fd[1] = fopen_maybe_compressed(deletedname, "w"))) {
	    fprintf(stderr, "Error, cannot open file %s for writing.\n",
		    deletedname);
	    exit(1);
	}
      else
      {
          fprintf(pd->fd[1], "# %" PRIu64 "\n", nrelmax - nrels);
      }
    }

    /* Write the header line for the file of remaining relations:
     * compute last index i such that ideals_weight[i] != 0
     */
    {
	uint64_t last_used = nprimemax - 1;
	while (ideals_weight[last_used] == 0)
	    last_used--;

	fprintf(pd->fd[0], "# %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", nrels,
		last_used + 1, nprimes);
    }

    /* second pass over relations in files */
    filter_rels(
            input_files,
            (filter_rels_callback_t) &thread_print, pd,
            EARLYPARSE_NEED_LINE,
            NULL, NULL);

    /* write final values to stdout */
    /* This output, incl. "Final values:", is required by the script */
    fprintf(stdout, "Final values:\nnrels=%" PRIu64 " nprimes=%" PRIu64 " "
            "excess=%" PRId64 "\nweight=%1.0f weight*nrels=%1.2e\n",
            nrels, nprimes, ((int64_t) nrels) - nprimes, pd->info.W,
            pd->info.W * (double) nrels);
    fflush(stdout);

    /* Free allocated stuff */
    if (ideals_weight != NULL)
	free(ideals_weight);
    ideals_weight = NULL;
    if (sum2_index != NULL)
	free(sum2_index);
    sum2_index = NULL;

    bit_vector_clear(rel_used);

    if (filelist)
        filelist_clear(input_files);

    fclose_maybe_compressed(pd->fd[0], purgedname);
    if (pd->fd[1]) fclose_maybe_compressed(pd->fd[1], deletedname);

    /* print usage of time and memory */
    print_timing_and_memory(wct0);

    param_list_clear(pl);

    return 0;
}
