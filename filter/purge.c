/* purge --- remove singletons

Copyright 2008, 2009, 2010, 2011, 2012 Alain Filbois, Francois Morain,
                                       Paul Zimmermann

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

/* References:
 * On the Number Field Sieve Integer Factorisation Algorithm,
   Stefania Cavallar, PhD Thesis, University of Leiden, 2002.
*/

/*
  This program works in two passes over the relation files:
  - the first pass loads in memory only rational primes >= minpr and algebraic
    ideals >= minpa, but stores all ideals in the hash table, and keeps a count
    of the number of non-stored ideals for each relation.
    By default minpr and minpa are taken as rlim and alim respectively.
    Then simultaneously singleton removal is performed, and heavy relations
    are discarded, until the final excess is 'keep'.
  - the second pass goes through the relations again, and dumps the remaining
    ones in the format needed by 'merge'.

  This program uses the following data structures:
  rel_used[i]    - non-zero iff relation i is kept (so far)
  rel_compact[i] - list of 'h' indices in H table of considered (p,r) for row i
                   (terminated by a -1 sentinel)
  ideals_weight [h] - number of occurrences of h in current relations

Exit value:
- 0 if enough relations
- 1 if an error occurred (then we should abort the factorization)
- 2 if not enough relations
*/

/*
  index_t : (32 + 32 * need64) bits. Signed.
  p_r_values_t : 32 if H.hm < 2^32-1, otherwise 64 bits.
  Note : sizeof(index_t)>=sizeof(p_r_values_t)
*/

#include "cado.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>
#ifdef HAVE_LIBGEN_H
#include <libgen.h>
#endif

#include "mod_ul.c"
#include "portability.h"
#include "utils.h"
#include "typedefs.h"

#include "filter_utils.h"

//#define STAT
//#define STAT_VALUES_COEFF //STAT must be defined. Interesting only DL

//#define USE_CAVALLAR_WEIGHT_FUNCTION

#define DEFAULT_NPASS 50
#define DEFAULT_KEEP 160
#define DEFAULT_REQUIRED_EXCESS 0.1
#define DEFAULT_NPT 4

/* Main variables */
char *argv0;                    /* = argv[0]; */

static index_t **rel_compact  = NULL; /* see above */
weight_t *ideals_weight = NULL;

static char ** fic;
static bit_vector rel_used, Tbv;
static index_t *sum2_index = NULL; /*sum of rows index for primes of weight 2*/
static index_t relsup, prisup;

static uint64_t nrelmax = 0,
                nprimemax = 0;
static int64_t keep = DEFAULT_KEEP; /* maximun final excess */
static unsigned int npass = DEFAULT_NPASS;
static double required_excess = DEFAULT_REQUIRED_EXCESS;
static unsigned int npt = DEFAULT_NPT;

static float w_ccc;
static uint8_t binfilerel; /* True (1) if a rel_used file must be read */
static uint8_t boutfilerel;/* True (1) if a rel_used file must be written */

#ifdef STAT
uint64_t __stat_weight;
#ifdef STAT_VALUES_COEFF
#define STAT_VALUES_COEFF_LEN 10
uint64_t __stat_nbcoeffofvalue[STAT_VALUES_COEFF_LEN+1];
#endif
#endif



static const unsigned char nbbits[256] = {
  0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8 };


/*****************************************************************************/


/* Delete a relation: set rel_used[i] to 0, update the count of primes
   in that relation.
   Warning: we only update the count of primes that we consider, i.e.,
   primes with index >= min_index.
*/
static index_t
delete_relation (index_t i)
{
  index_t *tab;
  index_t nremoveprimes = 0;
  HC_T *o;

  for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++)
  {
    o = &(ideals_weight[*tab]);
    ASSERT(*o);
    if (*o<UMAX(*o) && !(--(*o))) 
      nremoveprimes++;
  }

  /* We do not free rel_compact[i] as it is freed with my_malloc_free_all */
  bit_vector_clearbit(rel_used, (size_t) i);

  return nremoveprimes;
}



/*****************************************************************************/
/* Code for clique removal.
   A clique is a connected components of the relation R, where R(i1,i2) iff
   i1 and i2 share a prime of weight 2.
   We remove the heaviest cliques.
   Each ideal h contributes to the weight of the cliques depending on its
   weight (see weight_function_clique).
*/

typedef struct {
  float w;
  index_t i;
} comp_t;

static int
compare (const void *v1, const void *v2)
{
  comp_t *w1 = (comp_t*) v1;
  comp_t *w2 = (comp_t*) v2;

  return (w1->w >= w2->w) ? -1 : 1;
}

float
weight_function_clique (HC_T w)
{
#ifdef USE_CAVALLAR_WEIGHT_FUNCTION
  if (w >= 3)
    return ldexpf (1, -(w-1));
  else if (w == 2)
    return 0.25;
  else
    return 0.0;
#else
  if (w >= 3)
    return powf(0.8, (float) (w-2.0));
  else if (w == 2)
    return 0.125;
  else
    return 0.0;
#endif
}

/* Compute connected component of row i for the relation R(i1,i2) if rows
   i1 and i2 share a prime of weight 2.
   Return number of rows of components, and put in w the total weight. */
static index_t
compute_connected_component (index_t i)
{
  index_t *myrel_compact = rel_compact[i], h, k, n = 1;

  bit_vector_setbit(Tbv, (size_t) i); /* mark row as visited */
  while ((h = *myrel_compact++) != UMAX(h))
    {
      if (ideals_weight[h] == 2)
        {
          k = sum2_index[h] - i; /* other row where prime of index h appears */
          if (!bit_vector_getbit(Tbv, (size_t) k)) /* row k not visited yet */
            n += compute_connected_component (k);
        }
      /* we use the multiplier 5 here, so that the average weight (assumed to
         be 1) is in the middle of the Count[10] array */
      w_ccc += 5.0 * weight_function_clique (ideals_weight[h]);
    }
  return n;
}

/* Delete connected component of row i, assuming the bit-vector is set.
   Warning: we might have some H->hashcount[h] = 3, which is decreased
   to 2, but we don't want to treat that case. Thus we check in addition
   that sum2_index[h] <> 0, which only occurs when H->hashcount[h] = 2 initially. */
static index_t
delete_connected_component (index_t i, index_t *nprimes)
{
  index_t *myrel_compact = rel_compact[i], h, k, w = 1;

  bit_vector_clearbit(Tbv, (size_t) i); /* mark row as visited */
  /* bit i of rel_used is cleared in delete_relation below */
  while ((h = *myrel_compact++) != UMAX(h)) {
    if (ideals_weight[h] == 2 && sum2_index[h]) { /* first row that contains ideal of index h */
      k = sum2_index[h] - i; /* other row where prime of index h appears */
      if (bit_vector_getbit(Tbv, (size_t) k) == 1) /* row k was not visited yet */
  w += delete_connected_component (k, nprimes);
    }
  }
  *nprimes -= delete_relation (i);
  return w;
}

static void
cliques_removal (index_t target_excess, index_t *nrels, index_t *nprimes)
{
  int64_t excess = (((int64_t) *nrels) - *nprimes);
  index_t chunk;
  comp_t *tmp = NULL; /* (weight, index) */
  index_t *myrelcompact, i, h;
  double N = 0.0; /* number of rows */
  unsigned int wceil, j, ltmp = 0, alloctmp = 0xff;
#define MAX_WEIGHT 10
  unsigned int Count[MAX_WEIGHT], *pc; /* Count[w] is the # of components of weight >= w */

  if (excess <=  keep || excess <= target_excess)
    return;

  chunk = excess - target_excess;

  /* first collect sums for primes with weight 2, and compute total weight */
  memset (sum2_index, 0, nprimemax * sizeof(index_t));
  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, (size_t) i)) 
    {
      for (myrelcompact = rel_compact[i]; (h = *myrelcompact++) != UMAX(h); )
        if (ideals_weight[h] == 2) 
          sum2_index[h] += i;
      N += 1.0;
    }
  
  ASSERT_ALWAYS(N == (double) *nrels);

  /* now initialize bit table for relations used */
  bit_vector_neg (Tbv, rel_used);
  memset(Count, 0, sizeof(unsigned int) * MAX_WEIGHT);
  tmp = (comp_t *) malloc (alloctmp * sizeof (comp_t));
  ASSERT_ALWAYS (tmp != NULL);

  for (i = 0; i < nrelmax; i++)
  {
    if (!bit_vector_getbit(Tbv, (size_t) i)) 
    {
      w_ccc = 0.0;
      h = compute_connected_component (i);
      wceil = (unsigned int) w_ccc + 1;
      /* don't consider weight-0 components */
      if ((w_ccc > 0.0) && ((wceil >= MAX_WEIGHT) || (Count[wceil] < chunk)))
      {
        if (ltmp >= alloctmp)
        {
          alloctmp = 2 * alloctmp + 1;
          tmp = (comp_t *) realloc (tmp, alloctmp * sizeof (comp_t));
        }
      tmp[ltmp++] = (comp_t) { w_ccc, i };
      if (wceil > MAX_WEIGHT)
        wceil = MAX_WEIGHT;
      for (pc = &(Count[wceil]); pc-- != Count; (*pc)++);
      }
    }
  }
  qsort (tmp, ltmp, sizeof(comp_t), compare);

  for (j = 0; j < ltmp && *nrels > target_excess + *nprimes; j++)
    *nrels -= delete_connected_component (tmp[j].i, nprimes);

  fprintf (stderr, "    deleted %u heavier connected components at %2.2lf\n",
                     j, seconds ());

#if DEBUG >= 1
  fprintf (stderr, "    DEBUG: ltmp=%u chunk=%u target=%u\n",
                   ltmp, chunk, target_excess);
#endif

  free (tmp);
}

/*****************************************************************************/
/* Code for singletons removal.
   Exist in multithread if __sync_sub_and_fetch exists.
*/

#ifndef HAVE_SYNC_FETCH
static void
onepass_singleton_removal (index_t *nrels, index_t *nprimes)
{
  index_t *tab, i;
  index_t nremoverels = 0;
  index_t nremoveprimes = 0;

  for (i = 0; i < nrelmax; i++)
  {
    if (bit_vector_getbit(rel_used, (size_t) i))
    {
      for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++)
      {  
        if (ideals_weight[*tab] == 1) 
        {
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

#else /* ifndef HAVE_SYNC_FETCH */

typedef struct {
  unsigned int nb;
  pthread_t mt;
  index_t begin, end, sup_nrel, sup_npri;
} ti_t;
static ti_t *ti;

/* Hightest criticality for performance. I inline all myself. */
static void
onepass_thread_singleton_removal (ti_t *mti)
{
  index_t *tab, i;
  HC_T *o;
  bv_t j;

  mti->sup_nrel = mti->sup_npri = 0;
  for (i = mti->begin; i < mti->end; i++)
  {
    j = (((bv_t) 1) << (i & (BV_BITS - 1)));

    if (rel_used->p[i>>LN2_BV_BITS] & j)
      for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++)
        if (ideals_weight[*tab] == 1)
        {
          for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++)
          {
            o = &(ideals_weight[*tab]);
            ASSERT(*o);
            if (*o < UMAX(*o) && !__sync_sub_and_fetch(o, 1))
              (mti->sup_npri)++;
          }
          /* rel_compact[i] is not freed , it is freed with my_malloc_free_all*/
          rel_used->p[i>>LN2_BV_BITS] &= ~j;
          (mti->sup_nrel)++;
          break;
        }
  }
  pthread_exit(NULL);
}

static void
onepass_singleton_parallel_removal (unsigned int nb_thread, index_t *nrels,
                                    index_t *nprimes)
{
  pthread_attr_t attr;
  index_t pas, incr;
  unsigned int i;
  int err;

  ti = (ti_t *) malloc (nb_thread * sizeof (ti_t));
  ASSERT_ALWAYS (ti != NULL);
  ti[0].begin = 0;
  pas = (nrelmax / nb_thread) & ((index_t) ~(BV_BITS -1));
  incr = 0;
  for (i = 0, incr = 0; i < nb_thread - 1; ) {
    incr += pas;
    ti[i].nb = i;
    ti[i++].end = incr;
    ti[i].begin = incr;
  }
  ti[i].nb = i;
  ti[i].end = nrelmax;
  pthread_attr_init (&attr);
  pthread_attr_setstacksize (&attr, 1<<16);
  pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
  for (i = 0; i < nb_thread; i++)
    if ((err = pthread_create (&ti[i].mt, &attr, (void *) onepass_thread_singleton_removal, &ti[i])))
    {
    fprintf (stderr, "onepass_singleton_parallel_removal : pthread_create error 1: %d. %s\n", err, strerror (errno));
    exit (1);
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
#endif /* ifdef HAVE_SYNC_FETCH */

static void
remove_all_singletons (index_t *nrels, index_t *nprimes, int64_t *excess)
{
  index_t oldnrels;
  *excess = (((int64_t) *nrels) - *nprimes);
  fprintf (stderr, "  nrels=%"PRid" nprimes=%"PRid" excess=%"PRId64"\n",
                   *nrels, *nprimes, *excess);
  do 
  {
    oldnrels = *nrels;
#ifdef HAVE_SYNC_FETCH
    ASSERT(npt);
    onepass_singleton_parallel_removal(npt, nrels, nprimes);
#else
    onepass_singleton_removal(nrels, nprimes);
#endif
    *excess = (((int64_t) *nrels) - *nprimes);
    fprintf (stderr, "  new_nrels=%"PRid" new_nprimes=%"PRid" excess=%"PRId64""
                     " at %2.2lf\n", *nrels, *nprimes, *excess, seconds());
  } while (oldnrels != *nrels);
}

static void
singletons_and_cliques_removal (index_t *nrels, index_t *nprimes)
{
  index_t oldnrels = 0;
  int64_t oldexcess, excess, target_excess;
  unsigned int count;

  //First step of singletons removal
  remove_all_singletons(nrels, nprimes, &excess);

  if (excess <= 0) /* covers case nrel = nprimes = 0 */
  {
    fprintf(stderr, "number of relations <= number of ideals\n");
    exit (2);
  }

  if ((double) excess < required_excess * ((double) *nprimes))
  {
    fprintf(stderr, "(excess / nprimes) = %.2f < %.2f. See -required_excess "
                    "argument.\n", ((double) excess / (double) *nprimes), 
                    required_excess);
    exit (2);
  }

  index_t chunk = excess / npass;
  
  //npass pass of clique removal + singletons removal
  for (count = 0; count < npass && excess > 0; count++)
  {
    oldnrels = *nrels;
    oldexcess = excess;
    target_excess = excess - chunk;
    if (target_excess <  keep)
      target_excess = keep;
    fprintf (stderr, "Step %u on %u: target excess is %"PRId64"\n",
                     count, npass, target_excess);
    cliques_removal (target_excess, nrels, nprimes);

    remove_all_singletons(nrels, nprimes, &excess);
    fprintf (stderr, "  [each excess row deleted %2.2lf rows]\n",
                (double) (oldnrels - *nrels) / (double) (oldexcess - excess));
  }


  /* May need an extra pass of clique removal + singletons removal if excess is
    still larger than keep. It may happen due to the fact that each clique does
    not make the excess go down by one but can (rarely) left the excess 
    unchanged. */
  if (excess > keep)
  {
    oldnrels = *nrels;
    oldexcess = excess;
    target_excess = excess - chunk;
    target_excess = keep;

    fprintf (stderr, "Step extra: target excess is %"PRId64"\n", target_excess);
    cliques_removal (target_excess, nrels, nprimes);

    remove_all_singletons(nrels, nprimes, &excess);
    fprintf (stderr, "  [each excess row deleted %2.2lf rows]\n",
                (double) (oldnrels - *nrels) / (double) (oldexcess - excess));
  }
}

/* singleton removal when binary out file (boutfilerel) is requested */
static int
no_singleton(buf_rel_t *br)
{
  weight_t i;
  index_t h;

  for (i = 0; i < br->nb; i++)
  {
    if (ideals_weight[br->primes[i].h] == 1)
    {
      relsup++;

      for (i = 0; i < br->nb; i++)
      {
        h = br->primes[i].h;
        if (ideals_weight[h] != UMAX(weight_t) && !(--(ideals_weight[h])))
          prisup++;
      }
      return 0;
    }
  }
  return 1;
}


/*****************************************************************************/
/* I/O functions */

/* We don't use memory barrier nor (pre)processor orders for portability.
   So, if a ring buffer is empty, one problem exists if the producter
   produces one and the consumer takes it immediatly: the slot might be not
   complety written.
   Same problem exists when the buffer is full in the other sense.
   It's NOT a bug code, but the instructions reordonnancing of the
   optimiser compiler, which in the case of an empty buffer, may
   increase the producter counter BEFORE the end of the complete
   writing of the slot.
   The only << solution >> without synchronous barrier or (pre)processor
   order is a nanosleep after the end of the waiting loop.
   * empty buffer :
       if (A == B + 1) NANOSLEEP;
   * full buffer :
       if (A + 1 == B + SIZEBUF) NANOSLEEP;
   It's very dirty!
*/

/* Callback function called by prempt_scan_relations */

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
        NANOSLEEP;
      else if (cpt_rel_a == cpy_cpt_rel_b)
        pthread_exit(NULL);
    }

    j = (unsigned int) (cpy_cpt_rel_b & (SIZE_BUF_REL - 1));
    my_rel = &(arg->buf_data[j]);

    if (cpt_rel_a == cpy_cpt_rel_b + 1)
      NANOSLEEP;

    if (bit_vector_getbit(arg->rel_used, (size_t) my_rel->num))
    {
      arg->nprimes+=insert_rel_in_table_no_e (my_rel, arg->min_index,
                                     boutfilerel, rel_compact, ideals_weight);
#ifdef STAT
      arg->W += (double) my_rel->nb_above_min_index;
#endif
    }

    test_and_print_progress_now ();
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}

/* Callback function called by prempt_scan_relations */

static void *
thread_print(buf_arg_t *arg)
{
  unsigned int j, aff;
  unsigned long cpy_cpt_rel_b;
  buf_rel_t *my_rel;

  cpy_cpt_rel_b = cpt_rel_b;
  for ( ; ; )
  {
    while (cpt_rel_a == cpy_cpt_rel_b)
      if (!is_finish())
        NANOSLEEP;
      else if (cpt_rel_a == cpy_cpt_rel_b)
          pthread_exit(NULL);

    j = (unsigned int) (cpy_cpt_rel_b & (SIZE_BUF_REL - 1));
    my_rel = &(arg->buf_data[j]);

    if (cpt_rel_a == cpy_cpt_rel_b + 1)
      NANOSLEEP;

    aff = bit_vector_getbit(arg->rel_used, (size_t) my_rel->num);
    if (boutfilerel && aff)
    {
      if (!(no_singleton(my_rel)))
        bit_vector_clearbit(arg->rel_used, (size_t) my_rel->num);
    }
    else if (aff)
    {
      arg->W += (double) my_rel->nb;
      print_relation (arg->f_remaining, my_rel);
    }
    else if (!boutfilerel && arg->f_deleted != NULL)
      print_relation (arg->f_deleted, my_rel);

    test_and_print_progress_now ();
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}


/***************** utils functions for purge binary **************************/

void
read_rel_used_from_infile (const char *infilerel, size_t rel_used_nb_bytes,
                           index_t *nrels)
{
  FILE *in;
  void *p, *pl, *pf;
  index_t nr = 0;

  fprintf (stderr, "Loading rel_used file %s, %zu bytes\n", infilerel,
                   rel_used_nb_bytes);
  if (!(in = fopen_maybe_compressed (infilerel, "r")))
  {
    fprintf (stderr, "Error, cannot open file %s for reading.\n", infilerel);
    exit (1);
  }
  for (p = (void *) rel_used->p, pf = p + rel_used_nb_bytes; (pf - p); )
  {
    pl = p + fread(p, 1, pf - p, in);
    while (p != pl)
      nr += nbbits[*((uint8_t *) p++)];

    if ((p == pf) ^ !feof(in))
    {
      fprintf (stderr, "Error, the length of file %s is incorrect compared "
                       "with nrels (%"PRIu64").\n", infilerel, nrelmax);
      exit (1);
    }
  }
  fclose_maybe_compressed(in, infilerel);
  *nrels = nr;
  fprintf (stderr, "Number of used relations is %"PRid"\n", *nrels);
}



  /* Build the file list (ugly). It is the concatenation of all
   *  b s p
   * where:
   *    b is the basepath (empty if not given)
   *    s ranges over all subdirs listed in the subdirlist (empty if no
   *    such list)
   *    p ranges over all paths listed in the filelist.
   */
char **
filelist_from_file_with_subdirlist(const char *basepath, const char *filelist,
                                   const char *subdirlist)
{
  /* count the number of files in the filelist */
  int nfiles = 0;
  int nsubdirs = 0;
  char ** fl = filelist_from_file (NULL, filelist, 0);
  for(char ** p = fl ; *p ; p++, nfiles++);

  char ** sl = filelist_from_file (basepath, subdirlist, 1);
  for(char ** p = sl ; *p ; p++, nsubdirs++);

  fic = (char **) malloc ((nsubdirs * nfiles + 1) * sizeof (char *));
  ASSERT_ALWAYS(fic != NULL);

  char ** full = fic;
  for(char ** f = fl ; *f ; f++) {
    for(char ** s = sl ; *s ; s++, full++) {
      int ret = asprintf(full, "%s/%s", *s, *f);
      ASSERT_ALWAYS(ret >= 0);
    }
  }
  *full=NULL;
  filelist_clear(fl);
  filelist_clear(sl);
  return fic;
}

static void
usage (const char *argv0)
{
  fprintf (stderr, "Usage: %s [options] ", argv0);
  fprintf (stderr, "[ -filelist <fl> [-basepath <path>] [-subdirlist <sl>] ");
  fprintf (stderr, "| file1 ... filen ]\n");
  fprintf (stderr, "Mandatory command line options: \n");
  fprintf (stderr, "    -out outfile  - write remaining relations in outfile\n");
  fprintf (stderr, "    -nrels nnn    - number of initial relations\n");
  fprintf (stderr, "    -nprimes nnn  - number of prime ideals in renumber table\n");
  fprintf (stderr, "    -minindex nnn - purge primes with index >= nnn\n");
  fprintf (stderr, "\nOther command line options: \n");
  fprintf (stderr, "    -outdel file - output file for deleted relations\n");
  fprintf (stderr, "    -keep    nnn - prune if excess > nnn (default 160)\n");
  fprintf (stderr, "    -npthr   nnn - threads number for suppress singletons\n");
  fprintf (stderr, "    -inprel  file_rel_used - load active relations\n");
  fprintf (stderr, "    -outrel  file_rel_used - write active relations\n");
  fprintf (stderr, "    -npass   nnn - number of step of clique removal (default %d)\n", DEFAULT_NPASS);
  fprintf (stderr, "    -required_excess nnn - percentage of excess required at the end of the first singleton removal step (default %.2f)\n",
  DEFAULT_REQUIRED_EXCESS);
  fprintf (stderr, "    -path_antebuffer <dir> - where is antebuffer\n");
  exit (1);
}

void
purge_parse_npt_param (param_list pl)
{
  const char * snpt = param_list_lookup_string(pl, "npthr");
  if (snpt) 
  {
    char *p, oldp;
    if ((p = strchr(snpt, 'x'))) 
    {
      unsigned int x, y;
      oldp = *p;
      *p = 0;
      if (sscanf(snpt, "%u", &x) && sscanf(&p[1], "%u", &y))
        npt = x * y;
      else
      {
        *p = oldp;
        fprintf (stderr, "Malformed -npthr option: %s\n", snpt);
        usage(argv0);
      }
    } 
    else if (!sscanf(snpt, "%u", &npt))
    {
      fprintf (stderr, "Malformed -npthr option: %s\n", snpt);
      usage(argv0);
    }
  }
}

/*************************** main ********************************************/

int
main (int argc, char **argv)
{
  argv0 = argv[0];
  int k;
  size_t rel_used_nb_bytes;
  param_list pl;
  buf_arg_t buf_arg;
  buf_rel_t *buf_rel;
  uint64_t min_index = UMAX(uint64_t);
  index_t nrels, nprimes;
  size_t tot_alloc_bytes = 0, cur_alloc;

  double wct0 = wct_seconds ();
  fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
  for (k = 1; k < argc; k++)
    fprintf (stderr, " %s", argv[k]);
  fprintf (stderr, "\n");

  param_list_init(pl);

  argv++,argc--;

  /* read all command-line parameters */
  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) continue;
    /* Since we accept file names freeform, we decide to never abort
     * on unrecognized options */
    if (!strcmp(*argv, "--help")) usage(argv0);
    break;
  }

  /* read command-line parameters */
  param_list_parse_uint64(pl, "nrels", &nrelmax);
  param_list_parse_uint64(pl, "nprimes", &nprimemax);
  param_list_parse_int64(pl, "keep", &keep);
  param_list_parse_uint64(pl, "minindex", &min_index);
  purge_parse_npt_param (pl);
  param_list_parse_uint(pl, "npass", &npass);
  param_list_parse_double(pl, "required_excess", &required_excess);
  const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
  const char * filelist = param_list_lookup_string(pl, "filelist");
  const char * basepath = param_list_lookup_string(pl, "basepath");
  const char * subdirlist = param_list_lookup_string(pl, "subdirlist");
  const char * purgedname = param_list_lookup_string(pl, "out");
  const char * infilerel = param_list_lookup_string(pl, "inrel");
  const char * outfilerel = param_list_lookup_string(pl, "outrel");
  const char * deletedname = param_list_lookup_string(pl, "outdel");

  set_antebuffer_path (argv0, path_antebuffer);

  /* Check command-line parameters */
  binfilerel = (infilerel != 0);
  boutfilerel = (outfilerel != 0);

  if (param_list_warn_unused(pl)) {
    fprintf (stderr, "Unused options in command-line\n");
    usage(argv0);
  }
  if ((basepath || subdirlist) && !filelist) {
    fprintf(stderr, "-basepath / -subdirlist only valid with -filelist\n");
    usage(argv0);
  }

  if ((filelist != NULL) + (argc != 0) != 1) {
    fprintf(stderr, "Provide either -filelist or freeform file names\n");
    usage(argv0);
  }

  if (!filelist) // If no filelist was given, files are on the command-line
    fic = argv;
  else if (!subdirlist) //If no subdirlist was given, fic is easy to construct
    fic = filelist_from_file (basepath, filelist, 0);
  else //with subdirlist is a little bit trickier
    fic = filelist_from_file_with_subdirlist(basepath, filelist, subdirlist);

  if (nrelmax == 0)
  {
    fprintf (stderr, "Error, missing -nrels ... option (or nrels=0)\n");
    usage (argv0);
  }
  if (nprimemax == 0)
  {
    fprintf (stderr, "Error, missing -nprimes ... option (or nprimes=0)\n");
    usage (argv0);
  }
  if (min_index > nprimemax)
  {
    fprintf (stderr, "Error, missing -minindex ... option (or > nprimes)\n");
    usage (argv0);
  }
  /* If nrels or nprimes > 2^32, then we need index_t to be 64-bit */
  if (((nprimemax >> 32) != 0 || (nrelmax >> 32) != 0) && sizeof(index_t) < 8)
  {
    fprintf (stderr, "Error, -nrels or -nprimes is too large for a 32-bit "
                     "program\nSee #define index_size in typedefs.h\n");
    exit(1);
  }
  ASSERT_ALWAYS (min_index <= nprimemax);

  /* Printing relevant information */
  fprintf (stderr, "Weight function used during clique removal:\n"
                   "  0     1     2     3     4     5     6     7\n");
  for (k = 0; k < 8; k++)
    fprintf (stderr, "%0.3f ", weight_function_clique ((weight_t) k));
  fprintf (stderr, "\n");

  fprintf (stderr, "Number of relations is %"PRIu64"\n", nrelmax);
  fprintf (stderr, "Number of prime ideals below the two large prime bounds: "
                   "%"PRIu64"\n", nprimemax);

  /* Allocating memory */

  cur_alloc = nprimemax * sizeof (weight_t);
  ideals_weight = (weight_t *) malloc (cur_alloc);
  ASSERT_ALWAYS (ideals_weight != NULL);
  tot_alloc_bytes += cur_alloc;
  fprintf (stderr, "Allocated ideals_weight of %zuMb (total %zuMb so far)\n",
                    cur_alloc >> 20, tot_alloc_bytes >> 20);
  memset(ideals_weight, 0, cur_alloc);

  rel_used_nb_bytes = (nrelmax + 7) >> 3;
  bit_vector_init(rel_used, nrelmax);
  ASSERT_ALWAYS (rel_used->p != NULL);
  tot_alloc_bytes += rel_used_nb_bytes;
  fprintf (stderr, "Allocated rel_used of %zuMB (total %zuMB so far)\n",
                   rel_used_nb_bytes >> 20, tot_alloc_bytes >> 20);

  if (binfilerel)
    read_rel_used_from_infile(infilerel, rel_used_nb_bytes, &nrels);
  else
  {
    bit_vector_set(rel_used, 1);
    nrels = nrelmax;
  }
  /* For the last byte of rel_used, put the bits to 0 if it does not
   * correspond to a rel num */
  if (nrelmax & (BV_BITS - 1))
    rel_used->p[nrelmax>>LN2_BV_BITS] &= 
                                (((bv_t) 1)<<(nrelmax & (BV_BITS - 1))) - 1;

  if (!boutfilerel)
  {
    bit_vector_init(Tbv, nrelmax);
    ASSERT_ALWAYS (Tbv->p != NULL);
    tot_alloc_bytes += rel_used_nb_bytes;
    fprintf (stderr, "Allocated Tbv of %zuMB (total %zuMB so far)\n",
                     rel_used_nb_bytes >> 20, tot_alloc_bytes >> 20);

    cur_alloc = nprimemax * sizeof (index_t);
    sum2_index = (index_t *) malloc (cur_alloc);
    ASSERT_ALWAYS (sum2_index != NULL);
    tot_alloc_bytes += cur_alloc;
    fprintf (stderr, "Allocated sum2_index of %zuMb (total %zuMb so far)\n",
                     cur_alloc >> 20, tot_alloc_bytes >> 20);

    cur_alloc = nrelmax * sizeof (index_t *);
    rel_compact = (index_t **) malloc (cur_alloc);
    ASSERT_ALWAYS (rel_compact != NULL);
    tot_alloc_bytes += cur_alloc;
    fprintf (stderr, "Allocated rel_compact of %zuMb (total %zuMb so far)\n",
                     cur_alloc >> 20, tot_alloc_bytes >> 20);
  }

  cur_alloc = SIZE_BUF_REL * sizeof (buf_rel_t);
  buf_rel = (buf_rel_t *) malloc (cur_alloc);
  ASSERT_ALWAYS (buf_rel != NULL);
  tot_alloc_bytes += cur_alloc;
  fprintf (stderr, "Allocated buf_rel of %zuMb (total %zuMb so far)\n",
                    cur_alloc >> 20, tot_alloc_bytes >> 20);

  memset (&buf_arg, 0, sizeof(buf_arg_t));
  buf_arg.f_deleted = NULL;
  buf_arg.min_index = (index_t) min_index;
  buf_arg.buf_data = buf_rel;
  buf_arg.rel_used = rel_used;
  buf_arg.needed = NEEDED_HMIN;


  /**********************Start interessing stuff *****************************/
  if (!boutfilerel) 
    fprintf (stderr, "Pass 1, reading and storing ideals with index h >= "
                     "%"PRIu64"\n", min_index);
  else
    fprintf (stderr, "Pass 1, reading ideals with index h >= %"PRIu64"\n", 
                     min_index);

  /* first pass over relations in files */
  prempt_scan_relations (fic, &thread_insert, &buf_arg, NULL);
  nprimes = buf_arg.nprimes;
  
  tot_alloc_bytes += get_my_malloc_bytes();
  fprintf (stderr, "Allocated rel_compact[i] %zuMB (total %zuMB so far)\n",
                   get_my_malloc_bytes() >> 20, tot_alloc_bytes >> 20);
#ifdef STAT
  {
  size_t tmp = ((uint64_t) buf_arg.W + nrelmax) * sizeof(index_t);
  double ratio = 100.0 * (double) (((double) tmp) / 
                                   ((double) get_my_malloc_bytes() ));
  fprintf (stderr, "STAT: W_active=%1.0f\nSTAT: Should take %zuMB in memory, "
                   "take %zuMB (%.2f %%)\n", buf_arg.W, tmp >> 20,
                   get_my_malloc_bytes() >> 20, ratio);
  }
#endif

  if (buf_arg.nrels != nrelmax) 
  {
    fprintf (stderr, "Error, -nrels value should match the number of scanned "
                     "relations\nexpected %"PRIu64" relations, found %"PRid"\n",
                     nrelmax, buf_arg.nrels);
    exit (1);
  }

  if (!boutfilerel) 
  {
    singletons_and_cliques_removal (&nrels, &nprimes);
    if (nrels < nprimes)
    {
      fprintf(stderr, "number of relations < number of ideals\n");
      exit (2);
    }
    if (nrels == 0 || nprimes == 0)
    {
      fprintf(stderr, "number of relations or number of ideals is 0\n");
      exit (2);
    }


    /* free rel_compact[i] and rel_compact. We do not need it anymore */
    my_malloc_free_all ();
    tot_alloc_bytes -= get_my_malloc_bytes();
    fprintf (stderr, "Freed rel_compact[i] %zuMB (total %zuMB so far)\n",
                     get_my_malloc_bytes() >> 20, tot_alloc_bytes >> 20);

    free (rel_compact);
    size_t tmp = (nrelmax * sizeof(index_t*));
    tot_alloc_bytes -= tmp;
    fprintf (stderr, "Freed rel_compact %zuMB (total %zuMB so far)\n",
                     tmp >> 20, tot_alloc_bytes >> 20);
  }
  else
  {
    if (!binfilerel)
      fprintf (stderr, "   nrels=%"PRid" nprimes=%"PRid" excess=%"PRId64"\n",
                       nrels, nprimes, ((int64_t) nrels) - nprimes);
    else
      fprintf (stderr, "   nrels=%"PRid" (out of %"PRIu64") nprimes=%"PRid" "
                       "excess=%"PRId64"\n", nrels, nrelmax, nprimes, 
                       ((int64_t) nrels) - nprimes);
  }

  relsup = 0;
  prisup = 0;

  /* reread the relation files and convert them to the new coding */
  fprintf (stderr, "Storing remaining relations...\n");
  
  buf_arg.min_index = 0;
  if (!(buf_arg.f_remaining = fopen_maybe_compressed (purgedname, "w")))
  {
    fprintf (stderr, "Error, cannot open file %s for writing.\n", purgedname);
    exit (1);
  }
  if (deletedname != NULL)
    if (!(buf_arg.f_deleted = fopen_maybe_compressed (deletedname, "w")))
    {
      fprintf (stderr, "Error, cannot open file %s for writing.\n",deletedname);
      exit (1);
    }


  // compute last index i such that ideals_weight[i] != 0
  {
    index_t last_used = nprimemax - 1;
    while (ideals_weight[last_used] == 0)
      last_used--;

    fprintf (buf_arg.f_remaining, "# %"PRid" %"PRid"\n", nrels, last_used + 1);
  }
  
  /* second pass over relations in files */
  buf_arg.needed = NEEDED_ABH;
  prempt_scan_relations (fic, &thread_print, &buf_arg, NULL);


  if (!boutfilerel)
  {
    /* write final values to stdout */
    fprintf(stdout, "Final values:\nnrels=%"PRid" nprimes=%"PRid" "
                    "excess=%"PRId64"\nweight=%1.0f weight*nrels=%1.2e\n", 
                    nrels, nprimes, ((int64_t) nrels) - nprimes, buf_arg.W, 
                    buf_arg.W * (double) nrels);
    fflush (stdout);
  }
  else
  {
    FILE *out;
    void *p, *pf;

    /* For the last byte of rel_used, put the bits to 0 if it does not
    * correspond to a rel num */
    if (nrelmax & (BV_BITS - 1))
      rel_used->p[nrelmax>>LN2_BV_BITS] &= 
                                  (((bv_t) 1)<<(nrelmax & (BV_BITS - 1))) - 1;

    fprintf (stderr, "Deleted %"PRid" relations with singleton(s)\n"
                     "Number of primes deleted: %"PRid"\n"
                     "Writing rel_used file %s, %zu bytes\n",
                     relsup, prisup, outfilerel, rel_used_nb_bytes);

    if (!(out = fopen_maybe_compressed (outfilerel, "w"))) 
    {
      fprintf (stderr, "Error, cannot open file %s for writing.\n", outfilerel);
      exit (1);
    }
    p = (void *) rel_used->p;
    for (pf = p + rel_used_nb_bytes; p != pf; p += fwrite(p, 1, pf - p, out));
    
    fclose_maybe_compressed (out, outfilerel);
  }


  /* Free allocated stuff */
  if (ideals_weight != NULL)
    free(ideals_weight);
  ideals_weight = NULL;
  if (buf_rel != NULL)
    free(buf_rel);
  buf_rel = NULL;
  if (sum2_index != NULL)
    free(sum2_index);
  sum2_index = NULL;

  bit_vector_clear(rel_used);
  if (!boutfilerel)
    bit_vector_clear(Tbv);

  if (filelist)
    filelist_clear(fic);

  fclose_maybe_compressed (buf_arg.f_remaining, purgedname);
  if (buf_arg.f_deleted != NULL)
    fclose_maybe_compressed (buf_arg.f_deleted, deletedname);

  param_list_clear(pl);

  /* print usage of time and memory */
  print_timing_and_memory (wct0);

  return 0;
}
