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

#if (defined FOR_FFS) || (defined FOR_NFS_DL)
#define FOR_DL
#endif

#include "mod_ul.c"
#include "portability.h"
#include "utils.h"
#include "typedefs.h"
#ifdef FOR_DL
#include "utils_ffs.h"
#endif

#define DEBUG 1
//#define STAT_FFS

//#define USE_CAVALLAR_WEIGHT_FUNCTION

#define MAX_FILES 1000000
#define DEFAULT_NPASS 50
#define DEFAULT_KEEP 160
#define DEFAULT_REQUIRED_EXCESS 0.1

typedef struct {
  volatile unsigned int ok;
  unsigned int num, end;
} fr_t;


/* Main variables */
static char antebuffer[PATH_MAX];         /* "directory/antebuffer" or "cat" */
static char rep_cado[PATH_MAX];           /* directory of cado */

static index_t **rel_compact  = NULL; /* see above */

#define DECR_SATURATED_WEIGHT(o,n) if (*o < UMAX(*o) && !(--(*o))) n--;
weight_t *ideals_weight = NULL;

index_t *newindex;

static char ** fic;
static char *pmin, *pminlessone;
static FILE *ofile;     /* For the principal file output. */
static bit_vector rel_used, Tbv;
static relation_stream rs;
static index_t *sum; /* sum of row indices for primes with weight 2 */
static cado_poly pol;
static double wct0;
static double W; /* total weight of the matrix (used in second pass) */
static size_t tot_alloc_bytes;
static index_t nrel,
  nprimes,
  relsup,
  prisup,
  newnrel,
  newnprimes;
static uint64_t nrelmax = 0,
                nprimemax = 0,
                keep = DEFAULT_KEEP, /* default maximun final excess */
                min_index = 0;
static int raw = 0;
static unsigned int npt = 4;
static float w_ccc;
static uint8_t binfilerel; /* True (1) if a rel_used file must be read */
static uint8_t boutfilerel;/* True (1) if a rel_used file must be written */

#ifdef FOR_DL
static FILE *ofile2;
static int pipe_2;
#ifdef STAT_FFS
uint64_t __stat_count[11] = {0,0,0,0,0,0,0,0,0,0,0};
uint64_t __stat_nonzero = 0;
#endif
#endif

static const unsigned char ugly[256] = {
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  0,   1,   2,   3,   4,   5,   6,   7,   8,   9, 255, 255, 255, 255, 255, 255,
  255,  10,  11,  12,  13,  14,  15, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255,  10,  11,  12,  13,  14,  15, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };

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

/*
  1<<NFR = Number of computation threads.
  The best is 2 (NFR = 1).
  if you have 2 cores or less, you could try NFR=1.
  if you have >4 cores, or hyperthreading, AND very slow cores (< 2 Ghz): NFR=2
*/
#define NFR (1)
/* 1<<NNFR = Number of sentences computed of (find root + hashkey of each term)
   computed in one pass. If the number is too big, the buffer between these
   threaded computations and insertRelation (which cannot be parallelized,
   because the hashkey is not bijective) is too big (memory waste) and
   the pipe-line is slow to start;
   If the number is too small, the computation threads are too often in
   nanosleep to keep CPU.
   NNFR=8 to 16, the greatest is the fastest for fast processors & fast
   memory; 14 seems to be the best choice (maximal speed for medium memory use).
*/
#define NNFR (14)

/* Size of relations buffer between parsing & insertRelation.
   CAREFUL! T_REL must be greater (at least double) than
   (1<<(NFR+NNFR+(NFR==0))).
   Stocks the sentences precomputed but not insered. About
   64K sentences for the optimal.
*/
#define T_REL MAX((1<<(NFR+NNFR+1+(NFR==0))),(1<<6))

/* The realistic minimal non-CPU waiting with nanosleep is about
   10 to 40 µs (1<<13 for nanosleep).
   But all the I/O between the threads have been buffered,
   and a thread do a nanosleep only if its buffer is empty.
   So I use here ~2ms (1<<21) to optimize CPU scheduler.
   Max pause is about 4 to 8ms (1<<22, 1<<23); after the program
   is slow down.
*/
static const struct timespec wait_classical = { 0, 1<<21 };

#define NB_PRIMES_OPT 31

typedef struct {
  index_t h;
  exponent_t e;
} prime_t;

typedef struct {
  int64_t a;
  uint64_t b;
  prime_t *primes; /*if nb<=NB_PRIME_OPT, contains the address of primes_data*/
  prime_t primes_data[NB_PRIMES_OPT];
  weight_t nb;           /* number of primes */
  weight_t nb_alloc;     /* allocated space for primes */
  weight_t nb_above_min_index; /* nb of primes above min_index, must be <=nb */
  index_t num;          /* Relation number */
} buf_rel_t;


static buf_rel_t *buf_rel;

/* For the multithread sync */
static volatile unsigned long cpt_rel_a;
static volatile unsigned long cpt_rel_b;
static volatile unsigned int end_insertRelation = 0;

/* copied from utils/antebuffer.c */
#ifndef HAVE_NANOSLEEP
  int nanosleep(const struct timespec *req, struct timespec *rem) {
    if (rem == NULL) {
      /* Dummy to shut up the warning */
    }
#ifdef HAVE_USLEEP
    unsigned long usec = req->tv_sec * 1000000UL + req->tv_nsec / 1000UL;
    usleep(usec);
#else
    sleep(req->tv_sec);
#endif
    return 0;
  }
#endif

/********************** own memory allocation routines ***********************/

/* Rationale: calling one malloc() for each relation in insertNormalRelation()
   is expensive, since malloc() allocates some extra information to keep track
   of every memory blocks. Instead, we allocate memory in big blocks of size
   BLOCK_SIZE. */

#define BLOCK_SIZE (1<<20)  /* memory blocks are allocated of that # of index_t's */
/* relcompact_list is a list of blocks, each one of BLOCK_SIZE index_ts */
static index_t **relcompact_list = NULL;
static index_t *myrelcompact;
static unsigned int relcompact_used = BLOCK_SIZE; /* usage of current block */
static unsigned int relcompact_size = 0;  /* minimal size of relcompact_list */
static int relcompact_current = -1; /* index of current block */

/* return a pointer to an array of n (index_t) */
static index_t *
my_malloc (unsigned int n)
{
  index_t *ptr;

  if (relcompact_used + n > BLOCK_SIZE) {
    relcompact_used = 0;
    if (((unsigned int) (++relcompact_current)) == relcompact_size) {
      relcompact_size = relcompact_size ? (relcompact_size << 1) : (1<<16);
      if (!(relcompact_list = (index_t **) realloc(relcompact_list, relcompact_size * sizeof(index_t *)))) {
    fprintf (stderr, "my_malloc_int: realloc error : %s\n", strerror (errno));
    exit (1);
  }
      }
    SMALLOC(relcompact_list[relcompact_current], BLOCK_SIZE, "my_malloc_int 1");
    myrelcompact = relcompact_list[relcompact_current];
  }
  ptr = &(myrelcompact[relcompact_used]);
  relcompact_used += n;
  tot_alloc_bytes += n * sizeof (index_t);
  return ptr;
}

static void
my_malloc_free_all (void)
{
  while (relcompact_current >= 0) {
    SFREE (relcompact_list[relcompact_current]);
    relcompact_current--;
  }
  SFREE (relcompact_list);
  relcompact_list = NULL;
  relcompact_used = BLOCK_SIZE;
  relcompact_size = 0;
}

/*****************************************************************************/

/* Print the relation 'rel' in matrix format, i.e., a line of the form:

   a,b:t_1,t_2,...,t_k

   a (signed decimal) is a
   b (nonnegative decimal) is b
   t_1 ... t_k (hexadecimal) are the indices of the ideals (starting at 0)

   Return the weight of the relation.
*/

#define FFSCOPYDATA(E)       \
  t = p - op;                \
  for (j = (unsigned int) ((E) - 1); j--; p += t) memcpy (p, op, t)

static int
fprint_rel_row (FILE *file, buf_rel_t *my_buf_rel)
{
  char buf[1<<12], *p;
  unsigned int nb_coeff = 0;
  char *op;
  size_t t;
  unsigned int i, j;
  index_t h;
  exponent_t e;

  p = d64toa10(buf, my_buf_rel->a);
  *p++ = ',';
  p = u64toa10(p, my_buf_rel->b);
  *p++ = ':';


  for (i = 0; i < my_buf_rel->nb; i++)
  {
    e = my_buf_rel->primes[i].e;
    h = newindex[my_buf_rel->primes[i].h];

    op = p;
    p = u64toa16(p, (uint64_t) h);
    *p++ = ',';
    FFSCOPYDATA(e);
  }
  nb_coeff += i;

  *(--p) = '\n';
  p[1] = 0;
  fputs(buf, file);

  return nb_coeff;
}

 /* Write relation j from buffer to rel_compact
    We put in rel_compact only primes such that their index h is greater or
    equal to min_index
 */

static inline void
insertNormalRelation (unsigned int j)
{
  buf_rel_t *my_br;
  unsigned int i, itmp;
  index_t *my_tmp;
  index_t h;

  itmp = 0;
  my_br = &(buf_rel[j]);
  my_tmp = boutfilerel ? NULL : my_malloc(my_br->nb_above_min_index);

  for (i = 0; i < my_br->nb; i++)
  {
    h =  my_br->primes[i].h;
    if (ideals_weight[h] == 0)
    {
      ideals_weight[h] = 1;
      nprimes++;
    }
    else if (ideals_weight[h] != UMAX (weight_t))
      ideals_weight[h]++;

    if (!boutfilerel && h >= min_index)
      my_tmp[itmp++] = h;
  }

  if (!boutfilerel)
  {
    my_tmp[itmp] = UMAX(*my_tmp); /* sentinel */
    rel_compact[my_br->num] = my_tmp;
  }
}

/* Delete a relation: set rel_used[i] to 0, update the count of primes
   in that relation.
   Warning: we only update the count of primes that we consider, i.e.,
   primes with index >= min_index.
*/
static void
delete_relation (index_t i)
{
  index_t *tab;
  HC_T *o;

  for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++)
  {
    o = &(ideals_weight[*tab]);
    ASSERT(*o);
    DECR_SATURATED_WEIGHT(o, newnprimes);
  }

  /* We do not free rel_compact[i] as it is freed with my_malloc_free_all */
  bit_vector_clearbit(rel_used, (size_t) i);
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
          k = sum[h] - i; /* other row where prime of index h appears */
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
   that sum[h] <> 0, which only occurs when H->hashcount[h] = 2 initially. */
static index_t
delete_connected_component (index_t i)
{
  index_t *myrel_compact = rel_compact[i], h, k, w = 1;

  bit_vector_clearbit(Tbv, (size_t) i); /* mark row as visited */
  /* bit i of rel_used is cleared in delete_relation below */
  while ((h = *myrel_compact++) != UMAX(h)) {
    if (ideals_weight[h] == 2 && sum[h]) { /* first row that contains ideal of index h */
      k = sum[h] - i; /* other row where prime of index h appears */
      if (bit_vector_getbit(Tbv, (size_t) k) == 1) /* row k was not visited yet */
  w += delete_connected_component (k);
    }
  }
  delete_relation (i);
  return w;
}

static void
deleteHeavierRows (unsigned int npass)
/* int *nrel, int *nprimes, int nrelmax, int keep)
   &newnrel, &newnprimes, nrelmax, keep */
{
  static unsigned int count = 0;
  static index_t chunk;
  comp_t *tmp = NULL; /* (weight, index) */
  index_t *myrelcompact, i, h;
  double N = 0.0; /* number of rows */
  unsigned int wceil, j, ltmp = 0, alloctmp = 0xff;
  long target;
#define MAX_WEIGHT 10
  unsigned int Count[MAX_WEIGHT], *pc; /* Count[w] is the # of components of weight >= w */

  if (newnrel - newnprimes <= keep)
    return;
  if (!count)
    chunk = (newnrel - newnprimes) / npass;

  /* first collect sums for primes with weight 2, and compute total weight */
  MEMSETZERO(sum, nprimemax);
  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, (size_t) i)) {
      for (myrelcompact = rel_compact[i]; (h = *myrelcompact++) != UMAX(h); )
  if (ideals_weight[h] == 2) sum[h] += i;
      N += 1.0;
    }
  fprintf (stderr, "Step %u on %u: Matrix has %1.0f real (non null) rows\n",
                   count, npass, N);
  ASSERT_ALWAYS(N == (double) newnrel);

  /* now initialize bit table for relations used */
  bit_vector_init(Tbv, (size_t) nrelmax);
  bit_vector_neg (Tbv, rel_used);
  memset(Count, 0, sizeof(unsigned int) * MAX_WEIGHT);
  SMALLOC(tmp, alloctmp, "deleteHeavierRows 2");

  for (i = 0; i < nrelmax; i++)
    {
    if (!bit_vector_getbit(Tbv, (size_t) i)) {
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

  /* remove heaviest components, assuming each one decreases the excess by 1;
     we remove only part of the excess at each call of deleteHeavierRows,
     hoping to get "better" components to remove at the next call. */
  if (++count < npass)
    {
      target = ((long) newnrel) - newnprimes - chunk;
      if (target < (long) keep) target = keep;
    }
  else
    target = keep; /* enough steps */

#if DEBUG >= 1
  fprintf (stderr, "DEBUG: newnrel=%u newnprimes=%u\n"
                   "DEBUG: ltmp=%u chunk=%u target=%lu\n", newnrel,
                   newnprimes, ltmp, chunk, target);
#endif

  for (j = 0; j < ltmp && newnrel > target + newnprimes; j++)
    newnrel -= delete_connected_component (tmp[j].i);
  fprintf (stderr, "deleted %u heavier connected components at %2.2lf\n",
                     j, seconds ());
  bit_vector_clear(Tbv);
  free (tmp);
}

/*****************************************************************************/
/* Code for singletons removal.
   Exist in multithread if __sync_sub_and_fetch exists.
*/

#ifndef HAVE_SYNC_FETCH
static void
onepass_singleton_removal ()
{
  index_t *tab, i;

  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, (size_t) i))
      for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++)
  if (ideals_weight[*tab] == 1) {
    delete_relation(i);
    newnrel--;
    break;
  }
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
            //fprintf (stderr, "i=%u h=%x o=%d\n", i, *tab, *o);
            ASSERT(*o);
            if (*o < UMAX(*o) && !__sync_sub_and_fetch(o, 1))
              (mti->sup_npri)++;
          }
          /* rel_compact[i] = NULL; */
          rel_used->p[i>>LN2_BV_BITS] &= ~j;
          (mti->sup_nrel)++;
          break;
        }
  }
  pthread_exit(NULL);
}

static void
onepass_singleton_parallel_removal (unsigned int nb_thread)
{
  pthread_attr_t attr;
  index_t pas, incr;
  unsigned int i;
  int err;

  SMALLOC(ti, nb_thread, "onepass_singleton_parallel_removal :");
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
    newnrel -= ti[i].sup_nrel;
    newnprimes -= ti[i].sup_npri;
  }
  pthread_attr_destroy(&attr);
  SFREE(ti);
}
#endif /* ifdef HAVE_SYNC_FETCH */

static void
remove_singletons (unsigned int npass, double required_excess)
{
  index_t oldnewnrel = 0, oldtmpnewnrel = 0;
#if HR == 32
  int32_t oldexcess = 0, excess;
#else
  int64_t oldexcess = 0, excess;
#endif
  int count = 0;

  if (!binfilerel) newnrel = nrel;
  newnprimes = nprimes;
  excess = ((long) newnrel) - newnprimes;
  for ( ; newnrel != oldnewnrel || excess > (long) keep ; ) {
    /* delete heavy rows when we have reached a fixed point */
    if (newnrel == oldnewnrel) {
      /* check we have enough excess initially (at least required_excess) */
      if (count++ == 0 && (double) excess < required_excess*(double)newnprimes)
        {
          fprintf(stderr, "excess < %.2f * #primes. See -required_excess "
                          "argument.\n", required_excess);
          exit (2);
        }
      if (oldexcess > excess)
  fprintf (stderr, "   [each excess row deleted %2.2lf rows]\n",
     (double) (oldtmpnewnrel - newnrel) / (double) (oldexcess - excess));
      oldexcess = excess;
      oldtmpnewnrel = newnrel;
      deleteHeavierRows (npass);
    }
    oldnewnrel = newnrel;
#ifdef HAVE_SYNC_FETCH
    ASSERT(npt);
    onepass_singleton_parallel_removal(npt);
#else
    onepass_singleton_removal();
#endif
    excess = ((long) newnrel) - newnprimes;
    fprintf (stderr, "   new_nrows=%lu new_ncols=%lu (%ld) at %2.2lf\n",
       (unsigned long) newnrel, (unsigned long) newnprimes, (long) excess, seconds());
  }
  /* Warning: we might have empty rows, i.e., empty rel_compact[i] lists,
     if all primes in a relation are less than minpr or minpa. */
  nrel = newnrel;
  nprimes = newnprimes;
}

/*****************************************************************************/
/* This function renumbers used primes (those with H->hashcount[i] > 1)
   and puts the corresponding index in H->hashcount[i].

   At return, nprimes is the number of used primes.

   We locate used primes and do not try to do fancy things as sorting w.r.t.
   weight, since this will probably be done later on.
   All rows will be 1 more that needed -> subtract 1 in fprint...! */

static void
renumber (const char *sos)
{
  FILE *fsos = NULL;
  index_t i, nb = 0;
  static int count = 0;

  SMALLOC(newindex, nprimemax, "renumber 1");

  if (sos != NULL)
  {
    fprintf (stderr, "Output renumber table in file %s\n", sos);
    fsos = fopen_maybe_compressed (sos, "w");
    fprintf (fsos, "# each row contains 2 hexadecimal values: n i\n"
    "# n is the new ideal index (for merge and replay)\n"
    "# i is the ideal index value in the renumber file\n");
  }

  for (i = 0; i < nprimemax; i++)
  {
    if (ideals_weight[i])
    {
    /* Since we consider only primes >= minpr or minpa,
       smaller primes might appear only once here, thus we can't
       assert H->hashcount[i] > 1, but H->hashcount[i] = 1 should
       be rare if minpr/minpa are well chosen (not too large). */
#if DEBUG >= 1
      if (ideals_weight[i] == 1)
      {
        if (i < min_index)
          fprintf (stderr, "Warning: ");
        else
          fprintf (stderr, "WARNING (probably an error): ");

        fprintf (stderr, "singleton prime with index %lu\n", (unsigned long) i);
        count++;
      }
#endif

      if (fsos)
        fprintf(fsos, "%lx %lx\n", (unsigned long) nb, (unsigned long) i);

      newindex[i] = nb++;
    }
    else
    {
      newindex[i] = UMAX(index_t);
    }
  }
    if (fsos)
      fclose_maybe_compressed (fsos, sos);
    newnprimes = nb;
#if DEBUG >= 1
  fprintf (stderr, "Warning: %d singleton primes at the end of purge\n",count);
#endif
}

/* singleton removal when binary out (or in ??) file is requested */
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
static void
prempt_load (prempt_t prempt_data) {
  char **p_files = prempt_data->files, *pprod;
  FILE *f;
  size_t try_load, load;
  char *p, *l, *pmax = &(prempt_data->buf[PREMPT_BUF]);

  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);

  while (*p_files)
  {
    if (!(f = popen (*p_files, "r")))
    {
      fprintf (stderr, "prempt_load: popen error. %s\n", strerror (errno));
      exit (1);
    }
#if DEBUG >= 1
    /* Write the command that read the files (one file per line) */
    p = strchr(*p_files, '/');
    if (p)
    {
      *p = 0;
      fprintf (stderr, "%s\n", *p_files);
      *p = '/';
      for ( ; ; )
      {
        l = strchr(p, ' ');
        if (l)
        {
          *l = 0;
          fprintf(stderr, "   %-70s\n", p);
          *l = ' ';
          p = strchr(&(l[1]), '/');
          if (!p)
          {
            fprintf(stderr, "%s\n", &(l[1]));
            break;
          }
        }
        else
        {
          fprintf(stderr, "   %-70s\n", p);
          break;
        }
      }
    }
#endif

    pprod = (char *) prempt_data->pprod;
    for ( ; ; )
    {
      if ((pprod != prempt_data->pcons) &&
    (((((PREMPT_BUF + ((size_t) prempt_data->pcons))) - ((size_t) pprod)) &
      (PREMPT_BUF - 1)) <= PREMPT_ONE_READ))
        {
          nanosleep(&wait_classical, NULL);
        }
      else
      {
        try_load =
         MIN((PREMPT_BUF + ((size_t)pmin)) - ((size_t) pprod), PREMPT_ONE_READ);
        if ((load = fread (pprod, 1, try_load, f)))
        {
          pprod += load;
          if (pprod == pmax)
            pprod = pmin;
          prempt_data->pprod = pprod;
        }
        else if (feof(f))  // we go to the next batch of files
        {
          pclose(f);
          free(*p_files);
          *p_files++ = NULL;
          break;
        }
        else // error
        {
          fprintf (stderr, "prempt_load: load error (%s) from\n%s\n",
                           strerror (errno), *p_files);
          exit (1);
        }
      }
    }
  }
  prempt_data->end = 1;
  pthread_exit (NULL);
}

static inline void
#ifndef FOR_DL
relation_stream_get_fast (prempt_t prempt_data, unsigned int j)
#else
relation_stream_get_fast (prempt_t prempt_data, unsigned int j, int passtwo)
#endif
{
  buf_rel_t *mybufrel = &(buf_rel[j]);
  int64_t n;
  char *p;
  unsigned int nb_primes_read;
  unsigned long pr;
  unsigned char c, v;
  weight_t nb_above_index = 1; // count the -1 at the end of the relations
                                // in rel_compact
#ifdef FOR_FFS
  unsigned int basis_ab = 16;
#else
  unsigned int basis_ab = 10;
#endif

#define LOAD_ONE(P) { c = *P; P = ((size_t) (P - pminlessone) & (PREMPT_BUF - 1)) + pmin; }

  p = (char *) prempt_data->pcons;

  LOAD_ONE(p);
  if (c == '-') {
    mybufrel->a = -1;
    LOAD_ONE(p);
  }
  else
    mybufrel->a = 1;
  for (n = 0 ; (v = ugly[c]) < basis_ab ; ) {
#ifdef FOR_FFS
    n = (n << 4) + v;
#else
    n = n * 10 + v;
#endif
    LOAD_ONE(p);
  }
  ASSERT_ALWAYS(c == ',');
  mybufrel->a *= n;

  n = 0;
  LOAD_ONE(p);
  for ( ; (v = ugly[c]) < basis_ab ; ) {
#ifdef FOR_FFS
    n = (n << 4) + v;
#else
    n = n * 10 + v;
#endif
    LOAD_ONE(p);
  }
  ASSERT_ALWAYS(c == ':');
  mybufrel->b = n;

  nb_primes_read = 0;
  for ( c = 0 ; ; )
  {
    if (c == '\n')
      break;

    LOAD_ONE(p);
    for (pr = 0 ; (v = ugly[c]) < 16 ; )
    {
      pr = (pr << 4) + v;
      LOAD_ONE(p);
    }
    ASSERT_ALWAYS(c == ',' || c == '\n');

    if (nb_primes_read > 0 && mybufrel->primes[nb_primes_read-1].h == pr)
        mybufrel->primes[nb_primes_read-1].e++;
    else
    {
      if (mybufrel->nb_alloc == nb_primes_read)
      {
        mybufrel->nb_alloc += mybufrel->nb_alloc >> 1;
        if (nb_primes_read == NB_PRIMES_OPT)
        {
          prime_t *p = mybufrel->primes;
          SMALLOC(mybufrel->primes, mybufrel->nb_alloc, "realloc primes");
          memcpy (mybufrel->primes, p, NB_PRIMES_OPT * sizeof(prime_t));
        }
        else
          mybufrel->primes = (prime_t *)
            realloc (mybufrel->primes, mybufrel->nb_alloc * sizeof(prime_t));
        fprintf (stderr, "nb_alloc = %u\n", mybufrel->nb_alloc);
      }

      nb_above_index += (weight_t) (pr >= min_index);
      mybufrel->primes[nb_primes_read++] = (prime_t) { .h = pr, .e = 1};
    }
  }

  mybufrel->nb = nb_primes_read;
  mybufrel->nb_above_min_index = nb_above_index;

#ifdef STAT_FFS
  unsigned int i;
  for (i = 0; i < mybufrel->nb; i++)
  {
    if (bit_vector_getbit(rel_used, (size_t) mybufrel->num))
    {
      __stat_nonzero++;
      if (abs(mybufrel->rel.primes[i].e) > 10)
        __stat_count[0]++;
      else
        __stat_count[abs(mybufrel->rel.rp[i].e)]++;
    }
  }
#endif
}

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
       if (A == B + 1) nanosleep (&wait_classical, NULL);
   * full buffer :
       if (A + 1 == B + SIZEBUF) nanosleep (&wait_classical, NULL);
   It's very dirty!
*/

void
insertRelation()
{
  unsigned int j;
  unsigned long cpy_cpt_rel_b;

  cpy_cpt_rel_b = cpt_rel_b;
  for ( ; ; )
  {
    while (cpt_rel_a == cpy_cpt_rel_b)
    {
      if (!end_insertRelation)
        nanosleep (&wait_classical, NULL);
      //else if (cpt_rel_a == cpy_cpt_rel_b)
      else
        pthread_exit(NULL);
    }

    j = (unsigned int) (cpy_cpt_rel_b & (T_REL - 1));

    if (cpt_rel_a == cpy_cpt_rel_b + 1)
      nanosleep (&wait_classical, NULL);

    if (bit_vector_getbit(rel_used, (size_t) buf_rel[j].num))
      insertNormalRelation (j);

    if (relation_stream_disp_progress_now_p(rs))
      fprintf(stderr, "read useful %lu relations in %.1fs"
                      " -- %.1f MB/s -- %.1f rels/s\n",
                      rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}

static void
printrel()
{
  unsigned int j, aff;
  unsigned long cpy_cpt_rel_b;

  cpy_cpt_rel_b = cpt_rel_b;
  for ( ; ; )
  {
    while (cpt_rel_a == cpy_cpt_rel_b)
      if (!end_insertRelation)
        nanosleep (&wait_classical, NULL);
      else
        if (cpt_rel_a == cpy_cpt_rel_b)
          pthread_exit(NULL);

    j = (unsigned int) (cpy_cpt_rel_b & (T_REL - 1));
    if (cpt_rel_a == cpy_cpt_rel_b + 1)
      nanosleep (&wait_classical, NULL);

    aff = bit_vector_getbit(rel_used, (size_t) buf_rel[j].num);
    if (boutfilerel && aff)
    {
      if (!(no_singleton(&(buf_rel[j]))))
        bit_vector_clearbit(rel_used, (size_t) buf_rel[j].num);
    }
    else if (aff)
      W += (double) fprint_rel_row(ofile, &(buf_rel[j]));
#ifdef FOR_DL
    else if (!boutfilerel)
      fprint_relation_row(ofile2, &(buf_rel[j].rel));
#endif

    if (relation_stream_disp_progress_now_p(rs))
      fprintf(stderr, "re-read & print useful %lu relations in %.1fs"
        " -- %.1f MB/s -- %.1f rels/s\n", rs->nrels, rs->dt, rs->mb_s,
        rs->rels_s);
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}

/* Read all relations from file, and fills the rel_used and rel_compact arrays
   for each relation i:
   - rel_used[i] = 0 if the relation i is deleted
     rel_used[i] = 1 if the relation i is kept (so far)
   - rel_compact is an array, terminated by -1, of pointers to the entries
     in the hash table for the considered primes

     Trick: we only read relations for which rel_used[i]==1.
 */
static int
prempt_scan_relations_pass_one ()
{
  char *pcons, *pcons_old, *pcons_max, *p, **ff;
  pthread_attr_t attr;
  pthread_t thread_load, thread_relation;
  prempt_t prempt_data;
  unsigned long cpy_cpt_rel_a;
  unsigned int length_line, i, k;
  int err;
  char c;

  end_insertRelation = 0;

  SMALLOC (buf_rel, T_REL, "prempt_scan_relations_pass_one 1");
  MEMSETZERO(buf_rel, T_REL);
  for (i = T_REL; i--; ) {
    buf_rel[i].primes = buf_rel[i].primes_data;
    buf_rel[i].nb_alloc = NB_PRIMES_OPT;
  }

  cpt_rel_a = cpt_rel_b = 0;
  cpy_cpt_rel_a = cpt_rel_a;
  nprimes = 0;
  relation_stream_init (rs);
  rs->pipe = 1;
  length_line = 0;

  prempt_data->files = prempt_open_compressed_rs (antebuffer, fic);

  SMALLOC (prempt_data->buf, PREMPT_BUF, "prempt_scan_relations_pass_one 3");
  pmin = prempt_data->buf;
  pminlessone = pmin - 1;
  prempt_data->pcons = pmin;
  prempt_data->pprod = pmin;
  pcons_max = &(prempt_data->buf[PREMPT_BUF]);
  prempt_data->end = 0;
  pthread_attr_init (&attr);
  pthread_attr_setstacksize (&attr, 1<<16);
  pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
  if ((err = pthread_create (&thread_load, &attr, (void *) prempt_load, prempt_data)))
    {
    fprintf (stderr, "prempt_scan_relations_pass_one: pthread_create error 1: %d. %s\n", err, strerror (errno));
    exit (1);
    }
  if ((err = pthread_create (&thread_relation, &attr, (void *) insertRelation, NULL)))
    {
    fprintf (stderr, "prempt_scan_relations_pass_one: pthread_create error 2: %d. %s\n", err, strerror (errno));
    exit (1);
    }

  pcons = (char *) prempt_data->pcons;
  for ( ; ; )
    {
      rs->pos += length_line;
      length_line = 0;
      prempt_data->pcons = pcons;

      while (pcons == prempt_data->pprod)
        if (!prempt_data->end)
          nanosleep (&wait_classical, NULL);
        else if (pcons == prempt_data->pprod)
          goto end_of_files;

      if (pcons == prempt_data->pprod + sizeof(*pcons))
        nanosleep (&wait_classical, NULL);

      rs->lnum++;
      if (*pcons != '#')
      {
        while ((((PREMPT_BUF + ((size_t) prempt_data->pprod)) -
               ((size_t) pcons)) & (PREMPT_BUF - 1)) <= ((unsigned int) RELATION_MAX_BYTES) &&
     !prempt_data->end)
        nanosleep(&wait_classical, NULL);

        if (pcons > prempt_data->pprod)
        {
          p = &(pcons_max[-1]);
          c = *p;
          *p = '\n';
          pcons_old = pcons;
          while (*pcons++ != '\n');
            *p = c;

          length_line = (pcons - pcons_old);
          if (pcons <= p)
            goto testendline;
          pcons = pmin;
          if (c == '\n')
            goto testendline;
        }

        p = &(((char *) prempt_data->pprod)[-1]);
        c = *p;
        *p = '\n';
        pcons_old = pcons;
        while (*pcons++ != '\n');
        *p = c;
        length_line += (pcons - pcons_old);

      testendline:
        if (c != '\n' && pcons == prempt_data->pprod)
        {
          fprintf (stderr, "prempt_scan_relations_pass_one: "
            "the last line has not a carriage return\n");
          exit(1);
        }
        if (length_line > ((unsigned int) RELATION_MAX_BYTES))
        {
          fprintf (stderr, "prempt_scan_relations_pass_one: relation line size"
                           " (%u) is greater than RELATION_MAX_BYTES (%d)\n",
                           length_line, RELATION_MAX_BYTES);
          exit(1);
        }

        while (cpy_cpt_rel_a == cpt_rel_b + T_REL)
          nanosleep(&wait_classical, NULL);

        k = (unsigned int) (cpy_cpt_rel_a & (T_REL - 1));
        if (cpy_cpt_rel_a + 1 == cpt_rel_b + T_REL)
          nanosleep(&wait_classical, NULL);

        buf_rel[k].num = rs->nrels++;

#ifndef FOR_DL
        if (bit_vector_getbit(rel_used, (size_t) buf_rel[k].num))
          relation_stream_get_fast (prempt_data, k);
#else
        relation_stream_get_fast (prempt_data, k, 0);
#endif
        /* Delayed find root computation by block of 1<<NNFR */
        if (cpy_cpt_rel_a && !(k & ((1<<NNFR)-1)))
        {
          if (cpy_cpt_rel_a > (1<<(NFR+NNFR)))
          cpt_rel_a = cpy_cpt_rel_a - (1<<(NFR+NNFR));
        }
        cpy_cpt_rel_a++;
      }
      else
      {
        do
        {
          while (pcons == prempt_data->pprod)
          {
            if (!prempt_data->end)
              nanosleep (&wait_classical, NULL);
            else if (pcons == prempt_data->pprod)
            {
              fprintf (stderr, "prempt_scan_relations_pass_one: at the end of"
                               " files, a line without \\n ?\n");
              exit (1);
            }
          }

          p = ((pcons <= prempt_data->pprod) ? (char *) prempt_data->pprod
                                             : pcons_max) - 1;
          c = *p;
          *p = '\n';
          pcons_old = pcons;
          while (*pcons++ != '\n');
            *p = c;

          length_line += (pcons - pcons_old);
          err = (pcons > p && c != '\n');
          if (pcons == pcons_max)
            pcons = pmin;
        } while (err);
      }
  }

 end_of_files:
  cpt_rel_a = cpy_cpt_rel_a;
  while (cpy_cpt_rel_a != cpt_rel_b)
    nanosleep(&wait_classical, NULL);

  end_insertRelation = 1;
  pthread_join(thread_relation, NULL);
  /* if (pthread_tryjoin_np (thread_load, NULL)) */
  pthread_cancel(thread_load);
  pthread_join(thread_load, NULL);
  pthread_attr_destroy(&attr);

  free (prempt_data->buf);
  for (ff = prempt_data->files; *ff; free(*ff++));
  free (prempt_data->files);
  for (i = T_REL; i-- ; )
    if (buf_rel[i].nb_alloc != NB_PRIMES_OPT ) SFREE(buf_rel[i].primes);

  SFREE (buf_rel);

  relation_stream_trigger_disp_progress(rs);
  fprintf (stderr, "End of read: %lu relations in %.1fs -- %.1f MB/s -- %.1f rels/s\n",
           rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
  relation_stream_clear(rs);

  if (rs->nrels != nrelmax) {
    fprintf (stderr, "Error, -nrels value should match the number of scanned relations\nexpected %lu relations, found %lu\n", (unsigned long) nrelmax, rs->nrels);
    exit (1);
  }

  return 1;
}

/* ReRead all relations from files and output remaining relations
   (those with rel_used[i] <> 0) in ofile.

   For FFS ouput deleted relations (those with rel_used[i] == 0) in ofile2

   If raw is non-zero, output relations in CADO format
   (otherwise in format used by merge).

   Raw or not raw, a ring buffer is used to stock the relations.
   It maybe useless in some cases.
 */
static int
prempt_scan_relations_pass_two (const char *oname,
#ifdef FOR_DL
        const char *oname2,
#else
        bit_vector_srcptr rel_used,
#endif
        index_t nrows, index_t ncols, int raw)
{
  char *pcons, *pcons_old, *pcons_max, *p, **f;
  pthread_attr_t attr;
  pthread_t thread_load, thread_printrel;
  prempt_t prempt_data;
  unsigned long cpy_cpt_rel_a;
  unsigned int length_line, i, k;
  int err;
  char c;

  int pipe;

  ofile = fopen_maybe_compressed2(oname, "w", &pipe, NULL);
#ifdef FOR_DL
  ofile2 = fopen_maybe_compressed2(oname2, "w", &pipe_2, NULL);
#endif
  if (!raw)
    fprintf (ofile, "%lu %lu\n", (unsigned long) nrows, (unsigned long) ncols);
  fprintf (stderr, "Final pass:\n");
  end_insertRelation = 0;

  SMALLOC(buf_rel, T_REL, "prempt_scan_relations_pass_two 1");
  MEMSETZERO(buf_rel, T_REL);
  for (i = T_REL; i--; ) {
    buf_rel[i].nb_alloc = NB_PRIMES_OPT;
    buf_rel[i].primes = buf_rel[i].primes_data;
  }

  cpt_rel_a = cpt_rel_b = 0;
  cpy_cpt_rel_a = cpt_rel_a;
  relation_stream_init (rs);
  rs->pipe = 1;
  length_line = 0;

  prempt_data->files = prempt_open_compressed_rs (antebuffer, fic);

  SMALLOC(prempt_data->buf, PREMPT_BUF, "prempt_scan_relations_pass_two 3");
  pmin = prempt_data->buf;
  pminlessone = pmin - 1;
  prempt_data->pcons = pmin;
  prempt_data->pprod = pmin;
  pcons_max = &(prempt_data->buf[PREMPT_BUF]);
  prempt_data->end = 0;
  pthread_attr_init (&attr);
  pthread_attr_setstacksize (&attr, 1<<16);
  pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);

  if ((err = pthread_create (&thread_load, &attr, (void *) prempt_load, prempt_data)))
    {
    fprintf (stderr, "prempt_scan_relations_pass_two: pthread_create error 1: %d. %s\n", err, strerror (errno));
    exit (1);
    }
  W = 0.0;
  if ((err = pthread_create (&thread_printrel, &attr, (void *) printrel, NULL)))
    {
    fprintf (stderr, "prempt_scan_relations_pass_two: pthread_create error 2: %d. %s\n", err, strerror (errno));
    exit (1);
    }

  pcons = (char *) prempt_data->pcons;
  for ( ; ; )
    {
      rs->pos += length_line;
      length_line = 0;
      prempt_data->pcons = pcons;

      while (pcons == prempt_data->pprod)
  if (!prempt_data->end)
    nanosleep (&wait_classical, NULL);
  else
    if (pcons == prempt_data->pprod)
      goto end_of_files;
      if (pcons == prempt_data->pprod + sizeof(*pcons)) nanosleep (&wait_classical, NULL);

      rs->lnum++;
      if (*pcons != '#')
  {
    while ((((PREMPT_BUF + ((size_t) prempt_data->pprod)) -
       ((size_t) pcons)) & (PREMPT_BUF - 1)) <= ((unsigned int) RELATION_MAX_BYTES) &&
     !prempt_data->end)
        nanosleep(&wait_classical, NULL);

    if (pcons > prempt_data->pprod)
      {
        p = &(pcons_max[-1]);
        c = *p;
        *p = '\n';
        pcons_old = pcons;
        while (*pcons++ != '\n');
        *p = c;
        length_line = (pcons - pcons_old);
        if (pcons <= p)
    goto lineok;
        pcons = pmin;
        if (c == '\n')
    goto lineok;
      }
    p = &(((char *) prempt_data->pprod)[-1]);
    c = *p;
    *p = '\n';
    pcons_old = pcons;
    while (*pcons++ != '\n');
    *p = c;
    length_line += (pcons - pcons_old);
  lineok:
    while (cpy_cpt_rel_a == cpt_rel_b + T_REL) nanosleep(&wait_classical, NULL);
    k = (unsigned int) (cpy_cpt_rel_a & (T_REL - 1));
    if (cpy_cpt_rel_a + 1 == cpt_rel_b + T_REL) nanosleep(&wait_classical, NULL);
    buf_rel[k].num = rs->nrels++;

#ifndef FOR_DL
    if (bit_vector_getbit(rel_used, (size_t) buf_rel[k].num))
      relation_stream_get_fast (prempt_data, k);
#else
    relation_stream_get_fast (prempt_data, k, 1);
#endif
    /* Delayed find root computation by block of 1<<NNFR */
    if (cpy_cpt_rel_a && !(k & ((1<<NNFR)-1)))
      {
        if (cpy_cpt_rel_a > (1<<(NFR+NNFR)))
    cpt_rel_a = cpy_cpt_rel_a - (1<<(NFR+NNFR));
      }
    cpy_cpt_rel_a++;
  }
      else
  {
    do
      {
        while (pcons == prempt_data->pprod)
    {
      if (!prempt_data->end)
        nanosleep (&wait_classical, NULL);
    }
        p = ((pcons <= prempt_data->pprod) ? (char *) prempt_data->pprod
       : pcons_max) - 1;
        c = *p;
        *p = '\n';
        pcons_old = pcons;
        while (*pcons++ != '\n');
        *p = c;
        length_line += (pcons - pcons_old);
        err = (pcons > p && c != '\n');
        if (pcons == pcons_max) pcons = pmin;
      }
    while (err);
  }
    }

 end_of_files:
  cpt_rel_a = cpy_cpt_rel_a;
  while (cpy_cpt_rel_a != cpt_rel_b)
    nanosleep(&wait_classical, NULL);

  end_insertRelation = 1;
  pthread_join(thread_printrel, NULL);
  /* if (pthread_tryjoin_np (thread_load, NULL)) */
  pthread_cancel(thread_load);
  pthread_join(thread_load, NULL);
  pthread_attr_destroy(&attr);

  free (prempt_data->buf);
  for (f = prempt_data->files; *f; free(*f++));
  free (prempt_data->files);
  for (i = T_REL; i-- ; ) {
    if (buf_rel[i].nb_alloc != NB_PRIMES_OPT ) SFREE(buf_rel[i].primes);
  }
  SFREE (buf_rel);

  relation_stream_trigger_disp_progress(rs);
  fprintf (stderr, "End of re-read: %lu relations in %.1fs -- %.1f MB/s -- %.1f rels/s\n",
           rs->nrels, rs->dt, rs->mb_s, rs->rels_s);

#if defined FOR_FFS && defined STAT_FFS
  fprintf (stderr, "# of non zero coeff: %lu\n", __stat_nonzero);
  for (int i = 1; i <= 10 ; i++)
    fprintf (stderr, "# of coeffs of abs value %d: %lu(%.2f%%)\n", i,
             __stat_count[i], 100 * (double) __stat_count[i]/__stat_nonzero);
  fprintf (stderr, "# of coeffs of abs value > 10: %lu(%.2f%%)\n",
           __stat_count[0], 100 * (double) __stat_count[0]/__stat_nonzero);
#endif

  /* write excess to stdout */
  if (!raw)
    printf ("NROWS:%lu WEIGHT:%1.0f WEIGHT*NROWS=%1.2e\n",
                 (unsigned long) nrows, W, W * (double) nrows);
  printf ("EXCESS: %lu\n", ((long) nrows) - ncols);
  fflush (stdout);

  relation_stream_clear(rs);

  if (pipe)  pclose(ofile);  else fclose(ofile);
#ifdef FOR_DL
  if (pipe_2) pclose(ofile2); else fclose(ofile2);
#endif

  return 1;
}

static void
usage (const char *argv0)
{
  fprintf (stderr, "Usage: %s [options] ", argv0);
  fprintf (stderr, "[ -filelist <fl> [-basepath <path>] [-subdirlist <sl>] ");
  fprintf (stderr, "| file1 ... filen ]\n");
  fprintf (stderr, "Mandatory command line options: \n");
  fprintf (stderr, "       -poly polyfile - use polynomial in polyfile\n");
  fprintf (stderr, "       -out outfile   - write remaining relations in outfile\n");
  fprintf (stderr, "       -nrels nnn     - number of initial relations\n");
#ifndef FOR_DL
  fprintf (stderr, "\n    Other command line options: \n");
#endif
  fprintf (stderr, "       -outdel file - output file for deleted relations\n");
  fprintf (stderr, "       -sos sosfile - to keep track of the renumbering\n");
#ifdef FOR_DL
  fprintf (stderr, "\n    Other command line options: \n");
#endif
  fprintf (stderr, "       -keep    nnn - prune if excess > nnn (default 160)\n");
  fprintf (stderr, "       -minpa   nnn - purge alg. primes >= nnn (default alim)\n");
  fprintf (stderr, "       -minpr   nnn - purge rat. primes >= nnn (default rlim)\n");
  fprintf (stderr, "       -nprimes nnn - expected number of prime ideals\n");
  fprintf (stderr, "       -raw         - output relations in CADO format\n");
  fprintf (stderr, "       -npthr   nnn - threads number for suppress singletons\n");
  fprintf (stderr, "       -inprel  file_rel_used : load active relations\n");
  fprintf (stderr, "       -outrel  file_rel_used : write active relations\n");
  fprintf (stderr, "       -npass   nnn - number of step of clique removal (default %d)\n", DEFAULT_NPASS);
  fprintf (stderr, "       -required_excess nnn - percentage of excess required at the end of the first singleton removal step (default %.2f)\n",
  DEFAULT_REQUIRED_EXCESS);
  exit (1);
}

int
main (int argc, char **argv)
{
  char *argv0 = argv[0];

  int k;
  size_t rel_used_nb_bytes;
  param_list pl;
  unsigned int npass = DEFAULT_NPASS;
  double required_excess = DEFAULT_REQUIRED_EXCESS;

  set_rep_cado(argv[0], rep_cado);
  wct0 = wct_seconds ();
  fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
  for (k = 1; k < argc; k++)
    fprintf (stderr, " %s", argv[k]);
  fprintf (stderr, "\n");

  param_list_init(pl);

  param_list_configure_switch(pl, "raw", &raw);

  argv++,argc--;

  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) continue;
    /* Since we accept file names freeform, we decide to never abort
     * on unrecognized options */
    if (!strcmp(*argv, "--help")) usage(argv0);
    break;
  }

  param_list_parse_uint64(pl, "nrels", (uint64_t *) &nrelmax);
  param_list_parse_uint64(pl, "nprimes", (uint64_t *) &nprimemax);
  param_list_parse_uint64(pl, "keep", (uint64_t *) &keep);

  param_list_parse_uint64(pl, "minindex", (uint64_t *) &min_index);

  /* param_list_parse_uint(pl, "npthr", (unsigned int *) &npt); */
  const char * snpt = param_list_lookup_string(pl, "npthr");
  if (snpt) {
    char *p, oldp;
    if ((p = strchr(snpt, 'x'))) {
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
    } else
      if (!sscanf(snpt, "%u", &npt))
        {
          fprintf (stderr, "Malformed -npthr option: %s\n", snpt);
          usage(argv0);
        }
  }
  param_list_parse_uint(pl, "npass", &npass);
  param_list_parse_double(pl, "required_excess", &required_excess);
  const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
  const char * filelist = param_list_lookup_string(pl, "filelist");
  const char * basepath = param_list_lookup_string(pl, "basepath");
  const char * subdirlist = param_list_lookup_string(pl, "subdirlist");
  const char * purgedname = param_list_lookup_string(pl, "out");
  const char * sos = param_list_lookup_string(pl, "sos");
  const char * infilerel = param_list_lookup_string(pl, "inrel");
  const char * outfilerel = param_list_lookup_string(pl, "outrel");
#ifdef FOR_DL
  const char * deletedname = param_list_lookup_string(pl, "outdel");
#endif

  search_antebuffer (rep_cado, path_antebuffer, antebuffer);

  binfilerel = (infilerel != 0);
  boutfilerel = (outfilerel != 0);

  cado_poly_init (pol);

  const char * tmp;

  ASSERT_ALWAYS((tmp = param_list_lookup_string(pl, "poly")) != NULL);
#ifndef FOR_FFS
  cado_poly_read(pol, tmp);
#else
  ffs_poly_read(pol, tmp);
#endif

  if (param_list_warn_unused(pl)) {
    fprintf (stderr, "Unused options in command-line\n");
    usage(argv0);
  }

  if ((basepath || subdirlist) && !filelist) {
    fprintf(stderr, "-basepath / -subdirlist only valid with -filelist\n");
    usage(argv0);
  }

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

  /* If nrels or nprimes > 2^32, then we need index_t to be 64-bit */
  if (((nprimemax >> 32) != 0 || (nrelmax >> 32) != 0) && sizeof(index_t) < 8)
  {
    fprintf (stderr, "Error, -nrels or -nprimes is too large for a 32-bit "
                     "program\nSee #define index_size in macros.h\n");
    exit(1);
  }

#ifdef FOR_DL /* For FFS we need to remember the renumbering of primes*/
  if (sos == NULL)
    {
      fprintf (stderr, "Error, missing -sos option.\n");
      exit(1);
    }
  if (deletedname == NULL)
    {
      fprintf (stderr, "Error, missing -outdel option.\n");
      exit(1);
    }
#endif

  fprintf (stderr, "Weight function used during clique removal:\n"
                   "  0     1     2     3     4     5     6     7\n"
                   "%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n",
                   weight_function_clique(0), weight_function_clique(1),
                   weight_function_clique(2), weight_function_clique(3),
                   weight_function_clique(4), weight_function_clique(5),
                   weight_function_clique(6), weight_function_clique(7));

  fprintf (stderr, "Number of relations is %lu\n", nrelmax);
  fprintf (stderr, "Number of prime ideals below the two large prime bounds: "
                   "%lu\n", nprimemax);

  SMALLOC(ideals_weight, nprimemax, "ideals_weight");
  tot_alloc_bytes = nprimemax * sizeof(weight_t);
  /* %zu is the C99 modifier for size_t */
  fprintf (stderr, "Allocated ideals_weight of %luMb (total %zuMb so far)\n",
                   nprimemax >> 20, tot_alloc_bytes >> 20);

  rel_used_nb_bytes = (nrelmax + 7) >> 3;
  bit_vector_init(rel_used, nrelmax);
  if (!rel_used->p)
  {
    fprintf (stderr, "Error, cannot allocate memory (%lu bytes) for "
                     "rel_used.\n", (unsigned long) rel_used_nb_bytes);
    exit (1);
  }

  tot_alloc_bytes += rel_used_nb_bytes;
  fprintf (stderr, "Allocated rel_used of %luMB (total %zuMB so far)\n",
                   rel_used_nb_bytes >> 20, tot_alloc_bytes >> 20);

  if (binfilerel)
  {
    FILE *in;
    void *p, *pl, *pf;
    int pipe;

    fprintf (stderr, "Loading rel_used file %s, %lu bytes\n", infilerel,
                     (unsigned long) rel_used_nb_bytes);
    if (!(in = fopen_maybe_compressed2 (infilerel, "r", &pipe, NULL)))
    {
      fprintf (stderr, "Error, cannot open file %s for reading.\n", infilerel);
      exit (1);
    }
    newnrel = 0;
    for (p = (void *) rel_used->p, pf = p + rel_used_nb_bytes; (pf - p); )
    {
      pl = p + fread(p, 1, pf - p, in);
      while (p != pl)
        newnrel += nbbits[*((uint8_t *) p++)];

      if ((p == pf) ^ !feof(in))
      {
        fprintf (stderr, "Error, the length of file %s is incorrect compared "
                         "with nrels (%lu).\n", infilerel, nrelmax);
        exit (1);
      }
    }
    fclose(in);
    fprintf (stderr, "Number of used relations is %lu\n",
                     (unsigned long) newnrel);
  }
  else
    bit_vector_set(rel_used, 1);

  /* What does it do ? */
  if (nrelmax & (BV_BITS - 1))
    rel_used->p[nrelmax>>LN2_BV_BITS] &= (((bv_t) 1)<<(nrelmax & (BV_BITS - 1))) - 1;

  if (!boutfilerel)
  {
    SMALLOC(rel_compact, nrelmax, "rel_compact");
    tot_alloc_bytes += nrelmax * (sizeof (index_t *));
  fprintf (stderr, "Allocated rel_compact of %luMB (total %zuMB so far)\n",
           (nrelmax * sizeof (index_t *)) >> 20, tot_alloc_bytes >> 20);
  }

    SMALLOC(sum, nprimemax, "sum");


  /* Build the file list (ugly). It is the concatenation of all
   *  b s p
   * where:
   *    b is the basepath (empty if not given)
   *    s ranges over all subdirs listed in the subdirlist (empty if no
   *    such list)
   *    p ranges over all paths listed in the filelist.
   *
   * If files are provided directly on the command line, the basepath
   * and subdirlist arguments are ignored.
   */
  if (!filelist) {
    fic = argv;
  } else if (!subdirlist) {
    fic = filelist_from_file (basepath, filelist, 0);
  } else {
    /* count the number of files in the filelist */
    int nfiles = 0;
    int nsubdirs = 0;
    char ** fl = filelist_from_file (NULL, filelist, 0);
    for(char ** p = fl ; *p ; p++, nfiles++);

    char ** sl = filelist_from_file (basepath, subdirlist, 1);
    for(char ** p = sl ; *p ; p++, nsubdirs++);

    SMALLOC(fic, nsubdirs * nfiles + 1, "main 3");
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
  }

  nrel = nrelmax;

  /************************** first pass *************************************/

  fprintf (stderr, "Pass 1, filtering ideals with index h >= %lu\n",
                   (unsigned long) min_index);

  prempt_scan_relations_pass_one ();

  if (!binfilerel)
    fprintf (stderr, "   nrels=%lu, nprimes=%lu; excess=%ld\n",
       (unsigned long) nrel, (unsigned long) nprimes, ((long) nrel) - nprimes);
  else
    fprintf (stderr, "   nrels=%lu, nrels used=%lu, nprimes=%lu; excess=%ld\n",
       (unsigned long) nrel, (unsigned long) newnrel, (unsigned long) nprimes, ((long) newnrel) - nprimes);

  if (!boutfilerel) {
    remove_singletons (npass, required_excess);
    fprintf (stderr, "   nrel=%lu, nprimes=%lu; excess=%ld\n",
       (unsigned long) nrel, (unsigned long) nprimes, ((long) nrel) - nprimes);
    if (nrel <= nprimes) /* covers case nrel = nprimes = 0 */
      {
  fprintf(stderr, "number of relations <= number of ideals\n");
  exit (2);
      }
  }

  if (!boutfilerel) {
    my_malloc_free_all ();

    fprintf (stderr, "Freeing rel_compact array...\n");
    /* we do not use it anymore */
    free (rel_compact);

    /*************************** second pass ***********************************/

    /* we renumber the primes in order of apparition in the hashtable */
    renumber (sos);

    /* reread the relation files and convert them to the new coding */
    fprintf (stderr, "Storing remaining relations...\n");
  }
  relsup = 0;
  prisup = 0;

#ifndef FOR_DL
  /* reread (purgedname, fic, rel_used, nrel_new, nprimes_new, raw); */
  prempt_scan_relations_pass_two (purgedname, rel_used, nrel, nprimes, raw);
#else
  /* reread (purgedname, deletedname, fic, rel_used, nrel_new, nprimes_new,                                                                raw); */
  prempt_scan_relations_pass_two (purgedname, deletedname, nrel, nprimes, raw);
#endif

  if (boutfilerel) {
    FILE *out;
    int pipe;
    void *p, *pf;

    if (nrelmax & (BV_BITS - 1))
      rel_used->p[nrelmax>>LN2_BV_BITS] &= (((bv_t) 1)<<(nrelmax & (BV_BITS - 1))) - 1;

    fprintf (stderr, "Relations with at least one singleton found and suppress:%lu\nNumber of primes suppress : %lu\nWriting rel_used file %s, %lu bytes\n",
       (unsigned long) relsup, (unsigned long) prisup, outfilerel, (unsigned long) rel_used_nb_bytes);
    if (!(out = fopen_maybe_compressed2 (outfilerel, "w", &pipe, NULL))) {
      fprintf (stderr, "Purge main: rel_used file %s cannot be written.\n", outfilerel);
      exit (1);
    }
    for (p = (void *) rel_used->p, pf = p + rel_used_nb_bytes; p != pf;
   p += fwrite(p, 1, pf - p, out));
    fclose (out);
  }

  bit_vector_clear(rel_used);
  SFREE(sum);
  SFREE(newindex);
  cado_poly_clear (pol);

  if (filelist) filelist_clear(fic);

  param_list_clear(pl);

  print_timing_and_memory (wct0);

  return 0;
}
