/* purge --- remove singletons

Copyright 2008, 2009, 2010, 2011, 2012 Francois Morain, Paul Zimmermann

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
  This program works in several passes over the relation files (currently 3):
  - the first pass considers only rational primes >= minpr and algebraic
    primes >= minpa. This means that only those primes are loaded in memory,
    and only relations with singletons in that range can be removed. By
    default minpr and minpa are taken as rlim and alim respectively.
    Usually more than half of the relations are discarded in that pass.
  - the second pass takes minpr=minpa=0, i.e., all primes are now considered
    (but only for relations surviving the first pass). It removes the remaining
    relations having singletons (typically very few).
    Also, if the excess exceeds 'keep', it prunes the heavier relations.
  - the final pass goes through the relations again, and dumps the remaining
    ones in the format needed by 'merge'.

  Remark: more passes can be implemented to decrease the memory usage in the
  first one. Simply increase the initial minpr/minpa, and decrease them at
  each pass. The only requirement is that the pass before the final one
  must have minpr=minpa=0, to ensure no singleton remains for 'merge'.

  This program uses the following data structures:
  rel_used[i]    - non-zero iff relation i is kept (so far)
  rel_compact[i] - list of 'h' indices in H table of considered (p,r) for row i
                   (terminated by a -1 sentinel)
  h = getHashAddr (H, p, r) - index of prime ideal (p, r) in hash table H
                              (rational primes use r = -2)
  GET_HASH_P(H,h) - prime corresponding to index h
  GET_HASH_R(H,h) - root  corresponding to index h (-2 for rational prime)
  H->hashcount[h] - number of occurrences of (p, r) in current relations
                   (in the last pass, H->hashcount[] stores the index of the
                    prime ideal).
*/

#include "cado.h"
#include <gmp.h>
#include "mod_ul.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>

#include "utils.h"
#include "filter_matrix.h" /* for BURIED_MAX_DENSITY */

#define MAX_FILES 1000000
#define MAX_STEPS 10   /* maximal number of steps in each pass */
#define FINAL_BOUND 0  /* final bound for rational and algebraic ideals */

/* Main variables */
static hashtable_t H;
static int **rel_compact = NULL;
static int ret MAYBE_UNUSED;
static int nrel, nprimes = 0;
static unsigned long nrelmax = 0;
static int nrel_new, nprimes_new, Hsize, Hsizer, Hsizea;
static double excess = 1.01;    /* minimum initial excess */
static long keep = 160;    /* maximum final excess */
static long minpr = -1, minpa = -1; /* negative values mean use minpr=rlim and
				       minpa=alim in first pass */
static cado_poly pol;
static unsigned long tot_alloc, tot_alloc0;
static int need64 = 0; /* non-zero if large primes are > 2^32 */
static int final, pass = 0;
static int raw = 0;
static char ** fic;
static int noclique = 0;
static double wct0;
static bit_vector rel_used;
static relation_stream rs;
/* Size of relations buffer between parsing & insertRelation.
   VERY IMPORTANT. Seems best for 1<<7; try 1<<8 for possible
   light speedup if your CPU is faster than a Nehalem 3.6 Ghz.
   1<<7 -> about 10MB of buffering (only ~2% is really used) 
*/
#define T_REL (1<<7) 
			
static volatile unsigned long cpt_rel_a;
static volatile unsigned long cpt_rel_b;
static volatile unsigned int end_insertRelation = 0;
typedef struct {
  relation_t *rel;
  unsigned long *num_rel;
  unsigned int *ltmp_rel;
} __buf_rel_t;
static __buf_rel_t buf_rel[1];


/* I suppose the max speed is about 120MB/s. So to load 1 byte, I need 8ns.
   The block in file loading is PREMPT_ONE_READ, the block in others is
   the average sentence, about ~100 bytes -> ~800 ns.
   BUT it seems the minimal useful time for nanosleep is about 50-80
   micro seconds... with nanosleep (0,1<<13).
   With PREMPT_ONE_READ = 64K, PREMPT_ONE_READ<<3 = 512K, so the delai is
   about 500 µs. OK for nanosleep(0,PREMPT_ONE_READ<<3).
   For others, if I use nanosleep to keep minimal CPU, I have to wait at
   least 80 µs -> 100 sentences.
   So T_REL is 1<<7 to avoid an empty buffer!
   1<<6 is possible, because 120MB/s is optimistic.
*/
static const struct timespec 
wait_load = { 0, (PREMPT_ONE_READ<<3) }, wait_classical = { 0, 1<<13 };

/* Be careful. 1<<13 = 8µs; in fact, about 50-80 µs */
inline void attente_minimale_passive ()
{
  static const struct timespec wait_min = { 0, 1<<13 };
  nanosleep(&wait_min, NULL);
}

inline void attente_minimale_active ()
{
  /* about 1 (for a nehalem 3 Ghz) to 5 µs */
  unsigned int i = (1<<6);
  while (i--) 
    __asm__ volatile ( "\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
" : : );
}

/********************** own memory allocation routines ***********************/

/* Rationale: calling one malloc() for each relation in insertNormalRelation()
   is expensive, since malloc() allocates some extra information to keep track
   of every memory blocks. Instead, we allocate memory in big blocks of size
   BLOCK_SIZE. */

#define BLOCK_SIZE (1<<20)  /* memory blocks are allocated of that # of int's */
/* relcompact_list is a list of blocks, each one of BLOCK_SIZE ints */
static int **relcompact_list = NULL;
static unsigned long relcompact_used = BLOCK_SIZE; /* usage of current block */
static unsigned int relcompact_size = 0;  /* minimal size of relcompact_list */
static int relcompact_current = -1; /* index of current block */

/* return a pointer to an array of n ints */
static int*
my_malloc_int (unsigned long n)
{
  int *ptr;
  
  if (relcompact_used + n > BLOCK_SIZE) {
    relcompact_used = 0;
    if (((unsigned int) (++relcompact_current)) == relcompact_size) {
      relcompact_size = ((relcompact_size) ? (relcompact_size << 1) : (1<<16));
      if ((relcompact_list = (int **) realloc(relcompact_list, relcompact_size * sizeof(int *))) == NULL) {
	  fprintf (stderr, "my_malloc_int: realloc error : %s\n", strerror (errno));
	  exit (1);
	}
      }
    if ((relcompact_list[relcompact_current] = (int *) malloc (BLOCK_SIZE * sizeof (int))) == NULL) {
      fprintf (stderr, "my_malloc_int: malloc error : %s\n", strerror (errno));
      exit (1);
    }
  }
  ptr = relcompact_list[relcompact_current] + relcompact_used;
  relcompact_used += n;
  tot_alloc += n * sizeof (int);
  return ptr;
}

static void
my_malloc_free_all (void)
{
  while (relcompact_current >= 0)
    free (relcompact_list[relcompact_current--]);
  free (relcompact_list);
  relcompact_list = NULL;
  relcompact_used = BLOCK_SIZE;
  relcompact_size = 0;
}

/*****************************************************************************/

/* dirty trick to distinguish rational primes: we store -2 for their root */
static unsigned long minus2;

/* Adds in table_ind[] the indices of the rational primes in 'rel'.
   All primes are assumed to be different (and to appear with an odd exponent).
   nb_coeff is the current index in table_ind[] (input and output).
 */
static void
fprint_rat (int *table_ind, int *nb_coeff, relation_t *rel, hashtable_t *H)
{
    int i, nbc = *nb_coeff;

    for (i = 0; i < rel->nb_rp; i++)
      table_ind[nbc++] = H->hashcount[getHashAddr (H, rel->rp[i].p, minus2)];
    *nb_coeff = nbc;
}

/* Adds in table_ind[] the indices of the algebraic primes in 'rel'.
   All primes are assumed to be different (and to appear with an odd exponent).
   nb_coeff is the current index in table_ind[] (input and output).
 */
static void
fprint_alg (int *table_ind, int *nb_coeff, relation_t rel, hashtable_t *H)
{
    int i, nbc = *nb_coeff;

    for (i = 0; i < rel.nb_ap; i++)
      table_ind[nbc++] = H->hashcount[getHashAddr (H, rel.ap[i].p,
                                                      rel.ap[i].r)];
    *nb_coeff = nbc;
}

/* Adds a free relation in table_ind[]: a is the corresponding prime */
static void
fprint_free (int *table_ind, int *nb_coeff, relation_t rel, hashtable_t *H)
{
    long p = rel.a;
    int i, nbc = *nb_coeff, index;

    index = getHashAddr (H, p, minus2);
    ASSERT(H->hashcount[index] >= 0);
    table_ind[nbc++] = H->hashcount[index];
    for(i = 0; i < rel.nb_ap; i++)
      {
	index = getHashAddr (H, p, rel.ap[i].p);
        ASSERT_ALWAYS(H->hashcount[index] >= 0);
        table_ind[nbc++] = H->hashcount[index];
    }
    *nb_coeff = nbc;
}

/* Print the relation 'rel' in matrix format, i.e., a line of the form:

   i a b k t_1 t_2 ... t_k

   i (decimal) is the row index from the nodup file (starting at 0)
   a (signed decimal) is a
   b (signed decimal) is b
   k (decimal) is the number of rational and algebraic primes in the relation
   t_1 ... t_k (hexadecimal) are the indices of the primes (starting at 0)

   Return the weight of the relation.

   Assumes the remaining primes in 'rel' are those with an odd exponent,
   and are all different.

   WARNING: the primes in the input relation are not necessarily sorted.
*/
static int
fprint_rel_row (FILE *file, int irel, relation_t rel, hashtable_t *H)
{
  int i;
  int *table_ind;
  int nb_coeff;

  table_ind = (int*) malloc ((rel.nb_rp + rel.nb_ap) * sizeof (int));

  nb_coeff = 0;

  if (rel.b == 0) /* free relation */
      fprint_free (table_ind, &nb_coeff, rel, H);
  else
    {
      /* adds rational primes in table_ind */
      fprint_rat (table_ind, &nb_coeff, &rel, H);

      /* adds algebraic primes in table_ind */
      fprint_alg (table_ind, &nb_coeff, rel, H);
    }

  fprintf (file, "%d %"PRId64" %"PRIu64" %d",
          irel, rel.a, rel.b, nb_coeff);
  for (i = 0; i < nb_coeff; ++i)
    /* due to the +1 in renumber */
    fprintf (file, " " PURGE_INT_FORMAT, table_ind[i] - 1);
  fprintf (file, "\n");

  free (table_ind);

  return nb_coeff;
}

/* First we count the number of large primes; then we store all primes in
   the hash table, but not in the relations. This might end up with singletons
   here and there, but we don't care, since they will be dealt with in
   merge.

   Meaning of the different parameters:
   minpr, minpa: only ideals > minpr (resp. minpa) are considered on the
                 rational (resp. algebraic) side. This means that the output
                 might contain ideals <= minpr or minpa appearing only once.
*/
static inline void
insertNormalRelation (relation_t *rel, unsigned long num_rel, int *tmp)
{
  int itmp, j, h, ok;

#if 0
    reduce_exponents_mod2 (rel);

    /* first count number of "large" primes */
    for (j = 0; j < rel->nb_rp; j++)
      ltmp += ((long) rel->rp[j].p >= minpr); /* only consider primes >= minpr */
    for (j = 0; j < rel->nb_ap; j++)
      ltmp += ((long) rel->ap[j].p >= minpa); /* only consider primes >= minpa */

    /* ltmp is the number of considered primes in the relation.
       We might have ltmp=0 if all primes are less than minpr, maxpr. */
    tmp = my_malloc_int (ltmp + 1);
#endif
    itmp = 0; /* number of entries in tmp */

    /* trick: we use the same hash table for rational and algebraic
       primes, but use a fake root -2 for rational primes, which
       ensures there is no collision with algebraic primes */
    for (j = 0; j < rel->nb_rp; j++)
      {
        ok = (((long) rel->rp[j].p) >= minpr);
        if (ok || final)
          /* in the final pass, we need to insert all ideals, since we need
             to renumber them in the output file */
          h = hashInsert (&H, rel->rp[j].p, minus2);
        if (ok)
          {
            nprimes += (H.hashcount[h] == 1); /* new prime */
            tmp[itmp++] = h;
          }
      }

    for (j = 0; j < rel->nb_ap; j++)
      {
        ok = (((long) rel->ap[j].p) >= minpa);
        if (ok || final)
          {
            rel->ap[j].r = findroot (rel->a, rel->b, rel->ap[j].p);
            h = hashInsert (&H, rel->ap[j].p, rel->ap[j].r);
          }
        if (ok)
          {
            nprimes += (H.hashcount[h] == 1); /* new prime */
            tmp[itmp++] = h;
          }
      }

    tmp[itmp] = -1; /* sentinel */
    rel_compact[num_rel] = tmp;
}

/* The information is stored in the ap[].p part, which is odd, but convenient.
   rel->ap.p[0..deg[ contains the deg roots of f(x) mod p, leading to deg
   ideals for the factorization of (p). */
static inline void
insertFreeRelation (relation_t *rel, unsigned long num_rel, int *tmp)
{
  long p = (long) rel->a; /* rel->b == 0 */
    int j, h, itmp;

    /* the prime on the rational side is rel->a
       the prime ideal on the algebraic side are (rel->a, rel->ap[j].p) */
#if 0
    ltmp = 1 + (p >= minpr);
    if (p >= minpa)
      ltmp += rel->nb_ap;
    tmp = my_malloc_int (ltmp);
#endif
    itmp = 0;
    if ((p >= minpr) || final)
      h = hashInsert (&H, p, minus2);
    if (p >= minpr)
      {
        nprimes += (H.hashcount[h] == 1); /* new prime */
        tmp[itmp++] = h;
      }
    for (j = 0; j < rel->nb_ap; j++)
      {
        if ((p >= minpa) || final)
          h = hashInsert(&H, p, rel->ap[j].p);
        if (p >= minpa)
          {
            nprimes += (H.hashcount[h] == 1); /* new ideal */
            tmp[itmp++] = h;
          }
      }
    tmp[itmp++] = -1;
    rel_compact[num_rel] = tmp;
}

/* Read all relations from file, and fills the rel_used and rel_compact arrays
   for each relation i:
   - rel_used[i] = 0 if the relation i is deleted
     rel_used[i] = 1 if the relation i is kept (so far)
   - rel_compact is an array, terminated by -1, of pointers to the entries
     in the hash table for the considered primes

     Trick: we only read relations for which rel_used[i]=1.
 */
#if 0
static int
scan_relations (char **ficname, int *nprimes, hashtable_t *H,
                bit_vector_srcptr rel_used, int nrelmax, int **rel_compact,
		long minpr, long minpa, unsigned long *tot_alloc, int final)
{
    *nprimes = 0;

    ASSERT(rel_compact != NULL);

    relation_stream_init(rs);

    for( ; *ficname ; ficname++) {
	fprintf(stderr, "   %-70s\n", *ficname);
        relation_stream_openfile(rs, *ficname);
        for ( ; ; ) {
            int irel = rs->nrels;
            if (bit_vector_getbit (rel_used, irel) == 0) /* skipped relation */
              {
                if (relation_stream_get_skip (rs) < 0)
                  break; /* end of file */
              }
            else
              {
                if (relation_stream_get (rs, NULL, 0) < 0)
                  break;
                ASSERT_ALWAYS(rs->nrels <= nrelmax);
                if (rs->rel.b > 0)
                  insertNormalRelation (rel_compact, irel, nprimes,
                                H, &(rs->rel), minpr, minpa, tot_alloc, final);
                else
                  insertFreeRelation (rel_compact, irel, nprimes,
                                H, &(rs->rel), minpr, minpa, tot_alloc, final);
              }
            if (!relation_stream_disp_progress_now_p(rs))
              continue;

            fprintf(stderr,
                    "read %d relations in %.1fs"
                    " -- %.1f MB/s -- %.1f rels/s\n",
                    rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
        }
        relation_stream_closefile(rs);
    }
    fprintf(stderr,
            "read %d relations in %.1fs -- %.1f MB/s -- %.1f rels/s\n",
            rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
    if (rs->nrels != nrelmax) {
        fprintf (stderr, "Error, -nrels value should match the number of scanned relations\n");
        fprintf (stderr, "expected %d relations, found %d\n", nrelmax,
                 rs->nrels);
        exit (EXIT_FAILURE);
    }

    relation_stream_clear(rs);

    return 1;
}
#endif

/* Return a non-zero value iff some prime (ideal) in the array tab[] is single
   (tab[j] is the index in the hash table of the corresponding prime ideal).
*/
static int
has_singleton (int *tab, hashtable_t *H)
{
    int j;

    for (j = 0; tab[j] != -1; j++)
      if (H->hashcount[tab[j]] == 1)
        return 1;
    return 0;
}

/* Delete a relation: set rel_used[i] to 0, update the count of primes
   in that relation, and set rel_compact[i] to NULL.
   Warning: we only update the count of primes that we consider, i.e.,
   rational primes >= minpr and algebraic primes >= minpa.
*/
static void
delete_relation (int i, int *nprimes, hashtable_t *H,
                 bit_vector_ptr rel_used, int **rel_compact)
{
  int j, *tab = rel_compact[i];

  for (j = 0; tab[j] != -1; j++)
    {
      H->hashcount[tab[j]]--; /* remove one occurrence of prime 'j' */
      *nprimes -= (H->hashcount[tab[j]] == 0);
    }
  rel_compact[i] = NULL;
  bit_vector_clearbit(rel_used, i);
}

/* New pruning code, which optimizes the decrease of N*W where N is the number
   of rows, and W is the total weight. We consider the connected
   components of the relation R(i1,i2) iff i1 and i2 share a prime
   of weight 2. If we remove one component of n rows and total weight w,
   then we save w*N+n*W (neglecting 2nd order terms), thus we remove
   first the components with the largest value of n/N + w/W. */

/* Define the weight of a connected component: those with the largest cost
   are removed first in the pruning; w (unsigned long) is the component weight,
   W (double) is the total matrix weight, n (int) is the number of rows of the
   component, and N (double) is the number of rows of the matrix. This macro
   must return a double. */
#define COST_MODEL 0

#if COST_MODEL == 0
/* optimize the decrease of W*N */
#define COST(w,W,n,N) ((double) (w) / (W) + (double) (n) / (N))
#elif COST_MODEL == 1
/* optimize the decrease of W */
#define COST(w,W,n,N) ((double) (w) / (W))
#elif COST_MODEL == 2
/* optimize the decrease of N */
#define COST(w,W,n,N) ((double) (n) / (N))
#else
#error "Invalid cost model"
#endif

static int
compare (const void *v1, const void *v2)
{
    float w1 = *((float*) v1);
    float w2 = *((float*) v2);

    return (w1 >= w2) ? -1 : 1;
}

/* Compute connected component of row i for the relation R(i1,i2) if rows
   i1 and i2 share a prime of weight 2.
   Return number of rows of components, and put in w the total weight. */
static int
compute_connected_component (bit_vector_ptr T, uint32_t i, hashtable_t *H,
                             unsigned long *w, int **rel_compact,
                             uint32_t *sum, int *rel_weight)
{
  int j, h, n;
  uint32_t k;

  n = 1;         /* current row */
  bit_vector_setbit(T, i); /* mark row as visited */
  for (j = 0; (h = rel_compact[i][j]) != -1; j++)
    if (H->hashcount[h] == 2)
      {
        k = sum[h] - i; /* other row where prime of index h appears */
        if (!bit_vector_getbit(T, k)) /* row k was not visited yet */
          n += compute_connected_component (T, k, H, w, rel_compact, sum,
                                            rel_weight);
      }
  *w += rel_weight[i]; /* add weight of non-buried ideals in current row */
  return n;
}

/* Delete connected component of row i, assuming the bit-vector is set.
   Warning: we might have some H->hashcount[h] = 3, which is decreased
   to 2, but we don't want to treat that case. Thus we check in addition
   that sum[h] <> 0, which only occurs when H->hashcount[h] = 2 initially. */
static int
delete_connected_component (bit_vector_ptr T, uint32_t i, hashtable_t *H,
                            int **rel_compact, uint32_t *sum, int *nprimes,
                            bit_vector_ptr rel_used)
{
  int j, h, w = 1;
  uint32_t k;

  bit_vector_clearbit(T, i); /* mark row as visited */
  bit_vector_clearbit(rel_used, i);
  for (j = 0; (h = rel_compact[i][j]) != -1; j++)
    {
      if (H->hashcount[h] == 2 && sum[h] != 0)
        { /* first row that contains ideal of index h */
          k = sum[h] - i; /* other row where prime of index h appears */
          if (bit_vector_getbit(T, k) == 1) /* row k was not visited yet */
            w += delete_connected_component (T, k, H, rel_compact, sum,
                                             nprimes, rel_used);
        }
    }
  delete_relation (i, nprimes, H, rel_used, rel_compact);
  return w;
}

typedef struct {
  float w;
  uint32_t i;
} comp_t;

static void
deleteHeavierRows (hashtable_t *H, int *nrel, int *nprimes,
                   bit_vector_ptr rel_used, int **rel_compact, int nrelmax,
                   int keep, int *rel_weight)
{
  uint32_t *sum; /* sum of row indices for primes with weight 2 */
  int i, j, h, n, ltmp = 0;
  unsigned long w;
  double W = 0.0; /* total matrix weight */
  double N = 0.0; /* numebr of rows */
  comp_t *tmp = NULL; /* (weight, index) */
  int target;
  static int count = 0;

  if ((*nrel - *nprimes) <= keep)
    return;

  /* first collect sums for primes with weight 2, and compute total weight */
  sum = (uint32_t*) malloc (H->hashmod * sizeof (uint32_t));
  memset (sum, 0, H->hashmod * sizeof (uint32_t));
  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, i)) {
        for (j = 0; (h = rel_compact[i][j]) != -1; j++)
          if (H->hashcount[h] == 2)
            sum[h] += i;
        N += 1.0;
        W += rel_weight[i]; /* weight of non-buried ideals */
      }
  fprintf (stderr, "Matrix has %1.0f rows and weight %1.0f\n", N, W);
  ASSERT_ALWAYS(N == (double) *nrel);

  /* now initialize bit table for relations used */
  bit_vector T;
  bit_vector_init_set(T, nrelmax, 0);

  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, i) && bit_vector_getbit(T, i) == 0)
      {
        w = 0;
        n = compute_connected_component (T, i, H, &w, rel_compact, sum,
                                         rel_weight);
        ltmp ++;
        tmp = (comp_t*) realloc (tmp, ltmp * sizeof (comp_t));
        tmp[ltmp - 1].w = (float) COST(w,W,n,N);
        tmp[ltmp - 1].i = i;
      }

  qsort (tmp, ltmp, sizeof(comp_t), compare);

  /* remove heaviest components, assuming each one decreases the excess by 1;
     we remove only half of the excess at each call of deleteHeavierRows,
     hoping to get "better" components to remove at the next call. */
  target = (*nrel - *nprimes + keep) / 2;
  if (++count >= MAX_STEPS)
    target = keep; /* enough steps */
  for (i = 0; i < ltmp && (*nrel) - (*nprimes) > target; i ++)
    {
      /* we could remove only one row from the current connected component,
         the other would be removed by onepass_singleton_removal(),
         but this would take more passes, and would thus be less efficient */
      *nrel -= delete_connected_component (T, tmp[i].i, H, rel_compact, sum,
                                           nprimes, rel_used);
    }

  bit_vector_clear(T);
  free (sum);
  free (tmp);
}

static void
onepass_singleton_removal (int nrelmax, int *nrel, int *nprimes,
                           hashtable_t *H, bit_vector_ptr rel_used, int **rel_compact)
{
  int i;

  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, i) && has_singleton (rel_compact[i], H))
      {
        delete_relation (i, nprimes, H, rel_used, rel_compact);
        *nrel -= 1;
      }
}

static void
remove_singletons (int *nrel, int nrelmax, int *nprimes, hashtable_t *H,
                   bit_vector_ptr rel_used, int **rel_compact, int keep,
                   int final)
{
  int old, newnrel = *nrel, newnprimes = *nprimes, oldexcess, excess;
  int count = 0;
  int *rel_weight = NULL;

  if (final) /* we compute the weight of non-buried columns */
    {
      int i, j, h, w = 0, wsmall = 0;
      rel_weight = (int*) malloc (nrelmax * sizeof (int));
      for (i = 0; i < nrelmax; i++)
        if (bit_vector_getbit (rel_used, i) == 1) /* active relation */
          {
            rel_weight[i] = 0;
            for (j = 0; (h = rel_compact[i][j]) != -1; j++)
              {
                w ++;
                if ((double) H->hashcount[h] / (double) nrelmax < 1.0 / 
                    (double) BURIED_MAX_DENSITY)
                  rel_weight[i] ++;
              }
            wsmall += rel_weight[i];
          }
    }

  excess = newnrel - newnprimes;
  do
    {
      old = newnrel;
      oldexcess = excess;

      /* for the final pass, delete heavy rows */
      if (final)
        deleteHeavierRows (H, &newnrel, &newnprimes, rel_used, rel_compact,
                           nrelmax, keep, rel_weight);
      if (newnrel != old)
        fprintf (stderr, "deleted heavier relations: %d %d at %2.2lf\n",
                 newnrel, newnprimes, seconds ());

      onepass_singleton_removal (nrelmax, &newnrel, &newnprimes, H,
                                   rel_used, rel_compact);

      excess = newnrel - newnprimes;
      fprintf (stderr, "   new_nrows=%d new_ncols=%d (%d) at %2.2lf\n",
               newnrel, newnprimes, excess, seconds());

      if (final && oldexcess > excess)
        fprintf (stderr, "   [each excess row deleted %2.2lf rows]\n",
                 (double) (old - newnrel) / (double) (oldexcess - excess));

      count ++;
      if (final == 0 && count >= MAX_STEPS) /* for the final pass we want
                                               to remove all singletons
                                               to avoid warnings */
        break;
    }
  while (newnrel != old);

  /* Warning: we might have empty rows, i.e., empty rel_compact[i] lists,
     if all primes in a relation are less then minpr or minpa. */

  *nrel = newnrel;
  *nprimes = newnprimes;
  free (rel_weight);
}

/* This function renumbers used primes (those with H->hashcount[i] > 1)
   and puts the corresponding index in H->hashcount[i].

   At return, nprimes is the number of used primes.

   We locate used primes and do not try to do fancy things as sorting w.r.t.
   weight, since this will probably be done later on.
   All rows will be 1 more that needed -> subtract 1 in fprint...! */
static void
renumber (int *nprimes, hashtable_t *H, const char *sos)
{
    FILE *fsos = NULL;
    unsigned int i;
    int nb = 1; /* we start at 1 here, but subtract 1 in fprint_rel_row */

    if (sos != NULL)
      {
	fprintf (stderr, "Output renumber table in file %s\n", sos);
	fsos = gzip_open (sos, "w");
        fprintf (fsos, "# each row contains 3 hexadecimal values: i p r\n");
        fprintf (fsos, "# i is the ideal index value (starting from 0)\n");
        fprintf (fsos, "# p is the corresponding prime\n");
        fprintf (fsos, "# r is the corresponding root (-2=fffffffe on the rational side)\n");
      }
    for(i = 0; i < H->hashmod; i++)
      if (H->hashcount[i] == 0)
        H->hashcount[i] = -1;
      else
        {
          /* Since we consider only primes >= minpr or minpa,
             smaller primes might appear only once here, thus we can't
             assert H->hashcount[i] > 1, but H->hashcount[i] = 1 should
             be rare if minpr/minpa are well chosen (not too large). */
          static int count = 0;
          if (H->hashcount[i] == 1 && count ++ < 10)
            {
              if (GET_HASH_R(H,i) == minus2)
                fprintf (stderr, "Warning: singleton rational prime %"
                        PRIu64 "\n",
                         GET_HASH_P(H,i));
              else
                fprintf (stderr, "Warning: singleton algebraic ideal (%"
                        PRIu64",%"PRIu64")\n",
                         GET_HASH_P(H,i), GET_HASH_R(H,i));
            }
          H->hashcount[i] = nb++;
          if (fsos != NULL)
            fprintf(fsos, "%x %" PRIx64 " %" PRIx64 "\n",
                    H->hashcount[i] - 1, GET_HASH_P(H,i), GET_HASH_R(H,i));
	}
    if (fsos != NULL)
      gzip_close (fsos, sos);
    nb--;
    *nprimes = nb;
}

/* Read again the relation files ficname[0], ..., ficname[nbfic-1],
   and output remaining relations (those with rel_used[i] <> 0) in ofile.

   If raw is non-zero, output relations in CADO format
   (otherwise in format used by merge).
*/
static void
reread (const char *oname, char ** ficname, hashtable_t *H,
        bit_vector_srcptr rel_used, int nrows, int ncols, int raw)
{
  FILE *ofile;
  int ret MAYBE_UNUSED, nr = 0;
  double W = 0.0; /* total weight */
  int pipe;

  ofile = fopen_compressed_w(oname, &pipe, NULL);
  if (raw == 0)
    fprintf (ofile, "%d %d\n", nrows, ncols);
  fprintf (stderr, "Final pass:\n");

  relation_stream_init(rs);

  for ( ; *ficname ; ficname++) {
      relation_stream_openfile(rs, *ficname);
      fprintf(stderr, "   %-70s\n", *ficname);
      for ( ; ; ) {
          int irel = rs->nrels;
          if (bit_vector_getbit (rel_used, irel) == 0) /* skipped relation */
            {
              if (relation_stream_get_skip (rs) < 0)
                break; /* end of file */
            }
          else
            {
              if (relation_stream_get(rs, NULL, 0) < 0)
                break;
              // ASSERT_ALWAYS(rs->nrels <= nrelmax);
              if (raw == 0)
              {
                  if (rs->rel.b > 0)
                  {
                      reduce_exponents_mod2 (&rs->rel);
                      computeroots (&rs->rel);
                  }
                  W += (double) fprint_rel_row (ofile, irel, rs->rel, H);
              }
              else
                  fprint_relation_raw (ofile, &rs->rel);
              nr++;
              if (nr >= nrows)
                  ret = 0; /* we are done */
            }
          if (!relation_stream_disp_progress_now_p(rs))
              continue;

          fprintf(stderr,
                  "re-read %lu relations in %.1fs"
                  " -- %.1f MB/s -- %.1f rels/s\n",
                  rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
      }
      relation_stream_closefile(rs);
  }
  relation_stream_trigger_disp_progress(rs);
  fprintf(stderr,
          "re-read %lu relations in %.1fs"
          " -- %.1f MB/s -- %.1f rels/s\n",
          rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
  relation_stream_clear(rs);
  if (pipe) pclose(ofile); else fclose(ofile);

  /* write excess to stdout */
  if (raw == 0)
    printf ("NROWS:%d WEIGHT:%1.0f WEIGHT*NROWS=%1.2e\n",
            nr, W, W * (double) nr);
  printf ("EXCESS: %d\n", nrows - ncols);
  fflush (stdout);
}

static void
prempt_load (prempt_t prempt_data) {
  char **p_files = prempt_data->files, *pprod;
  FILE *f;
  size_t try_load, load;
  char *p, *l, *pmax = &(prempt_data->buf[PREMPT_BUF]);

  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);

  while (*p_files) {
    if (!(f = popen (*p_files, "r"))) {
      fprintf (stderr, "prempt_load: popen error. %s\n", strerror (errno));
      exit (1);
    }
    p = strstr(*p_files, "/");
    if (p) {
      *p = 0;
      fprintf (stderr, "%s\n", *p_files);
      *p = '/';
      for ( ; ; ) {
	l = strstr(p, " ");
	if (l) {
	  *l = 0;
	  fprintf(stderr, "   %-70s\n", p);
	  *l = ' ';
	  p = strstr(&(l[1]), "/");
	}
	else {
	  fprintf(stderr, "   %-70s\n", p);
	  break;
	}
      }
    }
    pprod = (char *) prempt_data->pprod;
    for ( ; ; )
      if ((pprod != prempt_data->pcons) &&
	  (((((PREMPT_BUF + ((size_t) prempt_data->pcons))) - ((size_t) pprod)) &
	    (PREMPT_BUF - 1)) <= PREMPT_ONE_READ)) {
	nanosleep(&wait_load, NULL);
      }
      else {
	try_load = MIN((PREMPT_BUF + ((size_t)prempt_data->buf)) - ((size_t) pprod), PREMPT_ONE_READ);
	if ((load = fread (pprod, 1, try_load, f))) {
	  pprod += load;
	  if (pprod == pmax) pprod = prempt_data->buf;
	  prempt_data->pprod = pprod;
	}
	else
	  if (feof(f)) {
	    pclose(f);
	    free(*p_files);
	    *p_files++ = NULL;
	    break;
	  }
	  else {
	    fprintf (stderr, "prempt_load: load error (%s) from\n%s\n", strerror (errno), *p_files);
	    exit (1);
	  }
      }
  }
  prempt_data->end = 1;
  pthread_exit (NULL);
}

static inline
unsigned int relation_stream_get_fast (prempt_t prempt_data, relation_t *rel)
{
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
  int64_t n;
  char *p, *pmax, *pmin;
  unsigned int k, i;
  unsigned long pr;
  unsigned char c, v;
  unsigned int ltmp;
  
#define LOAD_ONE(P) { c = *(P)++; if ((P) == pmax) (P) = pmin; }
  
  ltmp = 0;

  pmin = prempt_data->buf;
  p = (char *) prempt_data->pcons;
  pmax = &(prempt_data->buf[PREMPT_BUF]);

  LOAD_ONE(p);
  if (c == '-') {
    rel->a = -1;
    LOAD_ONE(p);
  }
  else
    rel->a = 1;
  for (n = 0 ; (v = ugly[c]) < 10 ; ) {
    n = n * 10 + v;
    LOAD_ONE(p);
  }
  ASSERT_ALWAYS(c == ',');
  rel->a *= n;
  
  n = 0;
  LOAD_ONE(p);
  for ( ; (v = ugly[c]) < 10 ; ) {
    n = n * 10 + v;
    LOAD_ONE(p);
  }
  ASSERT_ALWAYS(c == ':');
  rel->b = n;
  
  if (!rs->parse_only_ab) {
    /* Do something if we're also interested in primes */
    ltmp = 1;

    for ( k = 0, c = 0 ; ; ) {
    next_rat:
      if (c == ':') break;
      LOAD_ONE(p);
      for (pr = 0 ; (v = ugly[c]) < 16 ; ) {
	pr <<= 4;
	pr += v;
	LOAD_ONE(p);
      }
      ASSERT_ALWAYS(c == ',' || c == ':');
      if (pr) {
	for (i = k; i--; ) {
	  if (rel->rp[i].p == pr) {
	    rel->rp[i].e++;
	    goto next_rat;
	  }
	}
	rel->rp[k++] = (rat_prime_t) { .p = pr,.e = 1};
      }
    }
    rel->nb_rp = k;
    
    for ( k = 0 ; ; ) {
    next_alg:
      if (c == '\n') break;
      LOAD_ONE(p);
      for (pr = 0 ; (v = ugly[c]) < 16 ; ) {
	pr <<= 4;
	pr += v;
	LOAD_ONE(p);
      }
      ASSERT_ALWAYS(c == ',' || c == '\n');
      if (pr) {
	for (i = k; i--; ) {
	  if (rel->ap[i].p == pr) {
	    rel->ap[i].e++;
	    goto next_alg;
	  }
	}
	rel->ap[k++] = (alg_prime_t) { .p = pr,.r = -1,.e = 1};
      }
    }
    rel->nb_ap = k;

    if (rel->b > 0) {
      reduce_exponents_mod2 (rel);
      for (i = rel->nb_rp; i ; ltmp += ((long) rel->rp[--i].p >= minpr));
      for (i = rel->nb_ap; i ; ltmp += ((long) rel->ap[--i].p >= minpa));
    }
    else {
      ltmp += (long) rel->a >= minpr;
      if ((long) rel->a >= minpa) ltmp += rel->nb_ap;
    }
  }
  return ltmp;
}

void insertRelation() {
  unsigned int n;
  unsigned long cpy_cpt_rel_b;
  int *tmp;

  /*
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);
  */
  cpy_cpt_rel_b = cpt_rel_b;
  while (cpt_rel_a == cpy_cpt_rel_b) nanosleep (&wait_classical, NULL);
  for ( ; ; ) {
    n = (unsigned int) (cpy_cpt_rel_b & (T_REL - 1));
    tmp = my_malloc_int(buf_rel->ltmp_rel[n]);
    if (buf_rel->rel[n].b > 0)
      insertNormalRelation (&(buf_rel->rel[n]), buf_rel->num_rel[n], tmp);
    else
      insertFreeRelation (&(buf_rel->rel[n]), buf_rel->num_rel[n], tmp);
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
    if (relation_stream_disp_progress_now_p(rs))
      fprintf(stderr,
	      "read useful %lu relations in %.1fs"
	      " -- %.1f MB/s -- %.1f rels/s\n",
	      rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
    while (cpt_rel_a == cpy_cpt_rel_b) {
      if (!end_insertRelation) nanosleep (&wait_classical, NULL);
      else if (cpt_rel_a == cpy_cpt_rel_b) pthread_exit(NULL);
    }
  }
}

static int
prempt_scan_relations () {
  char *pcons, *pcons_old, *pcons_max, *p;
  pthread_attr_t attr;
  pthread_t thread_load, thread_relation;
  prempt_t prempt_data;
  unsigned long cpy_cpt_rel_a;
  unsigned int dispo, length_line, i;
  int err;
  char c;
  static int pass1 = 1;

  end_insertRelation = 0;
  if (((buf_rel->num_rel  = malloc (sizeof(*(buf_rel->num_rel))  * T_REL)) == NULL) ||
      ((buf_rel->ltmp_rel = malloc (sizeof(*(buf_rel->ltmp_rel)) * T_REL)) == NULL) ||
      ((buf_rel->rel      = malloc (sizeof(*(buf_rel->rel))      * T_REL)) == NULL)) {
    fprintf (stderr, "prempt_scan_relations: malloc error. %s\n", strerror (errno));
    exit (1);
  }
  memset(buf_rel->rel, 0, sizeof(*(buf_rel->rel)) * T_REL);
  for (i = T_REL ; i ; relation_provision_for_primes(&(buf_rel->rel[--i]), RELATION_MAX_BYTES>>2, RELATION_MAX_BYTES>>2));
  cpt_rel_a = cpt_rel_b = 0;
  cpy_cpt_rel_a = cpt_rel_a;
  nprimes = 0;
  ASSERT(rel_compact != NULL);
  relation_stream_init(rs);
  rs->pipe = 1;
  length_line = 0;
  
  prempt_data->files = prempt_open_compressed_rs (fic);
  if (!(prempt_data->buf = malloc (PREMPT_BUF))) {
    fprintf (stderr, "prempt_scan_relations: malloc error. %s\n", strerror (errno));
    exit (1);
  }
  prempt_data->pcons = prempt_data->buf;
  prempt_data->pprod = prempt_data->buf;
  pcons_max = &(prempt_data->buf[PREMPT_BUF]);
  prempt_data->end = 0;
  pthread_attr_init (&attr);
  pthread_attr_setstacksize (&attr, 1<<16);
  pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
  if ((err = pthread_create (&thread_load, &attr, (void *) prempt_load, prempt_data))) {
    fprintf (stderr, "prempt_scan_relations: pthread_create erreur 1: %d. %s\n", err, strerror (errno)); 
    exit (1);
  }
  if ((err = pthread_create (&thread_relation, &attr, (void *) insertRelation, NULL))) {
    fprintf (stderr, "prempt_scan_relations: pthread_create erreur 2: %d. %s\n", err, strerror (errno)); 
    exit (1);
  }
  
  pcons = (char *) prempt_data->pcons;
  for ( ; ; ) {
    ASSERT_ALWAYS(length_line <= ((unsigned int) RELATION_MAX_BYTES));
    rs->pos += length_line;
    length_line = 0;
    prempt_data->pcons = pcons;

    while (pcons == prempt_data->pprod) {						
      if (prempt_data->end) if (pcons == prempt_data->pprod) goto end_of_files;
      nanosleep (&wait_classical, NULL);
    }
    if (*pcons == '#')                           goto goto_endline;
    if (pass1)                                   goto valid_line;
    if (bit_vector_getbit (rel_used, rs->nrels)) goto valid_line;
    rs->nrels++;

  goto_endline:
    do {
      while (pcons == prempt_data->pprod) {
	if (prempt_data->end)
	  if (pcons == prempt_data->pprod) {
	    fprintf (stderr, "prempt_scan_relations: at the end of files, a line without \\n ?\n");
	    exit (1); 
	  }
	nanosleep (&wait_classical, NULL);
      }
      p = ((pcons <= prempt_data->pprod) ? (char *) prempt_data->pprod : pcons_max) - 1;
      c = *p;
      *p = '\n';
      pcons_old = pcons;
      while (*pcons++ != '\n');
      *p = c;
      length_line += (pcons - pcons_old);
      err = (pcons > p && c != '\n');
      if (pcons == pcons_max) pcons = prempt_data->buf;
    } while (err);
    rs->lnum++;
    continue;

  valid_line:
    /* Here, the next line is useful. But we have to read a complete line (with \n). So,
       we try to load >= RELATION_MAX_BYTES. After, we search the number of caracters,
       and write at the MIN(nb_caracters, PREMPT_BUF) '\n' to speed up the research. */
    while ((dispo = ((PREMPT_BUF + ((size_t) prempt_data->pprod)) - ((size_t) pcons)) & (PREMPT_BUF - 1)) <= ((unsigned int) RELATION_MAX_BYTES)) {
      if (prempt_data->end) {
	dispo = ((PREMPT_BUF + ((size_t) prempt_data->pprod)) - ((size_t) pcons)) & (PREMPT_BUF - 1);
	break;
      }
      nanosleep(&wait_classical, NULL);
    }
    do {
      if (pcons == prempt_data->pprod) {
	fprintf (stderr, "prempt_scan_relations: relation line size is greater than RELATION_MAX_BYTES (%d)\n", RELATION_MAX_BYTES);
	exit(1);
      }
      p = ((pcons < prempt_data->pprod) ? (char *) prempt_data->pprod : pcons_max) - 1;
      c = *p;
      *p = '\n';
      pcons_old = pcons;
      while (*pcons++ != '\n');
      *p = c;
      length_line += (pcons - pcons_old);
      err = (pcons > p && c != '\n');
      if (pcons == pcons_max) pcons = prempt_data->buf;
    } while (err);
    /* Ok, we have a complete line */
    while (cpy_cpt_rel_a > cpt_rel_b+T_REL-1) {
      nanosleep(&wait_classical, NULL);
    }
    buf_rel->num_rel [cpy_cpt_rel_a & (T_REL - 1)] = rs->nrels;
    buf_rel->ltmp_rel[cpy_cpt_rel_a & (T_REL - 1)] =
      relation_stream_get_fast (prempt_data, &(buf_rel->rel[cpy_cpt_rel_a & (T_REL - 1)]));
    cpy_cpt_rel_a++;
    cpt_rel_a = cpy_cpt_rel_a;
    rs->nrels++;
    rs->lnum++;
  }
 end_of_files:
  while (cpy_cpt_rel_a != cpt_rel_b) nanosleep(&wait_classical, NULL);
  end_insertRelation = 1;
  pthread_join(thread_relation, NULL);
  if (pthread_tryjoin_np (thread_load, NULL))
    pthread_cancel(thread_load);
  pthread_join(thread_load, NULL);
  pthread_attr_destroy(&attr);
  free(prempt_data->buf); prempt_data->buf = NULL;
  
  fprintf (stderr, "read %lu relations in %.1fs -- %.1f MB/s -- %.1f rels/s\n",
           rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
  if (rs->nrels != nrelmax) {
    fprintf (stderr, "Error, -nrels value should match the number of scanned relations\nexpected %lu relations, found %lu\n", nrelmax, rs->nrels);
    exit (EXIT_FAILURE);
  }
  
  for (i = T_REL ; i ; free(buf_rel->rel[--i].rp));
  free (buf_rel->rel);      buf_rel->rel      = NULL;
  free (buf_rel->num_rel);  buf_rel->num_rel  = NULL;
  free (buf_rel->ltmp_rel); buf_rel->ltmp_rel = NULL;

  pass1 = 0;
  
  return 1;
}

static void
usage (void)
{
  fprintf (stderr, "Usage: purge [options] -poly polyfile -out purgedfile -nrels nnn [-basepath <path>] [-subdirlist <sl>] [-filelist <fl>] file1 ... filen\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "       -excess  nnn - initial excess (default 1.01)\n");
  fprintf (stderr, "       -keep    nnn - prune if excess > nnn (default 160)\n");
  fprintf (stderr, "       -minpa   nnn - purge alg. primes >= nnn (default alim)\n");
  fprintf (stderr, "       -minpr   nnn - purge rat. primes >= nnn (default rlim)\n");
  fprintf (stderr, "       -nprimes nnn - expected number of prime ideals\n");
  fprintf (stderr, "       -sos sosfile - to keep track of the renumbering\n");
  fprintf (stderr, "       -raw         - output relations in CADO format\n");
  fprintf (stderr, "       -noclique n  - if n <> 0, don't do clique removal\n");
  exit (1);
}

/* estimate the number of primes <= B */
static int
approx_phi (long B)
{
  ASSERT_ALWAYS((double) B <= 53030236260.0); /* otherwise B/log(B) > 2^31 */
  return (B <= 1) ? 0 : (int) ((double) B / log ((double) B));
}

int
main (int argc, char **argv)
{
  int k;
  
  wct0 = wct_seconds ();
  fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
  for (k = 1; k < argc; k++)
    fprintf (stderr, " %s", argv[k]);
  fprintf (stderr, "\n");
  
  param_list pl;
  param_list_init(pl);
  
  param_list_configure_knob(pl, "raw", &raw);
  
  argv++,argc--;
  
  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
    /* Since we accept file names freeform, we decide to never abort
     * on unrecognized options */
    if (strcmp(*argv, "--help") == 0)
      usage();
    break;
  }
  
  param_list_parse_ulong(pl, "nrels", &nrelmax);
  param_list_parse_int(pl, "nprimes", &nprimes);
  param_list_parse_long(pl, "minpr", &minpr);
  param_list_parse_long(pl, "minpa", &minpa);
  param_list_parse_double(pl, "excess", &excess);
  param_list_parse_long(pl, "keep", &keep);
  param_list_parse_int(pl, "noclique", &noclique);

  const char * filelist = param_list_lookup_string(pl, "filelist");
  const char * basepath = param_list_lookup_string(pl, "basepath");
  const char * subdirlist = param_list_lookup_string(pl, "subdirlist");
  const char * purgedname = param_list_lookup_string(pl, "out");
  const char * sos = param_list_lookup_string(pl, "sos");
  
  cado_poly_init (pol);
  
  const char * tmp;
  
  ASSERT_ALWAYS((tmp = param_list_lookup_string(pl, "poly")) != NULL);
  cado_poly_read(pol, tmp);
  
  if (param_list_warn_unused(pl)) {
    usage();
  }
  
  if ((basepath || subdirlist) && !filelist) {
    fprintf(stderr, "-basepath / -subdirlist only valid with -filelist\n");
    usage();
  }
  
  if (nrelmax == 0)
    {
      fprintf (stderr, "Error, missing -nrels ... option (or nrels=0)\n");
      usage ();
    }
  
  /* On a 32-bit computer, even 1 << 32 would overflow. Well, we could set
     map[ra] = 2^32-1 in that case, but not sure we want to support 32-bit
     primes on a 32-bit computer... */
  need64 = (pol->rat->lpb >= 32) || (pol->alg->lpb >= 32);
  
  if (need64 && sizeof (long) < 8)
    {
      fprintf (stderr, "Error, too large LPBs for a 32-bit computer\n");
      usage();
    }
  
  minus2 = (need64) ? 18446744073709551614UL : 4294967294UL;
  
  if (minpr < 0)
    minpr = pol->rat->lim;
  if (minpa < 0)
    minpa = pol->alg->lim;
  
  fprintf (stderr, "Number of relations is %lu\n", nrelmax);
  if (nprimes > 0)
    Hsize = nprimes;
  else{
    /* Estimating the number of needed primes (remember that hashInit
       multiplies by a factor 1.5). */
    Hsizer = approx_phi (1L << pol->rat->lpb) - approx_phi (minpr);
    Hsizea = approx_phi (1L << pol->alg->lpb) - approx_phi (minpa);
    Hsize = Hsizer + Hsizea;
  }
  fprintf (stderr, "Estimated number of prime ideals: %d\n", Hsize);
  hashInit (&H, Hsize, 1, need64);
  tot_alloc0 = H.hashmod * H.size;
  
  bit_vector_init_set(rel_used, nrelmax, 1);
  tot_alloc0 += nrelmax;
  fprintf (stderr, "Allocated rel_used of %luMb (total %luMb so far)\n",
	   nrelmax >> 20,
	   tot_alloc0 >> 20);
  
  /* FIXME: we could allocate rel_compact[] at each pass according to the
     number of *remaining* relations. We could also manage so that it takes
     only 32 bits per entry on a 64-bit machine. */
  rel_compact = (int **) malloc (nrelmax * sizeof (int *));
  tot_alloc0 += nrelmax * sizeof (int*);
  /* %zu is the C99 modifier for size_t */
  fprintf (stderr, "Allocated rel_compact of %zuMb (total %luMb so far)\n",
	   (nrelmax * sizeof (int *)) >> 20,
	   tot_alloc0 >> 20);
  
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
    fic = filelist_from_file(basepath, filelist);
  } else {
    /* count the number of files in the filelist */
    int nfiles = 0;
    int nsubdirs = 0;
    char ** fl = filelist_from_file(NULL, filelist);
    for(char ** p = fl ; *p ; p++, nfiles++);
    
    char ** sl = filelist_from_file(basepath, subdirlist);
    for(char ** p = sl ; *p ; p++, nsubdirs++);
    
    fic = malloc((nsubdirs * nfiles + 1) * sizeof(char*));
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
  for ( ; ; )
    {
      final = minpr <= FINAL_BOUND && minpa <= FINAL_BOUND;
      tot_alloc = tot_alloc0;
      
      fprintf (stderr, "Pass %d, filtering ideals >= %ld on rat. side and "
	       "%ld on alg. side:\n", ++pass, minpr, minpa);
      /*
        ret = scan_relations (fic, &nprimes, &H, rel_used, nrelmax,
	rel_compact, minpr, minpa, &tot_alloc, final);
      */
      prempt_scan_relations ();
      
      fprintf (stderr, "   nrels=%d, nprimes=%d; excess=%d\n",
	       nrel, nprimes, nrel - nprimes);
      
      fprintf (stderr, "   Starting singleton removal...\n");
      nrel_new = nrel;
      nprimes_new = nprimes;
      /* if one pass only, usually the initial excess is negative, but after
	 a few steps of removing singletons, we get a positive excess */
      if (final && (nrel_new < nprimes_new * excess) && pass >= 2)
	{
	  fprintf (stderr, "Initial excess is below requested %ld, stopping.\n",
		   (long int)(nprimes_new * (excess - 1)));
	  exit (1); /* initial excess is too small */
	}
      if (final && noclique)
	{
	  fprintf (stderr, "Initial excess is attained, purge will be completed at the next filtering.\n");
	  return 0;
	}
      remove_singletons (&nrel_new, nrelmax, &nprimes_new, &H, rel_used,
			 rel_compact, keep, final);
      
      fprintf (stderr, "   nrel=%d, nprimes=%d; excess=%d\n",
	       nrel_new, nprimes_new, nrel_new - nprimes_new);
      
      if (nrel_new <= nprimes_new) /* covers case nrel = nprimes = 0 */
	exit (1); /* if the excess is <= here, it will be more negative
		     when we decrease minpr and minpa in the next pass */
      
      my_malloc_free_all ();
      if (final)
	break;
      else
	minpr = minpa = FINAL_BOUND;
      hashClear (&H); /* Reset all tables to 0. FIXME: in the 2nd pass,
			 we might reuse some of the information computed
			 in the 1st pass (roots of primes, etc). */
      nrel = nrel_new;
    }
  
  fprintf (stderr, "Freeing rel_compact array...\n");
  /* we do not use it anymore */
  free (rel_compact);
  
  /* we renumber the primes in order of apparition in the hashtable */
  fprintf (stderr, "Renumbering primes...\n");
  renumber (&nprimes_new, &H, sos);
  
  /* reread the relation files and convert them to the new coding */
  fprintf (stderr, "Storing remaining relations...\n");
  reread (purgedname, fic, &H, rel_used, nrel_new, nprimes_new, raw);
  
  hashFree (&H);
  bit_vector_clear(rel_used);
  cado_poly_clear (pol);
  
  if (filelist) filelist_clear(fic);
  
  param_list_clear(pl);
  
  print_timing_and_memory (wct0);
  
  return 0;
}
