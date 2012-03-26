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
#include <string.h>
#include <pthread.h>

#include "utils.h"
#include "hashpair.h"
#include "gzip.h"
#include "filter_matrix.h" /* for BURIED_MAX_DENSITY */

#define MAX_FILES 1000000
#define MAX_STEPS 10   /* maximal number of steps in each pass */
#define FINAL_BOUND 0  /* final bound for rational and algebraic ideals */
#define MAX_THREADS 8

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER; /* mutual exclusion lock */

/********************** own memory allocation routines ***********************/

/* Rationale: calling one malloc() for each relation in insertNormalRelation()
   is expensive, since malloc() allocates some extra information to keep track
   of every memory blocks. Instead, we allocate memory in big blocks of size
   BLOCK_SIZE. */

#define BLOCK_SIZE 1000000 /* memory blocks are allocated of that # of int's */

/* relcompact_list is a list of blocks, each one of BLOCK_SIZE ints */
int **relcompact_list = NULL;
int relcompact_current = -1; /* index of current block */
unsigned long relcompact_used = BLOCK_SIZE; /* usage of current block */

/* return a pointer to an array of n ints */
static int*
my_malloc_int (unsigned long n)
{
  int *ptr;

  if (relcompact_used + n > BLOCK_SIZE)
    {
      /* allocate a new block */
      relcompact_current ++;
      relcompact_list = (int**) realloc (relcompact_list, (relcompact_current + 1) * sizeof (int*));
      relcompact_list[relcompact_current] = (int*) malloc (BLOCK_SIZE * sizeof (int));
      relcompact_used = 0;
    }
  ptr = relcompact_list[relcompact_current] + relcompact_used;
  relcompact_used += n;
  return ptr;
}

/*****************************************************************************/

/* dirty trick to distinguish rational primes: we store -2 for their root */
static unsigned long minus2;

/* First we count the number of large primes; then we store all primes in
   the hash table, but not in the relations. This might end up with singletons
   here and there, but we don't care, since they will be dealt with in
   merge.

   Meaning of the different parameters:
   minpr, minpa: only ideals > minpr (resp. minpa) are considered on the
                 rational (resp. algebraic) side. This means that the output
                 might contain ideals <= minpr or minpa appearing only once.
*/
static void
insertNormalRelation (int **rel_compact, int irel,
                      int *nprimes, hashtable_t *H, relation_t *rel,
                      unsigned long minpr, unsigned long minpa,
                      unsigned long *tot_alloc, int final)
{
  int *tmp = NULL, ltmp = 0, i, j, h, ok;

    reduce_exponents_mod2 (rel);

    /* first count number of "large" primes */
    for (j = 0; j < rel->nb_rp; j++)
      ltmp += (rel->rp[j].p >= minpr); /* only consider primes >= minpr */
    for (j = 0; j < rel->nb_ap; j++)
      ltmp += (rel->ap[j].p >= minpa); /* only consider primes >= minpa */

    /* ltmp is the number of considered primes in the relation.
       We might have ltmp=0 if all primes are less than minpr, maxpr. */
    pthread_mutex_lock (&lock);
    tmp = my_malloc_int (ltmp + 1);
    pthread_mutex_unlock (&lock);
    *tot_alloc += (ltmp + 1) * sizeof (int);
    i = 0; /* number of entries in tmp */

    /* trick: we use the same hash table for rational and algebraic
       primes, but use a fake root -2 for rational primes, which
       ensures there is no collision with algebraic primes */
    for (j = 0; j < rel->nb_rp; j++)
      {
        ok = rel->rp[j].p >= minpr;
        if (ok || final)
          /* in the final pass, we need to insert all ideals, since we need
             to renumber them in the output file */
          h = hashInsert (H, rel->rp[j].p, minus2);
        if (ok)
          {
            *nprimes += (H->hashcount[h] == 1); /* new prime */
            tmp[i++] = h;
          }
      }

    for (j = 0; j < rel->nb_ap; j++)
      {
        ok = rel->ap[j].p >= minpa;
        if (ok || final)
          {
            rel->ap[j].r = findroot (rel->a, rel->b, rel->ap[j].p);
            h = hashInsert (H, rel->ap[j].p, rel->ap[j].r);
          }
        if (ok)
          {
            *nprimes += (H->hashcount[h] == 1); /* new prime */
            tmp[i++] = h;
          }
      }

    tmp[i] = -1; /* sentinel */
    rel_compact[irel] = tmp;
}

/* The information is stored in the ap[].p part, which is odd, but convenient.
   rel->ap.p[0..deg[ contains the deg roots of f(x) mod p, leading to deg
   ideals for the factorization of (p). */
static void
insertFreeRelation (int **rel_compact, int irel, int *nprimes,
                    hashtable_t *H, relation_t *rel,
                    unsigned long minpr, unsigned long minpa,
                    unsigned long *tot_alloc, int final)
{
    unsigned long p = rel->a; /* rel->b == 0 */
    int j, h, *tmp = NULL, itmp, ltmp;

    /* the prime on the rational side is rel->a
       the prime ideal on the algebraic side are (rel->a, rel->ap[j].p) */
    ltmp = 1 + (p >= minpr);
    if (p >= minpa)
      ltmp += rel->nb_ap;
    tmp = my_malloc_int (ltmp);
    *tot_alloc += ltmp * sizeof (int);
    itmp = 0;
    if ((p >= minpr) || final)
      h = hashInsert (H, p, minus2);
    if (p >= minpr)
      {
        *nprimes += (H->hashcount[h] == 1); /* new prime */
        tmp[itmp++] = h;
      }
    for (j = 0; j < rel->nb_ap; j++)
      {
        if ((p >= minpa) || final)
          h = hashInsert(H, p, rel->ap[j].p);
        if (p >= minpa)
          {
            *nprimes += (H->hashcount[h] == 1); /* new ideal */
            tmp[itmp++] = h;
          }
      }
    tmp[itmp++] = -1;
    rel_compact[irel] = tmp;
}

static unsigned long
scan_files_mod (char **ficname, bit_vector_srcptr rel_used, int **rel_compact,
                int *nprimes, hashtable_t *H, long minpr, long minpa,
                unsigned long *tot_alloc, int final, int thread)
{
  int i;
  relation_stream rs;
  unsigned long nrels = 0;
  char *c;

  /* first we go to file 'thread' */
  for (i = 0; (i < thread) && (ficname[0] != NULL); i++, ficname++);

  while (ficname[0] != NULL)
    {
      relation_stream_init (rs);
      relation_stream_openfile (rs, ficname[0]);
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
            if (rs->rel.b > 0)
              insertNormalRelation (rel_compact, irel, nprimes,
                                H, &(rs->rel), minpr, minpa, tot_alloc, final);
            else
              insertFreeRelation (rel_compact, irel, nprimes,
                                H, &(rs->rel), minpr, minpa, tot_alloc, final);
          }
      }
      for (c = ficname[0] + strlen (ficname[0]); (c > ficname[0]) &&
             (c[-1] != '/'); c--);
      fprintf (stderr, "%d: %-50s %d relations\n", thread, c, rs->nrels);
      relation_stream_closefile (rs);
      relation_stream_clear (rs);

      nrels += rs->nrels;

      /* go to next file */
      for (i = 0; (i < MAX_THREADS) && (ficname[0] != NULL); i++, ficname++);
    }

  return nrels;
}

/* thread structure */
typedef struct
{
  char **ficname;
  bit_vector_srcptr rel_used;
  int **rel_compact;
  int *nprimes;
  hashtable_t *H;
  long minpr;
  long minpa;
  unsigned long *tot_alloc;
  int final;
  int thread;
  unsigned long nrels;
} __tab_struct;
typedef __tab_struct tab_t[1];

void*
one_thread (void* args)
{
  tab_t *tab = (tab_t*) args;

  tab[0]->nrels = scan_files_mod (tab[0]->ficname, tab[0]->rel_used,
                                  tab[0]->rel_compact, tab[0]->nprimes,
                                  tab[0]->H, tab[0]->minpr, tab[0]->minpa,
                                  tab[0]->tot_alloc, tab[0]->final,
                                  tab[0]->thread);
  return NULL;
}

/* Read all relations from file, and fills the rel_used and rel_compact arrays
   for each relation i:
   - rel_used[i] = 0 if the relation i is deleted
     rel_used[i] = 1 if the relation i is kept (so far)
   - rel_compact is an array, terminated by -1, of pointers to the entries
     in the hash table for the considered primes

     Trick: we only read relations for which rel_used[i]=1.
 */
static int
scan_relations (char **ficname, int *nprimes, hashtable_t *H,
                bit_vector_srcptr rel_used, int nrelmax, int **rel_compact,
		long minpr, long minpa, unsigned long *tot_alloc, int final)
{
  int i = 0;
  tab_t *T;
  pthread_t tid[MAX_THREADS];
  unsigned long nrels = 0; /* total number of relations */

  *nprimes = 0;

  ASSERT(rel_compact != NULL);

  T = malloc (MAX_THREADS * sizeof (tab_t));
  for (i = 0; i < MAX_THREADS; i++)
    {
      T[i]->rel_used = rel_used;
      T[i]->rel_compact = rel_compact;
      T[i]->nprimes = nprimes;
      T[i]->H = H;
      T[i]->minpr = minpr;
      T[i]->minpa = minpa;
      T[i]->tot_alloc = tot_alloc;
      T[i]->final = final;
    }

  for (i = 0; i < MAX_THREADS; i++)
    {
      T[i]->ficname = ficname;
      T[i]->thread = i;
      pthread_create (&tid[i], NULL, one_thread, (void *) (T+i));
    }

  for (i = 0; i < MAX_THREADS; i++)
    pthread_join (tid[i], NULL);

  for (i = 0; i < MAX_THREADS; i++)
    {
      fprintf (stderr, "Thread %d read %lu relations\n", i, T[i]->nrels);
      nrels += T[i]->nrels;
    }

  free (T);
  
  fprintf (stderr, "Read a total of %lu relations\n", nrels);
  
  return 1;
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

typedef struct {
  float w;
  uint32_t i;
} comp_t;

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
    hashtable_t H;
    int **rel_compact = NULL;
    int ret MAYBE_UNUSED, k;
    int nprimes = 0;
    unsigned int nrelmax = 0;
    int Hsize, Hsizer, Hsizea;
    double excess = 1.01;    /* minimum initial excess */
    long keep = 160;    /* maximum final excess */
    long minpr = -1, minpa = -1; /* negative values mean use minpr=rlim and
                                    minpa=alim in first pass */
    cado_poly pol;
    unsigned long tot_alloc, tot_alloc0;
    int need64 = 0; /* non-zero if large primes are > 2^32 */
    int final = 0;
    int raw = 0;
    char ** fic;
    int noclique = 0;
    double wct0 = wct_seconds ();

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

    param_list_parse_uint(pl, "nrels", &nrelmax);
    param_list_parse_int(pl, "nprimes", &nprimes);
    param_list_parse_long(pl, "minpr", &minpr);
    param_list_parse_long(pl, "minpa", &minpa);
    param_list_parse_double(pl, "excess", &excess);
    param_list_parse_long(pl, "keep", &keep);
    param_list_parse_int(pl, "noclique", &noclique);

    const char * filelist = param_list_lookup_string(pl, "filelist");
    const char * basepath = param_list_lookup_string(pl, "basepath");
    const char * subdirlist = param_list_lookup_string(pl, "subdirlist");

    cado_poly_init (pol);

    const char * tmp;

    ASSERT_ALWAYS((tmp = param_list_lookup_string(pl, "poly")) != NULL);
    cado_poly_read(pol, tmp);

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

    fprintf (stderr, "Number of relations is %u\n", nrelmax);
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

    bit_vector rel_used;
    bit_vector_init_set(rel_used, nrelmax, 1);
    tot_alloc0 += nrelmax;
    fprintf (stderr, "Allocated rel_used of %uMb (total %luMb so far)\n",
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

    tot_alloc = tot_alloc0;

    ret = scan_relations (fic, &nprimes, &H, rel_used, nrelmax,
                          rel_compact, minpr, minpa, &tot_alloc, final);
    ASSERT (ret);

    fprintf (stderr, "Freeing rel_compact array...\n");
    /* we do not use it anymore */
    free (rel_compact);

    hashFree (&H);
    bit_vector_clear(rel_used);
    cado_poly_clear (pol);

    if (filelist) filelist_clear(fic);

    param_list_clear(pl);

    print_timing_and_memory (wct0);

    return 0;
}
