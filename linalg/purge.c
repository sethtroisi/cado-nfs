/*
 * Program: filter
 * Author : F. Morain
 * Purpose: filtering
 *
 * Algorithm: Cavallar++
 *
 */

#include <gmp.h>
#include "mod_ul.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>

#include "utils/utils.h"

#include "hashpair.h"
#include "files.h"
#include "gzip.h"

#define DEBUG 1
#define LP_ONLY 0

// Data structure for one algebraic prime and for a table of those.
typedef struct {
  unsigned long prime;
  unsigned long root;
} rootprime_t;

typedef struct {
  rootprime_t *tab;
  unsigned long allocated;
  unsigned long length;
} tab_rootprime_t;

// Data structure for a table of rational primes.
typedef struct {
  unsigned long *tab;
  unsigned long allocated;
  unsigned long length;
} tab_prime_t;

// Data structure for an (a,b) pair and for a table of those
typedef struct {
  long a;
  unsigned long b;
} abpair_t;

typedef struct {
  abpair_t *tab;
  unsigned long allocated;
  unsigned long length;
} tab_abpair_t;

// Table of relations.
typedef struct {
  relation_t *tab;
  unsigned long allocated;
  unsigned long length;
} tab_relation_t;

/********************** own memory allocation routines ***********************/

#define MY_ALLOC

#ifdef MY_ALLOC
/* Rationale: calling one malloc() for each relation in insertNormalRelation()
   is expensive, since malloc() allocates some extra information to keep track
   of every memory blocks. Instead, we allocate memory in big blocks of size
   BLOCK_SIZE. */

#define BLOCK_SIZE 1000000

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
      /* first shrink current block */
      if (relcompact_current >= 0)
        relcompact_list[relcompact_current] = (int*) realloc (relcompact_list[relcompact_current], relcompact_used * sizeof (int));
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

static void
my_malloc_free_all (void)
{
  while (relcompact_current >= 0)
    free (relcompact_list[relcompact_current--]);
  free (relcompact_list);
}
#endif

/*****************************************************************************/

static int
isBadPrime(/*unsigned long p, tab_prime_t bad_primes*/) {
#if 1
    return 0;
#else
  int i;
  for (i = 0; i < bad_primes.length; ++i) {
    if (p == bad_primes.tab[i])
      return 1;
  }
  return 0;
#endif
}

static void
fprint_rat(int *table_ind, int *nb_coeff, relation_t rel, hashtable_t *H)
{
    int index, old_index, i, parity, nbc = *nb_coeff;
    uint64_t minus2 = (uint64_t) (-2);

    index = getHashAddr(H, rel.rp[0].p, minus2);
    old_index = index;
    parity = 1;

    for (i = 1; i < rel.nb_rp; ++i) {
	index = getHashAddr(H, rel.rp[i].p, minus2);
	if (index == old_index) {
	    parity = 1 - parity;
	} else {
	    if (parity == 1) {
		if(H->hashcount[old_index] >= 0)
		    table_ind[nbc++] = H->hashcount[old_index];
	    }
	    old_index = index;
	    parity = 1;
	}
    }
    if (parity == 1)
	if(H->hashcount[index] >= 0)
	    table_ind[nbc++] = H->hashcount[index];
    *nb_coeff = nbc;
}

static void
fprint_alg(int *table_ind, int *nb_coeff, relation_t rel,
           /*tab_prime_t bad_primes,*/ hashtable_t *H)
{
    int i, index, old_index, parity, nbc = *nb_coeff;

    for(i = 0; i < rel.nb_ap; i++)
      if(!isBadPrime())
	    break;
    if(i >= rel.nb_ap){
	fprintf(stderr, "WARNING: i=%d >= nb_ap=%d in fprint_alg\n", i, rel.nb_ap);
	return; // rare, but...
    }
    index = getHashAddr(H, rel.ap[i].p, rel.ap[i].r);
    old_index = index;
    parity = 1;
    i++;
    for (;i < rel.nb_ap; ++i) {
      if (isBadPrime())
	    continue;
	index = getHashAddr(H, rel.ap[i].p, rel.ap[i].r);
	if (index == old_index) {
	    parity = 1 - parity;
	} else {
	    if (parity == 1) {
		if(H->hashcount[old_index] >= 0)
		    table_ind[nbc++] = H->hashcount[old_index];
	    }
	    old_index = index;
	    parity = 1;
	}
    }
    if (parity == 1)
	if(H->hashcount[index] >= 0)
	    table_ind[nbc++] = H->hashcount[index];
    *nb_coeff = nbc;
}

static void
fprint_free(int *table_ind, int *nb_coeff, relation_t rel, hashtable_t *H)
{
    long p = rel.a;
    int i, nbc = *nb_coeff, index;
    uint64_t minus2 = (uint64_t) (-2);

    index = getHashAddr(H, p, minus2);
    if(H->hashcount[index] >= 0)
	table_ind[nbc++] = H->hashcount[index];
    for(i = 0; i < rel.nb_ap; i++){
	index = getHashAddr(H, p, rel.ap[i].p);
	if(H->hashcount[index] >= 0)
	    table_ind[nbc++] = H->hashcount[index];
    }
    *nb_coeff = nbc;
}

/* Print relations in a matrix format:
   don't take into account bad primes and even powers of primes.
   WARNING: the primes in the input relation are not necessarily sorted.
*/
static void
fprint_rel_row (FILE *file, int irel, relation_t rel,
                /*tab_prime_t bad_primes,*/ hashtable_t *H)
{
  int i;
  int *table_ind;
  int nb_coeff;

  table_ind = (int*) malloc((rel.nb_rp + rel.nb_ap)*sizeof(int));

  nb_coeff = 0;

  if(rel.b == 0)
      fprint_free(table_ind, &nb_coeff, rel, H);
  else{
      if((rel.nb_rp == 0) && (rel.b > 0))
	  fprintf(stderr, "WARNING: nb_rp = 0\n");
      else
	  fprint_rat(table_ind, &nb_coeff, rel, H);

      if(rel.nb_ap == 0)
	  fprintf(stderr, "WARNING: nb_ap=0\n");
      else
        fprint_alg(table_ind, &nb_coeff, rel, /*bad_primes,*/ H);
  }
  if(nb_coeff == 0)
      fprintf(stderr, "nb_coeff[%d] == 0\n", irel);
#if 0 // TODO: hum
  if(nb_coeff > 0){
      // nb_coeff = 0 is rather strange, no?
#endif
      fprintf(file, "%d %d", irel, nb_coeff);
      for (i = 0; i < nb_coeff; ++i) {
	  fprintf(file, " ");
	  // due to the +1 in renumber
	  fprintf(file, PURGE_INT_FORMAT, table_ind[i]-1);
      }
      fprintf(file, "\n");
#if 0 // TODO: hum
  }
#endif

  free(table_ind);
}

static int
specialHashInsert(hashtable_t *H, long p, unsigned long r, int irel)
{
    int h = getHashAddr(H, p, r);

    if(H->hashcount[h] == 0){
        // new empty place
        SET_HASH_P(H,h,p);
        SET_HASH_R(H,h,r);
	H->hashcount[h] = -(irel+1); // trick!!!
    }
    else
	// more than second time
	H->hashcount[h] = irel+1; // overcautious!!!
    return h;
}

// First we count the number of large primes; then we store all primes in
// the hash table, but not in the relations. This might end up with singletons
// here and there, but we don't care, since they will be dealt with in
// merge.
/* Meaning of the different parameters:
   minpr, minpa: only ideals > minpr (resp. minpa) are considered on the
                 rational (resp. algebraic) side. This means that the output
                 might contain ideals <= minpr or minpa appearing only once.
   maxpr, maxpa: relations with ideals > maxpr (resp. maxpa) on the rational
                 (resp. algebraic) side are discarded.
*/
static void
insertNormalRelation(char *rel_used, int **rel_compact, int irel, int *nprimes,
                     hashtable_t *H, relation_t *rel,
                     long maxpr, long maxpa, unsigned long minpr,
                     unsigned long minpa, int final, unsigned long *tot_alloc)
{
    int *tmp = NULL, ltmp = 0, itmp, j, h, ok = 1;
    /* special values of the roots are:
       -2: for rational primes, so that they don't conflict with algebraic
           primes,
       -1: for algebraic primes dividing lc(f).
       We apply a mask so that -2 in 64-bit translates to -2 in 32-bits. */
    uint64_t mask = (H->need64) ? (uint64_t) (-1) : (uint32_t) 4294967295;
    uint64_t minus2 = mask & (uint64_t) -2;

    reduce_exponents_mod2 (rel);
    computeroots (rel);
    if(final){
	// first count number of "large" primes
	for (j = 0; j < rel->nb_rp; j++)
          if (rel->rp[j].p > minpr) /* only consider primes > minpr */
            {
              if (rel->rp[j].p > (unsigned long) maxpr)
                ok = 0; /* discard relations with too large primes */
              ltmp++;
	    }
	for (j = 0; j < rel->nb_ap; j++)
          if (rel->ap[j].p > minpa) /* only consider primes > minpr */
            {
              if (rel->ap[j].p > (unsigned long) maxpa)
                ok = 0; /* discard relations with too large primes */
              ltmp++;
            }
	if(ok == 0){
	    rel_used[irel] = 0;
	    rel_compact[irel] = NULL;
	    return;
	}
        /* ltmp is the number of considered primes in the relation */
	if (ltmp == 0)
	    // ff??
	    rel_used[irel] = -1; // and tmp won't be used anyhow
	else
          {
#ifdef MY_ALLOC
            tmp = my_malloc_int (ltmp + 1);
#else
	    tmp = (int *) malloc ((ltmp + 1) * sizeof(int));
#endif
            *tot_alloc += (ltmp + 1) * sizeof(int);
          }
	itmp = 0;
	for(j = 0; j < rel->nb_rp; j++){
          /* trick: we use the same hash table for rational and algebraic
             primes, but use a fake root -2 for rational primes, which
             ensures there is no collision with algebraic primes */
            h = hashInsert(H, rel->rp[j].p, minus2);
	    if(H->hashcount[h] == 1)
		// new prime
		*nprimes += 1;
	    if (rel->rp[j].p <= minpr) /* don't consider small primes */
		continue;
	    tmp[itmp++] = h;
	}
	for(j = 0; j < rel->nb_ap; j++){
            /* apply mask to map 64-bit -1 to 32-bit -1 if needed */
	    h = hashInsert(H, rel->ap[j].p, mask & rel->ap[j].r);
	    if(H->hashcount[h] == 1)
              // new prime
              *nprimes += 1;
	    if(rel->ap[j].p <= minpa) /* don't consider small primes */
		continue;
	    tmp[itmp++] = h;
	}
	if(tmp != NULL){
	    tmp[itmp] = -1; // sentinel
	    rel_compact[irel] = tmp;
	}
	else
	    rel_compact[irel] = NULL; // to be sure
    }
    else{
	// reduce mode only
	for(j = 0; j < rel->nb_rp; j++){
	    h = specialHashInsert(H, rel->rp[j].p, minus2, irel);
	    if(H->hashcount[h] < 0)
		// new prime
		*nprimes += 1;
	}
	for(j = 0; j < rel->nb_ap; j++){
	    h = specialHashInsert(H, rel->ap[j].p, mask & rel->ap[j].r, irel);
	    if(H->hashcount[h] < 0){
		// new prime
		*nprimes += 1;
	    }
	}
    }
}

// The information is stored in the ap[].p part, which is odd, but convenient.
// rel->ap.p[0..deg[ contains the deg roots of f(x) mod p, leading to deg
// ideals for the factorization of (p).
static void
insertFreeRelation (char *rel_used, int **rel_compact, int irel, int *nprimes,
                    hashtable_t *H, relation_t *rel, int final,
                    unsigned long *tot_alloc)
{
    long p = rel->a; // rel->b == 0
    int j, h, *tmp = NULL, itmp;
    uint64_t minus2 = (H->need64) ? (uint64_t) (-2) : (uint32_t) (-2);

    if(final){
#if LP_ONLY
	rel_used[irel] = -1;
	rel_compact[irel] = NULL; // to be sure
#else
	tmp = (int *)malloc((1 + rel->nb_ap + 1) * sizeof(int));
        *tot_alloc += (1 + rel->nb_ap + 1) * sizeof(int);
	itmp = 0;
#endif
	h = hashInsert(H, p, minus2);
	if(H->hashcount[h] == 1)
	    // new prime
	    *nprimes += 1;
#if LP_ONLY == 0
	tmp[itmp++] = h;
#endif
	for(j = 0; j < rel->nb_ap; j++){
	    h = hashInsert(H, p, rel->ap[j].p);
	    if(H->hashcount[h] == 1)
		// new prime
		*nprimes += 1;
#if LP_ONLY == 0
	    tmp[itmp++] = h;
#endif
	}
#if LP_ONLY == 0
	tmp[itmp++] = -1;
	rel_used[irel] = 1;
	rel_compact[irel] = tmp;
#endif
    }
    else{
	// reduce mode
	h = specialHashInsert(H, p, minus2, irel);
	if(H->hashcount[h] < 0)
	    // new prime
	    *nprimes += 1;
	for(j = 0; j < rel->nb_ap; j++){
	    h = specialHashInsert(H, p, rel->ap[j].p, irel);
	    if(H->hashcount[h] < 0)
		// new prime
		*nprimes += 1;
	}
    }
}

// nrel can decrease, for instance in case of duplicate rows.
static int
scan_relations_from_file (int *irel, int *nrel, char *rel_used,
                          int **rel_compact, int *nprimes, hashtable_t *H,
                          FILE *file,
			  long maxpr, long maxpa, long minpr, long minpa,
                          int final, unsigned long *tot_alloc)
{
    relation_t rel;
    int ret;

    while(1){
	ret = fread_relation (file, &rel);
	if(ret != 1){
	    break;
	}
	*irel += 1;
	if(!(*irel % 100000))
	    fprintf(stderr, "nrel = %d at %2.2lfs (memory %luMb)\n",
                    *irel, seconds (), *tot_alloc / 1000000);
	rel_used[*irel] = 1;
	if(rel.b > 0)
          insertNormalRelation (rel_used, rel_compact, *irel, nprimes, H, &rel,
                                maxpr, maxpa, minpr, minpa, final, tot_alloc);
	else
          insertFreeRelation (rel_used, rel_compact, *irel, nprimes, H, &rel,
                              final, tot_alloc);
	if(rel_used[*irel] <= 0)
	    // relation removed
	    *nrel -= 1;
	clear_relation(&rel);
    }
    return ret;
}

/* Read all relations from file, and fills the rel_used and rel_compact arrays
   for each relation i:
   - rel_used[i] = 0 or -1 if the relation has no considered prime
       (i.e., a prime in the interval [minpr, maxpr] on the rational side),
       otherwise rel_used[i] = 1
   - rel_compact is an array, terminated by -1, of pointers to the entries
     in the hash table for the considered primes
 */
static int
scan_relations (char *ficname[], int nbfic, int *nrel, int *nprimes,
                hashtable_t *H,
                char *rel_used, int **rel_compact, long maxpr, long maxpa,
		long minpr, long minpa, int final, unsigned long *tot_alloc)
{
    FILE *relfile;
    int ret = 0, irel = -1;
    int i;

    *nprimes = 0;

    ASSERT(nbfic > 0);
    ASSERT(!final || (final && rel_compact != NULL));
    for(i = 0; i < nbfic; i++){
	relfile = gzip_open (ficname[i], "r");
	if(relfile == NULL){
	    fprintf(stderr, "Pb opening file %s\n", ficname[i]);
	    exit(1);
	}
	fprintf(stderr, "Adding file %s\n", ficname[i]);
	ret = scan_relations_from_file (&irel, nrel, rel_used, rel_compact,
                                        nprimes, H, relfile,
					maxpr, maxpa, minpr, minpa,
                                        final, tot_alloc);
	if (ret == 0) {
	    fprintf(stderr, "Warning: error when reading file %s\n", ficname[i]);
	    break;
	}
	gzip_close (relfile, ficname[i]);
    }
    fprintf(stderr, "Scanned %d relations\n", irel+1);

    return (ret == -1);
}

/* Return a positive value of some prime (ideal) in the array tab[] is single,
   otherwise return -1: tab[j] is the index in the hash table of the
   corresponding prime ideal.
*/
static int
has_singleton (int *tab, hashtable_t *H)
{
    int j;

    for(j = 0; tab[j] != -1; j++)
	if(H->hashcount[tab[j]] == 1)
	    return tab[j];
    return -1;
}

static void
delete_relation (int i, int *nprimes, hashtable_t *H,
                 char *rel_used, int **rel_compact)
{
    int j, *tab = rel_compact[i];

    for(j = 0; tab[j] != -1; j++){
	H->hashcount[tab[j]]--;
	if(H->hashcount[tab[j]] == 0)
	    *nprimes -= 1;
    }
#ifndef MY_ALLOC
    free (rel_compact[i]);
#endif
    rel_compact[i] = NULL;
    rel_used[i] = 0;
}

#if 0
static int
removeSingletonIfAny(int *nprimes, hashtable_t *H, char *rel_used,
                     int **rel_compact, int i)
{
    int j, j0 = -1, *tab = rel_compact[i];

    for(j = 0; tab[j] != -1; j++)
	if(j0 == -1){
	    if(H->hashcount[tab[j]] == 1){
		// first singleton found
		do{
		    j0++;
		    H->hashcount[tab[j0]]--;
		} while(j0 < j);
		*nprimes -= 1; // for that first j
	    }
	}
	else{
	    // removal being done...
	    H->hashcount[tab[j]]--;
	    if(H->hashcount[tab[j]] == 0)
		*nprimes -= 1;
	}
    if(j0 > -1){
#ifndef MY_ALLOC
	free(rel_compact[i]);
#endif
	rel_compact[i] = NULL;
	rel_used[i] = 0;
    }
    return (j0 > -1);
}
#endif

static int
compare(const void *v1, const void *v2)
{
    int w1 = *((int*) v1);
    int w2 = *((int*) v2);

    return (w1 >= w2 ? -1 : 1);
}

/* Delete the heavier relations if the current excess is larger than 'keep'
   We sort rows w.r.t. their weight and decide which one to delete. */
static void
deleteHeavierRows (hashtable_t *H, int *nrel, int *nprimes, char *rel_used,
                   int **rel_compact, int nrelmax, int keep)
{
    int *tmp, i, j = 0, nl;

    if((*nrel - *nprimes) < keep)
	return;

    tmp = (int *)malloc(((*nrel) << 1) * sizeof(int));
    for(i = 0; i < nrelmax; i++){
	if(rel_used[i] > 0){
	    for(nl = 0; rel_compact[i][nl] != -1; nl++);
	    tmp[j++] = nl;
	    tmp[j++] = i;
	}
    }
    qsort(tmp, *nrel, 2 * sizeof(int), compare);
    // first stupid idea: just remove heavy rows
    for(i = 0; i < *nrel; i += 2){
	if((*nrel)-(*nprimes) < keep)
	    break;
	delete_relation (tmp[i+1], nprimes, H, rel_used, rel_compact);
	*nrel -= 1;
    }
    free(tmp);
}

static void
onepass_singleton_removal (int nrelmax, int *nrel, int *nprimes,
                           hashtable_t *H, char *rel_used, int **rel_compact)
{
    int i;

    for(i = 0; i < nrelmax; i++){
#if DEBUG >= 2
	if(!(i % 1000000))
	    fprintf(stderr, "onepass: (i, nrel, nprimes)=(%d, %d, %d) at %2.2lf\n",
		    i, *nrel, *nprimes, seconds());
#endif
	if(rel_used[i] > 0){
#if 1
	    int h = has_singleton (rel_compact[i], H);
	    if(h >= 0){
#  if DEBUG >= 2
		printf("h = %d is single -> (%ld, %ld)\n",
		       h, GET_HASH_P(H,h), GET_HASH_R(H,h));
#  endif
		delete_relation (i, nprimes, H, rel_used, rel_compact);
		*nrel -= 1;
	    }
#else // not tested yet!
	    if(removeSingletonIfAny(nprimes, H, rel_used, rel_compact, i))
		*nrel -= 1;
#endif
	}
    }
}

#if 0
static void
find(int nrel, int h, char *rel_used, int **rel_compact)
{
    int i;

    for(i = 0; i < nrel; i++)
	if(rel_used[i] > 0){
	    int *tab = rel_compact[i], j;
	    for(j = 0; tab[j] != -1; j++)
		if(tab[j] == h)
		    fprintf(stderr, "GATO_Z: %d\n", i);
	}
}
#endif

static void
remove_singletons (int *nrel, int nrelmax, int *nprimes, hashtable_t *H,
                   char *rel_used, int **rel_compact, int keep)
{
    int old, newnrel = *nrel, newnprimes = *nprimes, i;

    do{
	old = newnrel;
	deleteHeavierRows (H, &newnrel, &newnprimes, rel_used, rel_compact,
                           nrelmax, keep);
	if(newnrel != old)
	    fprintf (stderr, "deleted heavier relations: %d %d at %2.2lf\n",
                     newnrel, newnprimes, seconds ());
	onepass_singleton_removal (nrelmax, &newnrel, &newnprimes, H,
                                   rel_used, rel_compact);
	fprintf(stderr, "new_nrows=%d new_ncols=%d (%d) at %2.2lf\n",
		newnrel, newnprimes, newnrel-newnprimes, seconds());
    } while(newnrel != old);
    // clean empty rows
    for(i = 0; i < nrelmax; i++){
	if(rel_used[i] > 0){
	    if(rel_compact[i][0] == -1){
		fprintf(stderr, "Empty row: %d\n", i);
	    }
	}
    }
    *nrel = newnrel;
    *nprimes = newnprimes;
}

// we locate used primes and do not try to do fancy things as sorting w.r.t.
// weight, since this will probably be done later on.
// All rows will be 1 more that needed -> subtract 1 in fprint...!
static void
renumber(int *nprimes, hashtable_t *H, char *sos)
{
    FILE *fsos = NULL;
    unsigned int i;
    int nb;

    if(sos != NULL){
	fprintf(stderr, "Outputting renumber table in file %s\n", sos);
	fsos = gzip_open (sos, "w");
    }
    for(i = 0, nb = 1; i < H->hashmod; i++)
      if(isBadPrime() || (H->hashcount[i] == 0))
	    H->hashcount[i] = -1;
	else{
#if 1
	    if(H->hashcount[i] <= 1)
		fprintf(stderr, "WARNING: H->hashcount[%d] = %d\n",
			i, H->hashcount[i]);
#else
	    ASSERT(H->hashcount[i] > 1);
#endif
	    if(H->hashcount[i] > 0){
		H->hashcount[i] = nb++;
		if(fsos != NULL)
		    fprintf(fsos, "%d %"PRIi64" %"PRIu64"\n",
			    H->hashcount[i]-1, GET_HASH_P(H,i),
                            GET_HASH_R(H,i));
	    }
	}
    if(fsos != NULL)
      gzip_close (fsos, sos);
    nb--;
    fprintf(stderr, "nb = %d\n", nb);
    *nprimes = nb;
}

static void
reread(FILE *ofile, char *ficname[], unsigned int nbfic,
       /*tab_prime_t bad_primes,*/
       hashtable_t *H, char *rel_used, int nrows /*, int cols*/)
{
    FILE *file;
    relation_t rel;
    int irel = -1, ret, nr = 0;
    unsigned int i;

    for(i = 0; i < nbfic; i++){
	file = gzip_open (ficname[i], "r");
	if (file == NULL) {
	    fprintf(stderr,"problem opening file %s for reading\n",ficname[i]);
	    exit(1);
	}
	fprintf(stderr, "Adding file %s\n", ficname[i]);
	do{
	    ret = fread_relation (file, &rel);
	    if(ret == 1){
		irel++;
		if(!(irel % 100000))
		    fprintf(stderr, "nrel = %d at %2.2lf\n", irel, seconds());
		if(rel_used[irel]){ // covers -1 and > 0
		    if(rel.b > 0){
			reduce_exponents_mod2 (&rel);
			computeroots (&rel);
		    }
		    fprint_rel_row(ofile, irel, rel, /*bad_primes,*/ H);
		    nr++;
		    if(nr >= nrows){
			ret = 0;
			break;
		    }
		}
		clear_relation(&rel);
	    }
	} while(ret == 1);
	gzip_close (file, ficname[i]);
    }
}

static void
reduce(char **ficname, unsigned int nbfic, hashtable_t *H, char *rel_used,
       int nrels, int nprimes)
{
    FILE *relfile;
    unsigned int i;
    int nr = 0, newnrels = 0, newnprimes = nprimes;
    char str[1024];

    // locate singletons
    for(i = 0; i < H->hashmod; i++)
	if(H->hashcount[i] < 0){
	    rel_used[-H->hashcount[i]-1] = 0; // that trick again!
	    newnprimes--;
	}
    // I/O
    for(i = 0; i < nbfic; i++){
	relfile = gzip_open (ficname[i], "r");
	if(relfile == NULL){
	    fprintf(stderr, "Pb opening file %s\n", ficname[i]);
	    exit(1);
	}
	fprintf(stderr, "Extracting and printing data from file %s\n",
		ficname[i]);
	while(fgets(str, 1024, relfile)){
	    if(rel_used[nr++]){
		newnrels++;
		printf("%s", str);
	    }
	}
	gzip_close (relfile, ficname[i]);
    }
    fprintf(stderr, "old_nrels=%d old_nprimes=%d\n", nrels, nprimes);
    fprintf(stderr, "new_nrels=%d new_nprimes=%d\n", newnrels, newnprimes);
    fprintf(stderr, "-nrels %d -nprimes %d\n", newnrels, newnprimes);
}

static void
usage (char *argv[])
{
  fprintf (stderr, "Usage: %s [options] -poly polyfile -purged purgedfile -nrels nnn file1 ... filen\n", argv[0]);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "       -nonfinal    - perform only one singleton pass\n");
  fprintf (stderr, "       -keep    nnn - stop when excess <= nnn (default -1)\n");
  fprintf (stderr, "       -maxpa   nnn - discard rels with alg. primes <= nnn (default 2^lpba)\n");
  fprintf (stderr, "       -maxpr   nnn - discard rels with rat. primes <= nnn (default 2^lpbr)\n");
  fprintf (stderr, "       -minpa   nnn - purge alg. primes >= nnn only (default 0)\n");
  fprintf (stderr, "       -minpr   nnn - purge rat. primes >= nnn only (default 0)\n");
  fprintf (stderr, "       -nprimes nnn - number of prime ideals\n");
  fprintf (stderr, "       -sos sosfile - to keep track of the renumbering\n");
}

/* Estimate the number of primes <= B. The 0.85 factor accounts for the
   fact that the combinatorial explosion happens before B/log(B). */
static int
approx_phi (long B)
{
  return (B <= 1) ? 0 : (int) (0.85 * (double) B / log ((double) B));
}

int
main(int argc, char **argv)
{
    FILE *purgedfile = NULL;
    tab_prime_t bad_primes;
    hashtable_t H;
    char **fic, *polyname = NULL, *sos = NULL, *purgedname = NULL;
    unsigned int nfic;
    char *rel_used;
    int **rel_compact = NULL;
    int ret, k;
    int nrel, nprimes = 0, final = 1;
    unsigned long int nrelmax = 0, i;
    int nrel_new, nprimes_new, Hsize, Hsizer, Hsizea;
    long maxpr = 0, maxpa = 0, keep = -1; // maximum value for nrows-ncols
    long minpr = 0, minpa = 0;
    cado_poly pol;
    unsigned long tot_alloc;
    int need64 = 0; /* non-zero if large primes are > 2^32 */

    fprintf (stderr, "%s.r%s", argv[0], REV);
    for (k = 1; k < argc; k++)
      fprintf (stderr, " %s", argv[k]);
    fprintf (stderr, "\n");

    while(argc > 1 && argv[1][0] == '-'){
	if(argc > 2 && strcmp (argv[1], "-poly") == 0){
	    polyname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 2 && strcmp (argv[1], "-nrels") == 0){
	    nrelmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 2 && strcmp (argv[1], "-nprimes") == 0){
	    nprimes = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 2 && strcmp (argv[1], "-maxpr") == 0){
	    maxpr = atol(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 2 && strcmp (argv[1], "-maxpa") == 0){
	    maxpa = atol(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 2 && strcmp (argv[1], "-minpr") == 0){
	    minpr = atol(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 2 && strcmp (argv[1], "-minpa") == 0){
	    minpa = atol(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 2 && strcmp (argv[1], "-keep") == 0){
	    keep = atol(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 1 && strcmp (argv[1], "-nonfinal") == 0){
	    final = 0;
	    argc--;
	    argv++;
	}
	else if(argc > 1 && strcmp (argv[1], "-sos") == 0){
	    sos = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 1 && strcmp (argv[1], "-purged") == 0){
	    purgedname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
    }
    if (keep <= 0)
	keep = nrelmax;

    if (argc == 1) {
        usage (argv);
	exit(1);
    }

    if (polyname == NULL)
      {
        fprintf (stderr, "Error, missing -poly ... option\n");
        exit (1);
      }

    if (nrelmax == 0)
      {
        fprintf (stderr, "Error, missing -nrels ... option\n");
        exit (1);
      }

    fic = argv+1;
    nfic = argc-1;
    //	fic = extractFic(&nfic, &nrelmax, argv[3]);
    cado_poly_init(pol);
    cado_poly_read(pol, polyname);

    if (maxpr == 0)
	maxpr = 1L << pol[0].lpbr;
    if (maxpa == 0)
	maxpa = 1L << pol[0].lpba;
    need64 = (maxpr >> 32) || (maxpa >> 32);

    fprintf(stderr, "Number of relations is %lu\n", nrelmax);
    if(nprimes > 0)
	Hsize = nprimes;
    else{
      /* Estimating the number of needed primes (remember that hashInit
         multiplies by a factor 1.5). */
      Hsizer = approx_phi (maxpr) - approx_phi (minpr);
      Hsizea = approx_phi (maxpa) - approx_phi (minpa);
      Hsize = Hsizer + Hsizea;
    }
    fprintf (stderr, "initializing hash tables with Hsize=%d...\n", Hsize);
    hashInit (&H, Hsize, 1, need64);
    tot_alloc = H.hashmod * H.size;

    rel_used = (char *) malloc (nrelmax * sizeof (char));
    tot_alloc += nrelmax * sizeof (char);
    fprintf (stderr, "Allocated rel_used of %luMb (total %luMb so far)\n",
             (nrelmax * sizeof (char)) / 1000000,
             tot_alloc / 1000000);
    if (final)
      {
        /* FIXME: the rel_compact alone uses a lot of memory. For 100M
           relations (which is reasonable for a 155-digit number with large
           prime bounds 2^30), on a 64-bit machine it uses 0.8Gb!!! */
	rel_compact = (int **) malloc (nrelmax * sizeof (int *));
        tot_alloc += nrelmax * sizeof (int*);
        fprintf (stderr, "Allocated rel_compact of %luMb (total %luMb so far)\n",
                 (nrelmax * sizeof (int *)) / 1000000,
                 tot_alloc / 1000000);
      }

    bad_primes.allocated = 100;
    bad_primes.length = 0;
    bad_primes.tab = (unsigned long *) malloc (bad_primes.allocated
                                               * sizeof(unsigned long));

    fprintf(stderr, "reading file of relations...\n");
    nrel = nrelmax;
    ret = scan_relations (fic, nfic, &nrel, &nprimes, &H,
                          rel_used, rel_compact,
                          maxpr, maxpa, minpr, minpa, final, &tot_alloc);
    ASSERT (ret);

    fprintf (stderr, "nrel(useful)=%d, nprimes=%d (expected %d)\n",
             nrel, nprimes, Hsize);

#if DEBUG >= 2
    for(i = 0; i < hashmod; i++)
	if(hashcount[i] == 1)
	    printf("H[%ld, %ld] = 1\n", hashtab_p[i], hashtab_r[i]);
#endif

    if(final == 0)
	reduce (fic, nfic, &H, rel_used, nrelmax, nprimes);
    else{
	fprintf(stderr, "starting singleton removal...\n");
	nrel_new = nrel;
	nprimes_new = nprimes;
	remove_singletons (&nrel_new, nrelmax, &nprimes_new, &H, rel_used,
                           rel_compact, keep);

	fprintf(stderr, "\nbadprimes = \n");
	for (i = 0; i < bad_primes.length; ++i){
	    fprintf(stderr, "%lu ", bad_primes.tab[i]);
	    nprimes_new--;
	}
	fprintf(stderr, "\n");

	fprintf(stderr, "nrel=%d, nprimes=%d; excess=%d\n",
		nrel_new, nprimes_new, nrel_new - nprimes_new);

	if (nrel_new <= nprimes_new) /* covers the case nrel = nprimes = 0 */
	    exit(1);

	// we renumber the primes in order of apparition in the hashtable
        fprintf (stderr, "Renumbering primes...\n");
	renumber (&nprimes_new, &H, sos);

#if DEBUG >= 2
	hashCheck(nprimes_new);
#endif

        fprintf (stderr, "Freeing rel_compact array...\n");
	// we do not use it anymore
#ifdef MY_ALLOC
        my_malloc_free_all ();
#else
	for(i = 0; i < nrelmax; i++)
	    if(rel_compact[i] != NULL)
		free(rel_compact[i]);
#endif
	free(rel_compact);

	// now, we reread the file of relations and convert it to the new coding...
        fprintf (stderr, "Storing remaining relations...\n");
	purgedfile = gzip_open (purgedname, "w");
	fprintf (purgedfile, "%d %d\n", nrel_new, nprimes_new);
	reread (purgedfile, fic, nfic, /*bad_primes,*/ &H, rel_used, nrel_new
               /*, nprimes_new*/);
	gzip_close (purgedfile, purgedname);
	// write excess to stdout
	printf("EXCESS: %d\n", nrel_new - nprimes_new);
    }
    free(bad_primes.tab);
    free(rel_used);

    return 0;
}
