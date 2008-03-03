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
#include <assert.h>
#include <math.h>
#include <string.h>

#define WANT_ASSERT

#include "utils/utils.h"

#include "hashpair.h"
#include "files.h"

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

void
fprint_rat(int *table_ind, int *nb_coeff, relation_t rel, hashtable_t *H)
{
    int index, old_index, i, parity, nbc = *nb_coeff;

    index = getHashAddr(H, rel.rp[0].p, -2);
    old_index = index;
    parity = 1;
    
    for (i = 1; i < rel.nb_rp; ++i) {
	index = getHashAddr(H, rel.rp[i].p, -2);
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

void
fprint_alg(int *table_ind, int *nb_coeff, relation_t rel,
           /*tab_prime_t bad_primes,*/ hashtable_t *H)
{
    int i, index, old_index, parity, nbc = *nb_coeff;

#if 0
    if(rel.nb_ap == 0)
	fprintf(stderr, "WARNING: nb_ap == 0 in fprint_alg\n");
    i = 0;
    while (isBadPrime(rel.ap[i].p, bad_primes))
	i++;
    if(i >= rel.nb_ap)
	fprintf(stderr, "WARNING: i=%d >= nb_ap=%d in fprint_alg\n", i, rel.nb_ap);
#else
    for(i = 0; i < rel.nb_ap; i++)
      if(!isBadPrime(/*rel.ap[i].p, bad_primes*/))
	    break;
    if(i >= rel.nb_ap){
	fprintf(stderr, "WARNING: i=%d >= nb_ap=%d in fprint_alg\n", i, rel.nb_ap);
	return; // rare, but...
    }
#endif
    index = getHashAddr(H, rel.ap[i].p, rel.ap[i].r);
    old_index = index;
    parity = 1;
    i++;
    for (;i < rel.nb_ap; ++i) {
      if (isBadPrime(/*rel.ap[i].p, bad_primes*/))
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

void
fprint_free(int *table_ind, int *nb_coeff, relation_t rel, hashtable_t *H)
{
    long p = rel.a;
    int i, nbc = *nb_coeff, index;

    index = getHashAddr(H, p, -2);
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
void 
fprint_rel_row (FILE *file, int irel, relation_t rel,
                /*tab_prime_t bad_primes,*/ hashtable_t *H)
{
  int i;
  int *table_ind;
  int nb_coeff;

  //  fprintf(file, "%ld %lu ", rel.a, rel.b); // old stuff
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

int
specialHashInsert(hashtable_t *H, long p, unsigned long r, int irel)
{
    int h = getHashAddr(H, p, r);

    if(H->hashcount[h] == 0){
        // new empty place
        H->hashtab_p[h] = p;
        H->hashtab_r[h] = r;
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
void
insertNormalRelation(char *rel_used, int **rel_compact, int irel, int *nprimes,
                     /* tab_prime_t *bad_primes, */ hashtable_t *H,
                     relation_t *rel, unsigned long rlim, unsigned long alim,
                     long maxpr, long maxpa, int final)
{
    int *tmp = NULL, ltmp = 0, itmp, j, h, ok = 1;

#if LP_ONLY == 0
    rlim = alim = 0; // ohhhhhhhhhhhhh!
#endif
    reduce_exponents_mod2(rel);
    computeroots(rel);
    if(final){
	// first count number of "large" primes
	for(j = 0; j < rel->nb_rp; j++)
	    if(rel->rp[j].p > rlim){
              if (rel->rp[j].p > (unsigned long) maxpr)
		    ok = 0;
		ltmp++;
	    }
	for(j = 0; j < rel->nb_ap; j++)
	    if(rel->ap[j].p > alim){
              if(rel->ap[j].p > (unsigned long) maxpa)
		    ok = 0;
		ltmp++;
	    }
	if(ok == 0){
	    rel_used[irel] = 0;
	    rel_compact[irel] = NULL;
	    return;
	}
	if(ltmp == 0)
	    // ff??
	    rel_used[irel] = -1; // and tmp won't be used anyhow
	else
	    tmp = (int *)malloc((ltmp + 1) * sizeof(int));
	itmp = 0;
	for(j = 0; j < rel->nb_rp; j++){
	    h = hashInsert(H, rel->rp[j].p, -2);
	    if(H->hashcount[h] == 1)
		// new prime
		*nprimes += 1;
	    if(rel->rp[j].p <= rlim)
		continue;
	    tmp[itmp++] = h;
	}
	for(j = 0; j < rel->nb_ap; j++){
	    h = hashInsert(H, rel->ap[j].p, rel->ap[j].r);
	    if(H->hashcount[h] == 1){
		// new prime
		*nprimes += 1;
#if 0 /* commented out, waiting for the FIXME to be solved */
                /* FIXME: the following comparison can never be true, since
                   rel->ap[j].r is an unsigned long!!! */
		if(rel->ap[j].r == -1){
		    // bad prime
		    if (bad_primes->allocated == bad_primes->length) {
			bad_primes->allocated += 100;
			bad_primes->tab = (unsigned long *)realloc((void *)bad_primes->tab, (bad_primes->allocated)*sizeof(unsigned long));
			ASSERT (bad_primes->tab != NULL);
		    }
		    bad_primes->tab[bad_primes->length] = rel->ap[j].p;
		    bad_primes->length++;
		}
#endif
	    }
	    if(rel->ap[j].p <= alim)
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
	    h = specialHashInsert(H, rel->rp[j].p, -2, irel);
	    if(H->hashcount[h] < 0)
		// new prime
		*nprimes += 1;
	}
	for(j = 0; j < rel->nb_ap; j++){
	    h = specialHashInsert(H, rel->ap[j].p, rel->ap[j].r, irel);
	    if(H->hashcount[h] < 0){
		// new prime
		*nprimes += 1;
#if 0 /* commented out: waiting for FIXME to be solved */
                /* FIXME: the following comparison can never be true, since
                   rel->ap[j].r is an unsigned long!!! */
		if(rel->ap[j].r == -1){
		    // bad prime
		    if (bad_primes->allocated == bad_primes->length) {
			bad_primes->allocated += 100;
			bad_primes->tab = (unsigned long *)realloc((void *)bad_primes->tab, (bad_primes->allocated)*sizeof(unsigned long));
			ASSERT (bad_primes->tab != NULL);
		    }
		    bad_primes->tab[bad_primes->length] = rel->ap[j].p;
		    bad_primes->length++;
		}
#endif
	    }
	}
    }
}

// The information is stored in the ap[].p part, which is odd, but convenient.
// rel->ap.p[0..deg[ contains the deg roots of f(x) mod p, leading to deg
// ideals for the factorization of (p).
void
insertFreeRelation(char *rel_used, int **rel_compact, int irel, int *nprimes, hashtable_t *H, relation_t *rel, int final)
{
    long p = rel->a; // rel->b == 0
    int j, h, *tmp = NULL, itmp;

    if(final){
#if LP_ONLY
	rel_used[irel] = -1;
	rel_compact[irel] = NULL; // to be sure
#else
	tmp = (int *)malloc((1 + rel->nb_ap + 1) * sizeof(int));
	itmp = 0;
#endif
	h = hashInsert(H, p, -2);
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
	h = specialHashInsert(H, p, -2, irel);
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
int
scan_relations_from_file (int *irel, int *nrel, char *rel_used,
                          int **rel_compact, int *nprimes,
                          /* tab_prime_t *bad_primes, */ hashtable_t *H,
                          FILE *file, unsigned long rlim,
                          unsigned long alim, long maxpr, long maxpa,
                          int final)
{
    relation_t rel;
    int ret;

    while(1){
	ret = fread_relation (file, &rel);
	if(ret != 1)
	    break;
	*irel += 1;
	if(!(*irel % 100000))
	    fprintf(stderr, "nrel = %d at %2.2lf\n", *irel, seconds());
	rel_used[*irel] = 1;
	if(rel.b > 0)
          insertNormalRelation (rel_used, rel_compact, *irel, nprimes,
                                /*bad_primes,*/ H, &rel, rlim, alim,
                                maxpr, maxpa, final);
	else
	    insertFreeRelation(rel_used,rel_compact,*irel,nprimes,H,&rel,final);
	if(rel_used[*irel] <= 0)
	    // relation removed
	    *nrel -= 1;
	clear_relation(&rel);
    }
    return ret;
}

// Read all relations from file.
int
scan_relations (char *ficname[], int nbfic, int *nrel, int *nprimes,
                /*tab_prime_t *bad_primes,*/ hashtable_t *H,
                char *rel_used, int **rel_compact, unsigned long rlim,
                unsigned long alim, long maxpr, long maxpa, int final)
{
    FILE *relfile;
    int ret, irel = -1;
    int i;
    
    *nprimes = 0;

    ASSERT(nbfic > 0);
    ASSERT(!final || (final && rel_compact != NULL));
    for(i = 0; i < nbfic; i++){
	relfile = fopen(ficname[i], "r");
	if(relfile == NULL){
	    fprintf(stderr, "Pb opening file %s\n", ficname[i]);
	    exit(1);
	}
	fprintf(stderr, "Adding file %s\n", ficname[i]);
	ret = scan_relations_from_file (&irel, nrel, rel_used, rel_compact,
                                        nprimes, /* bad_primes,*/ H,
                                        relfile, rlim, alim, maxpr, maxpa,
                                        final);
	if (ret == 0) {
	    fprintf(stderr, "Warning: error when reading file %s\n", ficname[i]);
	    break;
	}
	fclose(relfile);
    }
    fprintf(stderr, "Scanned %d relations\n", irel+1);

    return (ret == -1);
}

int
has_singleton(int *tab, hashtable_t *H)
{
    int j;

    for(j = 0; tab[j] != -1; j++)
	if(H->hashcount[tab[j]] == 1)
	    return tab[j];
    return -1;
}

void delete_relation(int i, /*int hs,*/ int *nprimes, hashtable_t *H,
                     char *rel_used, int **rel_compact)
{
    int j, *tab = rel_compact[i];

    for(j = 0; tab[j] != -1; j++){
	H->hashcount[tab[j]]--;
	if(H->hashcount[tab[j]] == 0)
	    *nprimes -= 1;
    }
    free(rel_compact[i]);
    rel_compact[i] = NULL;
    rel_used[i] = 0;
}

int
removeSingletonIfAny(int *nprimes, hashtable_t *H, char *rel_used, int **rel_compact, int i)
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
	free(rel_compact[i]);
	rel_compact[i] = NULL;
	rel_used[i] = 0;
    }
    return (j0 > -1);
}

int
compare(const void *v1, const void *v2)
{
    int w1 = *((int*) v1);
    int w2 = *((int*) v2);

    return (w1 >= w2 ? -1 : 1);
}

// We sort rows w.r.t. their weight and decide which one to delete.
void
deleteHeavierRows(hashtable_t *H, int *nrel, int *nprimes, char *rel_used, int **rel_compact, int nrelmax, int keep)
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
	delete_relation(tmp[i+1], /*0,*/ nprimes, H, rel_used, rel_compact);
	*nrel -= 1;
    }
    free(tmp);
}

void
onepass_singleton_removal(int nrelmax, int *nrel, int *nprimes, hashtable_t *H, char *rel_used, int **rel_compact /*, int keep*/)
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
	    int h = has_singleton(rel_compact[i], H);
	    if(h >= 0){
#  if DEBUG >= 2
		printf("h = %d is single -> (%ld, %ld)\n", 
		       h, H->hashtab_p[h], H->hashtab_r[h]);
#  endif
		delete_relation(i, /*h,*/ nprimes, H, rel_used, rel_compact);
		*nrel -= 1;
	    }
#else // not tested yet!
	    if(removeSingletonIfAny(nprimes, H, rel_used, rel_compact, i))
		*nrel -= 1;
#endif
	}
    }
}

void
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

void
remove_singletons(int *nrel, int nrelmax, int *nprimes, hashtable_t *H, char *rel_used, int **rel_compact, int keep)
{
    int old, newnrel = *nrel, newnprimes = *nprimes, i;

    do{
	old = newnrel;
	deleteHeavierRows(H,&newnrel,&newnprimes,rel_used,rel_compact,nrelmax,keep);
	fprintf(stderr,"dHR: %d %d at %2.2lf\n",newnrel,newnprimes,seconds());
	onepass_singleton_removal(nrelmax, &newnrel, &newnprimes, H, rel_used, rel_compact /*, keep*/);
	fprintf(stderr, "new/old = %d/%d; newnprimes=%d at %2.2lf\n",
		newnrel, old, newnprimes, seconds());
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
void
renumber(int *nprimes, /*tab_prime_t bad_primes,*/ hashtable_t *H)
{
  unsigned int i;
  int nb;

    for(i = 0, nb = 1; i < H->hashmod; i++)
      if(isBadPrime(/*H->hashtab_p[i], bad_primes*/) || (H->hashcount[i] == 0))
	    H->hashcount[i] = -1;
	else{
#if 1
	    if(H->hashcount[i] <= 1)
		fprintf(stderr, "WARNING: H->hashcount[%d] = %d\n",
			i, H->hashcount[i]);
#else
	    ASSERT(H->hashcount[i] > 1);
#endif
	    if(H->hashcount[i] > 0)
		H->hashcount[i] = nb++;
	}
    nb--;
    fprintf(stderr, "nb = %d\n", nb);
    *nprimes = nb;
}

void
reread(char *ficname[], unsigned int nbfic, /*tab_prime_t bad_primes,*/
       hashtable_t *H, char *rel_used, int nrows /*, int cols*/)
{
    FILE *file;
    relation_t rel;
    int irel = -1, ret, nr = 0;
    unsigned int i;

    for(i = 0; i < nbfic; i++){
	file = fopen(ficname[i], "r");
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
		    fprint_rel_row(stdout, irel, rel, /*bad_primes,*/ H);
		    nr++;
		    if(nr >= nrows){
			ret = 0;
			break;
		    }
		}
		clear_relation(&rel);
	    }
	} while(ret == 1);
	fclose(file);
    }
}

void
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
	relfile = fopen(ficname[i], "r");
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
	fclose(relfile);
    }
    fprintf(stderr, "old_nrels=%d old_nprimes=%d\n", nrels, nprimes);
    fprintf(stderr, "new_nrels=%d new_nprimes=%d\n", newnrels, newnprimes);
    fprintf(stderr, "-nrels %d -nprimes %d\n", newnrels, newnprimes);
}

int
main(int argc, char **argv)
{
    tab_prime_t bad_primes;
    hashtable_t H;
    char **fic, *polyname = NULL;
    unsigned int nfic;
    char *rel_used;
    int **rel_compact = NULL;
    int ret, k;
    int nrel, nprimes = 0, final = 1;
    unsigned int nrelmax = 0, i;
    int nrel_new, nprimes_new, Hsize, Hsizer, Hsizea;
    long maxpr = 0, maxpa = 0, keep = -1; // maximum value for nrows-ncols
    cado_poly pol;
    
    if (argc == 1) {
	fprintf(stderr, "usage: %s [filename]\n", argv[0]);
	fprintf(stderr, "  stdin input is not yet available, sorry.\n");
	exit(1);
    } 
    if (argc < 4) {
	fprintf(stderr, "usage: %s poly nrel file1 ... filen\n", argv[0]);
	fprintf(stderr, "  if no filename is given, takes input on stdin\n");
	exit(1);
    }

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
    }
    if(keep == -1)
	keep = nrelmax;

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
    read_polynomial(pol, polyname);
    if(maxpr == 0)
	maxpr = 1L << pol[0].lpbr;
    if(maxpa == 0)
	maxpa = 1L << pol[0].lpba;
    fprintf(stderr, "Number of relations is %u\n", nrelmax);
    if(nprimes > 0)
	Hsize = nprimes;
    else{
	// estimating the number of primes
#if 1
	Hsizer = maxpr / ((int)log((double)maxpr));
	Hsizea = maxpa / ((int)log((double)maxpa));
	Hsize = Hsizer + Hsizea;
#else
	// TODO: use maxpr and maxpa
	Hsizer = (1<<pol[0].lpbr)/((int)(pol[0].lpbr * log(2.0)));
	Hsizea = (1<<pol[0].lpba)/((int)(pol[0].lpba * log(2.0)));
	Hsize = (Hsizer > Hsizea ? Hsizer : Hsizea);
#endif
    }
    fprintf(stderr, "initializing hash tables with Hsize=%d...\n", Hsize);
    hashInit(&H, Hsize);

    rel_used = (char *)malloc(nrelmax * sizeof(char));
    if (final)
	rel_compact = (int **)malloc(nrelmax * sizeof(int *));
    
    bad_primes.allocated = 100;
    bad_primes.length = 0;
    bad_primes.tab = (unsigned long *)malloc(bad_primes.allocated*sizeof(unsigned long));

    fprintf(stderr, "reading file of relations...\n");
    nrel = nrelmax;
    ret = scan_relations (fic, nfic, &nrel, &nprimes, /*&bad_primes,*/ &H,
                          rel_used, rel_compact, pol[0].rlim,
                          pol[0].alim, maxpr, maxpa, final);
    ASSERT (ret);
    
    fprintf(stderr, "nrel(useful)=%d, nprimes=%d (%d)\n",nrel,nprimes,Hsize);

#if DEBUG >= 2
    for(i = 0; i < hashmod; i++)
	if(hashcount[i] == 1)
	    printf("H[%ld, %ld] = 1\n", hashtab_p[i], hashtab_r[i]);
#endif

    if(final == 0)
	reduce(fic, nfic, &H, rel_used, nrelmax, nprimes);
    else{
	fprintf(stderr, "starting singleton removal...\n");
	nrel_new = nrel;
	nprimes_new = nprimes;
	remove_singletons(&nrel_new,nrelmax,&nprimes_new,&H,rel_used,rel_compact,keep);
	
	fprintf(stderr, "\nbadprimes = \n");
	for (i = 0; i < bad_primes.length; ++i){
	    fprintf(stderr, "%lu ", bad_primes.tab[i]);
	    nprimes_new--;
	}
	fprintf(stderr, "\n");
	
	fprintf(stderr, "nrel=%d, nprimes=%d; excess=%d\n", 
		nrel_new, nprimes_new, nrel_new - nprimes_new);

	if(nrel_new < nprimes_new)
	    exit(1);
	
	// we have to renumber the primes in increasing weight order
	renumber(&nprimes_new, /*bad_primes,*/ &H);
	
#if DEBUG >= 2
	hashCheck(nprimes_new);
#endif
	
	// we do not use it anymore
	for(i = 0; i < nrelmax; i++)
	    if(rel_compact[i] != NULL)
		free(rel_compact[i]);
	free(rel_compact);
    
	// now, we reread the file of relations and convert it to the new coding...
	printf("%d %d\n", nrel_new, nprimes_new);
	reread(fic, nfic, /*bad_primes,*/ &H, rel_used, nrel_new
               /*, nprimes_new*/);
    }

    free(bad_primes.tab);
    free(rel_used);

    return 0;
}
