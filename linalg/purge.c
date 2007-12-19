/* 
 * Program: filter
 * Author : F. Morain
 * Purpose: filtering
 * 
 * Algorithm: Cavallar++
 *
 */


#include "cado.h"
#include <gmp.h>
#include "mod_ul.c"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "utils/utils.h"
#include "hashpair.h"

#include <string.h>

#define DEBUG 1

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

// Data strucutre for an (a,b) pair and for a table of those
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
isBadPrime(unsigned long p, tab_prime_t bad_primes) {
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
fprint_alg(int *table_ind, int *nb_coeff, relation_t rel, tab_prime_t bad_primes, hashtable_t *H)
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
	if(!isBadPrime(rel.ap[i].p, bad_primes))
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
	if (isBadPrime(rel.ap[i].p, bad_primes))
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
fprint_rel_row (FILE *file, int irel, relation_t rel, tab_prime_t bad_primes, hashtable_t *H)
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
	  fprint_alg(table_ind, &nb_coeff, rel, bad_primes, H);
  }
  if(nb_coeff == 0)
      fprintf(stderr, "nb_coeff[%d] == 0\n", irel);
#if 0 // TODO: hum
  if(nb_coeff > 0){
      // nb_coeff = 0 is rather strange, no?
#endif
      fprintf(file, "%d %d", irel, nb_coeff);
      for (i = 0; i < nb_coeff; ++i) {
	  fprintf(file, " %d", table_ind[i]-1); // due to the +1 in renumber
      }
      fprintf(file, "\n");
#if 0 // TODO: hum
  }
#endif

  free(table_ind);
}

/* reduces exponents mod 2, and discards primes with even exponent */
void
reduce_exponents_mod2 (relation_t *rel)
{
  int i, j;

  if(rel->nb_rp == 0)
      fprintf(stderr, "WARNING: nb_rp = 0 in reduce_exponents_mod2\n");

  for (i = j = 0; i < rel->nb_rp; i++)
    {
      rel->rp[i].e &= 1;
      if (rel->rp[i].e != 0)
        {
          rel->rp[j].p = rel->rp[i].p;
          rel->rp[j].e = 1;
          j ++;
        }
    }
  if(j == 0)
      fprintf(stderr, "WARNING: j_rp=0 in reduce_exponents_mod2\n");
  else
      rel->rp = (rat_prime_t*) realloc (rel->rp, j * sizeof (rat_prime_t));
  rel->nb_rp = j;

  if(rel->nb_ap == 0)
      fprintf(stderr, "WARNING: nb_ap = 0 in reduce_exponents_mod2\n");

  for (i = j = 0; i < rel->nb_ap; i++)
    {
      rel->ap[i].e &= 1;
      if (rel->ap[i].e != 0)
        {
          rel->ap[j].p = rel->ap[i].p;
          rel->ap[j].e = 1;
          j ++;
        }
    }
  if(j == 0)
      fprintf(stderr, "WARNING: j_ap = 0 in reduce_exponents_mod2\n");
  else
      rel->ap = (alg_prime_t*) realloc (rel->ap, j * sizeof (alg_prime_t));
  rel->nb_ap = j;
}

void
insertNormalRelation(int **rel_compact, int irel, int *nprimes, tab_prime_t *bad_primes, hashtable_t *H, relation_t *rel)
{
    int *tmp, itmp, j, h;

    reduce_exponents_mod2(rel);
    computeroots(rel);
    tmp = (int *)malloc((rel->nb_rp + rel->nb_ap + 1) * sizeof(int));
    itmp = 0;
    for(j = 0; j < rel->nb_rp; j++){
	h = hashInsert(H, rel->rp[j].p, -2);
	tmp[itmp++] = h;
	if(H->hashcount[h] == 1)
	    // new prime
	    *nprimes += 1;
	}
    for(j = 0; j < rel->nb_ap; j++){
	h = hashInsert(H, rel->ap[j].p, rel->ap[j].r);
	tmp[itmp++] = h;
	if(H->hashcount[h] == 1){
	    // new prime
	    *nprimes += 1;
	    if(rel->ap[j].r == -1){
		// bad prime
		if (bad_primes->allocated == bad_primes->length) {
		    bad_primes->allocated += 100;
		    bad_primes->tab = (unsigned long *)realloc((void *)bad_primes->tab, (bad_primes->allocated)*sizeof(unsigned long));
		    assert (bad_primes->tab != NULL);
		}
		bad_primes->tab[bad_primes->length] = rel->ap[j].p;
		bad_primes->length++;
	    }
	}
    }
    tmp[itmp] = -1; // sentinel
    rel_compact[irel] = tmp;
}

// The information is stored in the ap[].p part, which is odd, but convenient.
// rel->ap.p[0..deg[ contains the deg roots of f(x) mod p, leading to deg
// ideals for the factorization of (p).
void
insertFreeRelation(int **rel_compact, int irel, int *nprimes, hashtable_t *H, relation_t *rel)
{
    long p = rel->a; // rel->b == 0
    int j, *tmp, itmp, h;

    tmp = (int *)malloc((1 + rel->nb_ap + 1) * sizeof(int));
    itmp = 0;
    h = hashInsert(H, p, -2);
    tmp[itmp++] = h;
    if(H->hashcount[h] == 1)
	// new prime
	*nprimes += 1;
    for(j = 0; j < rel->nb_ap; j++){
	h = hashInsert(H, p, rel->ap[j].p);
	tmp[itmp++] = h;
	if(H->hashcount[h] == 1)
	    // new prime
	    *nprimes += 1;
    }
    tmp[itmp] = -1; // sentinel
    rel_compact[irel] = tmp;
}

// nrel can decrease, for instance in case of duplicate rows.
int scan_relations_from_file(int *irel, int *nrel, char *rel_used, int **rel_compact, int *nprimes, tab_prime_t *bad_primes, hashtable_t *H, hashtable_t *Hab, FILE *file)
{
    relation_t rel;
    int ret, hab;

    while(1){
	ret = fread_relation (file, &rel);
	if(ret != 1)
	    break;
	*irel += 1;
	if(!(*irel % 100000))
	    fprintf(stderr, "irel = %d\n", *irel);
	hab = hashInsert(Hab, rel.a, rel.b);
	if(Hab->hashcount[hab] > 1){
	    fprintf(stderr,"(%ld, %ld) appears more than once?\n",rel.a,rel.b);
	    rel_used[*irel] = 0;
	    rel_compact[*irel] = NULL;
	    *nrel -= 1;
	    continue;
	}
	rel_used[*irel] = 1;
	if(rel.b > 0)
	    insertNormalRelation(rel_compact,*irel,nprimes,bad_primes,H,&rel);
	else
	    insertFreeRelation(rel_compact, *irel, nprimes, H, &rel);
    }
    return ret;
}

// Read all relations from file.
int
scan_relations(char *ficname[], int nbfic, int *nrel, int *nprimes, tab_prime_t *bad_primes, hashtable_t *H, hashtable_t *Hab, char *rel_used, int **rel_compact)
{
    FILE *relfile;
    int ret, i, irel = -1;
    
    *nprimes = 0;

    for(i = 0; i < nbfic; i++){
	relfile = fopen(ficname[i], "r");
	if(relfile == NULL){
	    fprintf(stderr, "Pb opening file %s\n", ficname[i]);
	    exit(1);
	}
	fprintf(stderr, "Adding file %s\n", ficname[i]);
	ret = scan_relations_from_file(&irel, nrel, rel_used, rel_compact, nprimes, bad_primes, H, Hab, relfile);
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

void delete_relation(int i, int hs, int *nprimes, hashtable_t *H, char *rel_used, int **rel_compact)
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

void
onepass_singleton_removal(int nrelmax, int *nrel, int *nprimes, hashtable_t *H, char *rel_used, int **rel_compact)
{
    int i, h;

    for(i = 0; i < nrelmax; i++)
	if(rel_used[i]){
	    h = has_singleton(rel_compact[i], H);
	    if(h >= 0){
#if DEBUG >= 2
		printf("h = %d is single -> (%ld, %ld)\n", 
		       h, H->hashtab_p[h], H->hashtab_r[h]);
#endif
		delete_relation(i, h, nprimes, H, rel_used, rel_compact);
		*nrel -= 1;
	    }
	}
}

void
find(int nrel, int h, char *rel_used, int **rel_compact)
{
    int i;

    for(i = 0; i < nrel; i++)
	if(rel_used[i]){
	    int *tab = rel_compact[i], j;
	    for(j = 0; tab[j] != -1; j++)
		if(tab[j] == h)
		    fprintf(stderr, "GATO_Z: %d\n", i);
	}
}

void
remove_singletons(int *nrel, int nrelmax, int *nprimes, hashtable_t *H, char *rel_used, int **rel_compact)
{
    int old, newnrel = *nrel, newnprimes = *nprimes, i;

    do{
	old = newnrel;
	onepass_singleton_removal(nrelmax, &newnrel, &newnprimes, H, rel_used, rel_compact);
	fprintf(stderr, "new/old = %d/%d; %d\n", newnrel, old, newnprimes);
    } while(newnrel != old);
    // clean empty rows
    for(i = 0; i < nrelmax; i++){
	if(rel_used[i]){
	    if(rel_compact[i][0] == -1){
		fprintf(stderr, "Empty row: %d\n", i);
	    }
	}
    }
    *nrel = newnrel;
    *nprimes = newnprimes;
}

// we locate used primes and do not try to do fancy things as sorting w.r.t.
// weight, since this will probably done later on.
// All rows will be 1 more that needed -> subtract 1 in fprint...!
void
renumber(int *nprimes, tab_prime_t bad_primes, hashtable_t *H)
{
    int i, nb;

    for(i = 0, nb = 1; i < H->hashmod; i++)
	if(isBadPrime(H->hashtab_p[i], bad_primes) || (H->hashcount[i] == 0))
	    H->hashcount[i] = -1;
	else{
#if 1
	    if(H->hashcount[i] <= 1)
		fprintf(stderr, "WARNING: H->hashcount[%d] = %d\n",
			i, H->hashcount[i]);
#else
	    assert(H->hashcount[i] > 1);
#endif
	    if(H->hashcount[i] > 0)
		H->hashcount[i] = nb++;
	}
    nb--;
    fprintf(stderr, "nb = %d\n", nb);
    *nprimes = nb;
}

void
reread(char *ficname[], int nbfic, tab_prime_t bad_primes, hashtable_t *H, char *rel_used, int nrows, int cols)
{
    FILE *file;
    relation_t rel;
    int irel = -1, ret, i, nr = 0;

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
		    fprintf(stderr, "irel = %d\n", irel);
		if(rel_used[irel]){
		    if(rel.b > 0){
			reduce_exponents_mod2 (&rel);
			computeroots (&rel);
		    }
		    fprint_rel_row(stdout, irel, rel, bad_primes, H);
		    nr++;
		    if(nr >= nrows){
			ret = 0;
			break;
		    }
		}
	    }
	} while(ret == 1);
	fclose(file);
    }
}

int
main (int argc, char **argv)
{
    tab_prime_t bad_primes;
    hashtable_t H, Hab;
    FILE *relfile;
    char *rel_used;
    int **rel_compact;
    int ret;
    int i, nrelmax, nrel, nprimes, nrel_new, nprimes_new;
    
    fprintf (stderr, "%s revision %s\n", argv[0], REV);
    
    if (argc == 1) {
	relfile = stdin;
	fprintf(stderr, "usage: %s [filename]\n", argv[0]);
	fprintf(stderr, "  stdin input is not yet available, sorry.\n");
	exit(1);
    } 
    if (argc < 3) {
	fprintf(stderr, "usage: %s nrel file1 ... filen\n", argv[0]);
	fprintf(stderr, "  if no filename is given, takes input on stdin\n");
	exit(1);
    }
    nrelmax = atoi(argv[1]);
    
    fprintf(stderr, "initializing hash tables...\n");
    hashInit(&H, nrelmax);
    hashInit(&Hab, nrelmax);

    rel_used = (char *)malloc(nrelmax * sizeof(char));
    rel_compact = (int **)malloc(nrelmax * sizeof(int *));
    
    bad_primes.allocated = 100;
    bad_primes.length = 0;
    bad_primes.tab = (unsigned long *)malloc(bad_primes.allocated*sizeof(unsigned long));

    fprintf(stderr, "reading file of relations...\n");
    nrel = nrelmax;
    ret = scan_relations(argv+2, argc-2, &nrel, &nprimes, &bad_primes, &H, &Hab, rel_used, rel_compact);
    assert (ret);
    
    fprintf(stderr, "nrel(useful)=%d, nprimes=%d\n", nrel, nprimes);
    
#if DEBUG >= 2
    for(i = 0; i < hashmod; i++)
	if(hashcount[i] == 1)
	    printf("H[%ld, %ld] = 1\n", hashtab_p[i], hashtab_r[i]);
#endif
    
    fprintf(stderr, "starting singleton removal...\n");
    nrel_new = nrel;
    nprimes_new = nprimes;
    remove_singletons(&nrel_new, nrelmax, &nprimes_new, &H, rel_used, rel_compact);
    
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
    renumber(&nprimes_new, bad_primes, &H);
    
#if DEBUG >= 2
    hashCheck(nprimes_new);
#endif
    
#if 0
    // TODO: do more fancy things, like killing too heavy rows, etc.
    // if just doing this, then singletons may appear again, etc.
    if((nrel_new - nprimes_new) >= 500)
	nrel_new = nprimes_new + 500;
#endif

    // now, we reread the file of relations and convert it to the new coding...
    printf("%d %d\n", nrel_new, nprimes_new);
    reread(argv+2, argc-2, bad_primes, &H, rel_used, nrel_new, nprimes_new);
    
    free(bad_primes.tab);
    free(rel_used);
    for(i = 0; i < nrelmax; i++)
	if(rel_compact[i] != NULL)
	    free(rel_compact[i]);
    free(rel_compact);
    
    return 0;
}
