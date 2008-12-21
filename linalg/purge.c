/* purge --- remove singletons

Copyright 2008 Francois Morain, Paul Zimmermann

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

#include <gmp.h>
#include "mod_ul.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "utils/utils.h"

#include "hashpair.h"
#include "files.h"
#include "gzip.h"

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

/*****************************************************************************/

/* dirty trick to distinguish rational primes: we store -2 for their root */
static unsigned long minus2;

/* Adds in table_ind[] the indices of the rational primes in 'rel'.

   All primes are assumed to be different (and to appear with an odd exponent).

   nb_coeff is the current index in table_ind[] (input and output).
 */
static void
fprint_rat (int *table_ind, int *nb_coeff, relation_t rel, hashtable_t *H)
{
    int i, nbc = *nb_coeff;

    for (i = 0; i < rel.nb_rp; i++)
      table_ind[nbc++] = H->hashcount[getHashAddr (H, rel.rp[i].p, minus2)];
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
        ASSERT(H->hashcount[index] >= 0);
        table_ind[nbc++] = H->hashcount[index];
    }
    *nb_coeff = nbc;
}

/* Print the relation 'rel' in matrix format, i.e., a line of the form:

   i k t_1 t_2 ... t_k

   i (decimal) is the row index from the nodup file (starting at 0)
   k (decimal) is the number of rational and algebraic primes in the relation
   t_1 ... t_k (hexadecimal) are the indices of the primes (starting at 0)

   Assumes the remaining primes in 'rel' are those with an odd exponent,
   and are all different.

   WARNING: the primes in the input relation are not necessarily sorted.
*/
static void
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
      fprint_rat (table_ind, &nb_coeff, rel, H);

      /* adds algebraic primes in table_ind */
      fprint_alg (table_ind, &nb_coeff, rel, H);
    }

  fprintf (file, "%d %d", irel, nb_coeff);
  for (i = 0; i < nb_coeff; ++i)
    /* due to the +1 in renumber */
    fprintf (file, " " PURGE_INT_FORMAT, table_ind[i] - 1);
  fprintf (file, "\n");

  free (table_ind);
}

/* this function is used in 'reduce' mode only */
static int
specialHashInsert (hashtable_t *H, long p, unsigned long r, int irel)
{
    int h = getHashAddr(H, p, r);

    if(H->hashcount[h] == 0){
      /* new empty place */
        SET_HASH_P(H,h,p);
        SET_HASH_R(H,h,r);
	H->hashcount[h] = -(irel+1); /* trick!!! */
    }
    else
	// more than second time
      H->hashcount[h] = irel+1; /* overcautious!!! */
    return h;
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
static void
insertNormalRelation (int **rel_compact, int irel,
                      int *nprimes, hashtable_t *H, relation_t *rel,
                      unsigned long minpr, unsigned long minpa,
                      int final, unsigned long *tot_alloc)
{
    int *tmp = NULL, ltmp = 0, i, j, h;

    reduce_exponents_mod2 (rel);
    computeroots (rel); /* FIXME: we might compute roots only for algebraic
                           primes > minpa */
    if (final){
      /* first count number of "large" primes */
	for (j = 0; j < rel->nb_rp; j++)
          ltmp += (rel->rp[j].p >= minpr); /* only consider primes >= minpr */
	for (j = 0; j < rel->nb_ap; j++)
          ltmp += (rel->ap[j].p >= minpa); /* only consider primes >= minpa */

        /* ltmp is the number of considered primes in the relation.
           We might have ltmp=0 if all primes are less than minpr, maxpr. */
        tmp = my_malloc_int (ltmp + 1);
        *tot_alloc += (ltmp + 1) * sizeof (int);
	i = 0; /* number of entries in tmp */

        /* trick: we use the same hash table for rational and algebraic
           primes, but use a fake root -2 for rational primes, which
           ensures there is no collision with algebraic primes */
	for (j = 0; j < rel->nb_rp; j++)
          {
            /* we need to hash-insert also primes < minpr, to get a correct
               count of distinct primes */
            h = hashInsert (H, rel->rp[j].p, minus2);
            *nprimes += (H->hashcount[h] == 1); /* new prime */
	    if (rel->rp[j].p >= minpr)
              tmp[i++] = h;
          }

	for (j = 0; j < rel->nb_ap; j++)
          {
	    h = hashInsert (H, rel->ap[j].p, rel->ap[j].r);
            *nprimes += (H->hashcount[h] == 1); /* new prime */
	    if (rel->ap[j].p >= minpa)
              tmp[i++] = h;
          }

        tmp[i] = -1; /* sentinel */
        rel_compact[irel] = tmp;
    }
    else{
      /* reduce mode only */
	for(j = 0; j < rel->nb_rp; j++){
	    h = specialHashInsert(H, rel->rp[j].p, minus2, irel);
	    if(H->hashcount[h] < 0)
	      /* new prime */
		*nprimes += 1;
	}
	for(j = 0; j < rel->nb_ap; j++){
	    h = specialHashInsert(H, rel->ap[j].p, rel->ap[j].r, irel);
	    if(H->hashcount[h] < 0){
	      /* new prime */
		*nprimes += 1;
	    }
	}
    }
}

/* The information is stored in the ap[].p part, which is odd, but convenient.
   rel->ap.p[0..deg[ contains the deg roots of f(x) mod p, leading to deg
   ideals for the factorization of (p). */
static void
insertFreeRelation (char *rel_used, int **rel_compact, int irel, int *nprimes,
                    hashtable_t *H, relation_t *rel, int final,
                    unsigned long *tot_alloc)
{
    long p = rel->a; /* rel->b == 0 */
    int j, h, *tmp = NULL, itmp;

    if(final){
	tmp = (int *)malloc((1 + rel->nb_ap + 1) * sizeof(int));
        *tot_alloc += (1 + rel->nb_ap + 1) * sizeof(int);
	itmp = 0;
	h = hashInsert(H, p, minus2);
	if(H->hashcount[h] == 1)
	  /* new prime */
	    *nprimes += 1;
	tmp[itmp++] = h;
	for(j = 0; j < rel->nb_ap; j++){
	    h = hashInsert(H, p, rel->ap[j].p);
	    if(H->hashcount[h] == 1)
	      /* new prime */
		*nprimes += 1;
	    tmp[itmp++] = h;
	}
	tmp[itmp++] = -1;
	rel_used[irel] = 1;
	rel_compact[irel] = tmp;
    }
    else{
      /* reduce mode */
	h = specialHashInsert(H, p, minus2, irel);
	if(H->hashcount[h] < 0)
	  /* new prime */
	    *nprimes += 1;
	for(j = 0; j < rel->nb_ap; j++){
	    h = specialHashInsert(H, p, rel->ap[j].p, irel);
	    if(H->hashcount[h] < 0)
	      /* new prime */
		*nprimes += 1;
	}
    }
}

/* nrel can decrease, for instance in case of duplicate rows. */
static int
scan_relations_from_file (int *irel, int *nrel, char *rel_used,
                          int **rel_compact, int *nprimes, hashtable_t *H,
                          FILE *file, long minpr, long minpa,
                          int final, unsigned long *tot_alloc)
{
    relation_t rel;
    int ret;

    while(1){
        /* read rational and algebraic primes and put them in rel
           (identical primes are merged and the corresponding exponent
            computed, but the root for alg. primes is not computed) */
	ret = fread_relation (file, &rel);
	if(ret != 1){
	    break;
	}
	*irel += 1;
	if(!(*irel % 100000))
	    fprintf(stderr, "nrel = %d at %2.2lfs (memory %luMb)\n",
                    *irel, seconds (), *tot_alloc >> 20);
	rel_used[*irel] = 1;
	if(rel.b > 0)
          insertNormalRelation (rel_compact, *irel, nprimes, H, &rel,
                                minpr, minpa, final, tot_alloc);
	else
          insertFreeRelation (rel_used, rel_compact, *irel, nprimes, H, &rel,
                              final, tot_alloc);
	if(rel_used[*irel] <= 0)
	  /* relation removed */
	    *nrel -= 1;
	clear_relation(&rel);
    }
    return ret;
}

/* Read all relations from file, and fills the rel_used and rel_compact arrays
   for each relation i:
   - rel_used[i] = 0 or -1 if the relation has no considered prime
       (i.e., a prime in the interval [minpr, oo] on the rational side),
       otherwise rel_used[i] = 1
   - rel_compact is an array, terminated by -1, of pointers to the entries
     in the hash table for the considered primes
 */
static int
scan_relations (char *ficname[], int nbfic, int *nrel, int *nprimes,
                hashtable_t *H, char *rel_used, int **rel_compact,
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
                                        nprimes, H, relfile, minpr, minpa,
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
    rel_compact[i] = NULL;
    rel_used[i] = 0;
}

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

    if ((*nrel - *nprimes) <= keep)
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
    /* first stupid idea: just remove heavy rows */
    for(i = 0; i < *nrel; i += 2){
	if ((*nrel)-(*nprimes) <= keep)
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
	if(rel_used[i] > 0){
	    int h = has_singleton (rel_compact[i], H);
	    if(h >= 0){
		delete_relation (i, nprimes, H, rel_used, rel_compact);
		*nrel -= 1;
	    }
	}
    }
}

static void
remove_singletons (int *nrel, int nrelmax, int *nprimes, hashtable_t *H,
                   char *rel_used, int **rel_compact, int keep)
{
  int old, newnrel = *nrel, newnprimes = *nprimes;

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

    /* Warning: we might have empty rows, i.e., empty rel_compact[i] lists,
       if all primes in a relation are less then minpr or minpa. */

    *nrel = newnrel;
    *nprimes = newnprimes;
}

/* we locate used primes and do not try to do fancy things as sorting w.r.t.
   weight, since this will probably be done later on.
   All rows will be 1 more that needed -> subtract 1 in fprint...! */
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
      if(H->hashcount[i] == 0)
	    H->hashcount[i] = -1;
	else{
	    ASSERT(H->hashcount[i] > 1);
	    H->hashcount[i] = nb++;
	    if(fsos != NULL)
	      fprintf(fsos, "%d %" PRIx64 " %" PRIx64 "\n",
		      H->hashcount[i] - 1, GET_HASH_P(H,i), GET_HASH_R(H,i));
	}
    if(fsos != NULL)
      gzip_close (fsos, sos);
    nb--;
    *nprimes = nb;
}

static void
reread(FILE *ofile, char *ficname[], unsigned int nbfic,
       hashtable_t *H, char *rel_used, int nrows)
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
		if(rel_used[irel]){ /* covers -1 and > 0 */
		    if(rel.b > 0){
			reduce_exponents_mod2 (&rel);
			computeroots (&rel);
		    }
		    fprint_rel_row (ofile, irel, rel, H);
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

    /* locate singletons */
    for(i = 0; i < H->hashmod; i++)
	if(H->hashcount[i] < 0){
	  rel_used[-H->hashcount[i]-1] = 0; /* that trick again! */
	    newnprimes--;
	}
    /* I/O */
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
  fprintf (stderr, "Usage: %s [options] -poly polyfile -out purgedfile -nrels nnn file1 ... filen\n", argv[0]);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "       -nonfinal    - perform only one singleton pass\n");
  fprintf (stderr, "       -keep    nnn - stop when excess <= nnn (default -1)\n");
  fprintf (stderr, "       -minpa   nnn - purge alg. primes >= nnn (default alim)\n");
  fprintf (stderr, "       -minpr   nnn - purge rat. primes >= nnn (default rlim)\n");
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
main (int argc, char **argv)
{
    FILE *purgedfile = NULL;
    hashtable_t H;
    char **fic, *polyname = NULL, *sos = NULL, *purgedname = NULL;
    unsigned int nfic;
    char *rel_used;
    int **rel_compact = NULL;
    int ret, k;
    int nrel, nprimes = 0, final = 1;
    unsigned int nrelmax = 0;
    int nrel_new, nprimes_new, Hsize, Hsizer, Hsizea;
    long keep = -1; /* maximum value for nrows-ncols */
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
	else if(argc > 2 && strcmp (argv[1], "-minpr") == 0){
	    minpr = atol (argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 2 && strcmp (argv[1], "-minpa") == 0){
	    minpa = atol (argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if(argc > 2 && strcmp (argv[1], "-keep") == 0){
	    keep = atol (argv[2]);
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
	else if(argc > 1 && strcmp (argv[1], "-out") == 0){
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

    fic = argv + 1;
    nfic = argc - 1;
    cado_poly_init(pol);
    cado_poly_read(pol, polyname);

    /* On a 32-bit computer, even 1 << 32 would overflow. Well, we could set
       map[ra] = 2^32-1 in that case, but not sure we want to support 32-bit
       primes on a 32-bit computer... */
    need64 = (pol[0].lpbr >= 32) || (pol[0].lpba >= 32);

    if (need64 && sizeof (long) < 8)
      {
        fprintf (stderr, "Error, too large primes for a 32-bit computer\n");
        exit (1);
      }

    minus2 = ULONG_MAX - 1;

    if (minpr == 0)
      minpr = pol[0].rlim;
    if (minpa == 0)
      minpa = pol[0].alim;

    fprintf (stderr, "Number of relations is %u\n", nrelmax);
    if (nprimes > 0)
	Hsize = nprimes;
    else{
      /* Estimating the number of needed primes (remember that hashInit
         multiplies by a factor 1.5). */
      Hsizer = approx_phi (1L << pol[0].lpbr) - approx_phi (minpr);
      Hsizea = approx_phi (1L << pol[0].lpba) - approx_phi (minpa);
      Hsize = Hsizer + Hsizea;
    }
    fprintf (stderr, "initializing hash tables with Hsize=%d...\n", Hsize);
    hashInit (&H, Hsize, 1, need64);
    tot_alloc = H.hashmod * H.size;

    rel_used = (char *) malloc (nrelmax);
    tot_alloc += nrelmax;
    fprintf (stderr, "Allocated rel_used of %uMb (total %luMb so far)\n",
             nrelmax >> 20,
             tot_alloc >> 20);
    if (final)
      {
        /* FIXME: the rel_compact alone uses a lot of memory. For 100M
           relations (which is reasonable for a 155-digit number with large
           prime bounds 2^30), on a 64-bit machine it uses 0.8Gb!!! */
	rel_compact = (int **) malloc (nrelmax * sizeof (int *));
        tot_alloc += nrelmax * sizeof (int*);
        /* %zu is the C99 modifier for size_t */
        fprintf (stderr, "Allocated rel_compact of %zuMb (total %luMb so far)\n",
                 (nrelmax * sizeof (int *)) >> 20,
                 tot_alloc >> 20);
      }

    fprintf(stderr, "Reading file of relations...\n");
    nrel = nrelmax;
    ret = scan_relations (fic, nfic, &nrel, &nprimes, &H, rel_used,
                          rel_compact, minpr, minpa, final, &tot_alloc);
    ASSERT (ret);

    fprintf (stderr, "nrel(useful)=%d, nprimes=%d (expected %d)\n",
             nrel, nprimes, Hsize);

    if(final == 0)
	reduce (fic, nfic, &H, rel_used, nrelmax, nprimes);
    else{
	fprintf(stderr, "starting singleton removal...\n");
	nrel_new = nrel;
	nprimes_new = nprimes;
	remove_singletons (&nrel_new, nrelmax, &nprimes_new, &H, rel_used,
                           rel_compact, keep);

	fprintf(stderr, "nrel=%d, nprimes=%d; excess=%d\n",
		nrel_new, nprimes_new, nrel_new - nprimes_new);

	if (nrel_new <= nprimes_new) /* covers the case nrel = nprimes = 0 */
	    exit(1);

	/* we renumber the primes in order of apparition in the hashtable */
        fprintf (stderr, "Renumbering primes...\n");
	renumber (&nprimes_new, &H, sos);

        fprintf (stderr, "Freeing rel_compact array...\n");
	/* we do not use it anymore */
        my_malloc_free_all ();
	free(rel_compact);

	/* now, we reread the file of relations and convert it to the new
	   coding... */
        fprintf (stderr, "Storing remaining relations...\n");
	purgedfile = gzip_open (purgedname, "w");
	fprintf (purgedfile, "%d %d\n", nrel_new, nprimes_new);
	reread (purgedfile, fic, nfic, &H, rel_used, nrel_new);
	gzip_close (purgedfile, purgedname);
	/* write excess to stdout */
	printf("EXCESS: %d\n", nrel_new - nprimes_new);
    }
    free(rel_used);

    return 0;
}
