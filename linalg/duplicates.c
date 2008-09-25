/* 
 * Program: filter
 * Author : F. Morain
 * Purpose: removing duplicate relations
 * 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <gmp.h>

#include "utils/utils.h"

#include "hashpair.h"
#include "files.h"
#include "gzip.h"

#define STR_LEN_MAX 1024

#define AGRESSIVE_LAMBDA 10

/* If AGRESSIVE_MODE is 0, we use a hash table with open adressing
   and linear probing, each entry containing a 64-bit value identifying
   uniquely the corresponding (a,b) pair, at least for |a|, |b| < 2^53. */
#define AGRESSIVE_MODE 0
typedef struct {
  unsigned long hashmod;
  uint64_t     *hashtab;
} Hashtable_t;

/* If AGRESSIVE_MODE is 1, we use a bit-table of length AGRESSIVE_LAMBDA.
   For AGRESSIVE_LAMBDA=10, the proportion of false duplicates is about 5%.
   For AGRESSIVE_LAMBDA=50, it is about 1%.
*/
#if AGRESSIVE_MODE == 1
typedef struct {
    unsigned long hashmod;
    int len;
    unsigned long *tab;
} smallhash_t;

/* number of bits in an unsigned long */
#define NBITS (sizeof(unsigned long) << 3)
#endif

// TODO: read line by line and do not parse/print relations
// but only print strbuf when new.

// stolen from relation.c with the same semantics
int
fread_buf(char str[STR_LEN_MAX], FILE *file)
{
  int c, i;
    
  // skip spaces and commented lines
  do {
    // skip spaces
    do {
      c = fgetc(file);
      if (c == EOF)
	return -1;
    } while (isspace(c));
    // skip commented lines
    if (c == '#') {
      do {
	c = fgetc(file);
	if (c == EOF)
	  return -1;
      } while (c != '\n');
    } else {
      ungetc(c, file);
      break;
    }
  } while (1);
  
  // copy line into str
  i = 0;
  do {
    c = fgetc(file);
    if (c == EOF)
      return -1;
    str[i++] = c;
    if (i == STR_LEN_MAX) {
      fprintf(stderr, "warning: line too long\n");
      return 0;
    }
  } while (c != '\n');

  str[i] = '\0';
  return 1;
}

// No comment.
void
get_ab(long *a, unsigned long *b, char str[STR_LEN_MAX])
{
    int ia, ib;

    for(ia = 0; str[ia] != ','; ia++);
    str[ia] = ' ';
    for(ib = ia+1; str[ib] != ':'; ib++);
    str[ib] = ' ';
    sscanf(str, "%ld %lu", a, b);
    str[ia] = ',';
    str[ib] = ':';
}

#if AGRESSIVE_MODE == 0
/* Search the first entry in the hash table that is either 0 or H,
   starting from index h. If H was not already in the hash table, add it.
   Return non-zero if H was not already in the table.
*/
int
is_ab_new (Hashtable_t *Hab, unsigned int h, uint64_t H)
{
  while (Hab->hashtab[h] != 0 && Hab->hashtab[h] != H)
    {
      h ++;
      if (UNLIKELY(h == Hab->hashmod))
        h = 0;
    }
  if (Hab->hashtab[h] == 0)
    {
      Hab->hashtab[h] = H;
      return 1;
    }
  else
    return 0;
}
#else
int
is_ab_new(smallhash_t *Hab, long a, unsigned long b)
{
    int hab = getInitialAddressMod((unsigned long)a, b, Hab->HC0, Hab->HC1, Hab->hashmod);
    int i0, i1;

    i0 = hab / NBITS;
    i1 = hab - (i0 * NBITS);
    if(Hab->tab[i0] & ((1UL)<<i1))
	return 0;
    Hab->tab[i0] |= ((1UL)<<i1);
    return 1;
}
#endif

/* irel is the total number of relations read */
// We implement a trick suggested by PZ.
#if AGRESSIVE_MODE == 0
int
remove_duplicates_from_file(FILE *out, int *irel, unsigned int *nrels, Hashtable_t *Hab, int slice, int slice0, FILE *file, int verbose)
#else
int
remove_duplicates_from_file(FILE *out, int *irel, unsigned int *nrels, smallhash_t *Hab, int slice, int slice0, FILE *file, int verbose)
#endif
{
    int ret;
    long a;
    unsigned long b;
    unsigned long file_duplicates = 0;         /* duplicates in this file */
    static unsigned long total_duplicates = 0; /* duplicates in all files */
    char str[STR_LEN_MAX];
    uint64_t HAB;
    unsigned int hab, mask = (((unsigned)1)<<slice)-1;

    while(1){
	ret = fread_buf(str, file);
	if(ret != 1)
	    break;
	*irel += 1;
	if(!(*irel % 100000))
	    fprintf(stderr, "nrel = %d fdup = %lu at %2.2lf\n",
		    *irel, file_duplicates, seconds());
	get_ab(&a, &b, str);
        HAB = getInitialAddress ((uint64_t) a, (uint64_t) b,
                                 (uint64_t) MAGIC_HC0, (uint64_t) MAGIC_HC1);
	hab = HAB % Hab->hashmod;
	if((slice > 0) && ((hab & mask) != (unsigned)slice0))
	    continue;
	if (is_ab_new(Hab, hab, HAB)){
	    *nrels += 1;
	    fprintf(out, "%s", str);
	}
	else{
	    if (file_duplicates ++ < 10 && verbose)
              fprintf(stderr, "(%ld, %lu) appears more than once\n", a, b);
	}
    }
    total_duplicates += file_duplicates;
    fprintf(stderr, "Found %lu duplicates in this file (total %lu, %1.0f%%)\n",
	    file_duplicates, total_duplicates,
            100.0 * (double) total_duplicates / (double) *irel);
    return ret;
}

// Read all relations from file.
#if AGRESSIVE_MODE == 0
int
remove_duplicates(FILE *out, char *ficname[], int nbfic, unsigned int *nrels, Hashtable_t *Hab, int slice, int slice0, int verbose)
#else
int
remove_duplicates(FILE *out, char *ficname[], int nbfic, unsigned int *nrels, smallhash_t *Hab, int slice, int slice0, int verbose)
#endif
{
    FILE *relfile;
    int ret = 0, irel = -1;
    int i;
    
    ASSERT(nbfic > 0);
    for(i = 0; i < nbfic; i++){
	relfile = fopen(ficname[i], "r");
	if(relfile == NULL){
	    fprintf(stderr, "Pb opening file %s\n", ficname[i]);
	    exit(1);
	}
	fprintf(stderr, "Adding file %s\n", ficname[i]);
	ret = remove_duplicates_from_file(out, &irel, nrels, Hab, slice, slice0, relfile, verbose);
	if(ret == 0) {
	    fprintf(stderr, "Warning: error when reading file %s\n", ficname[i]);
	    break;
	}
	fclose(relfile);
    }
    fprintf(stderr, "Scanned %d relations\n", irel+1);
    
    return (ret == -1);
}

int
main(int argc, char **argv)
{
    FILE *outfile = NULL;
#if AGRESSIVE_MODE == 0
    Hashtable_t Hab;
#else
    smallhash_t Hab;
#endif
    char **fic;
    char *outname = NULL;
    unsigned int nfic;
    int ret, k, slice = 0, slice0;
    unsigned int nrelsmax = 0, nrels, Hsize;
    int verbose = 0;
    
    if(argc == 1) {
	fprintf(stderr, "usage: %s [filename]\n", argv[0]);
	fprintf(stderr, "  stdin input is not yet available, sorry.\n");
	exit(1);
    } 
    if(argc < 4) {
	fprintf(stderr, "usage: %s nrel file1 ... filen\n", argv[0]);
	fprintf(stderr, "  if no filename is given, takes input on stdin\n");
	exit(1);
    }

    fprintf (stderr, "%s.r%s", argv[0], REV);
    for (k = 1; k < argc; k++)
      fprintf (stderr, " %s", argv[k]);
    fprintf (stderr, "\n");
    
    while(argc > 1 && argv[1][0] == '-'){
	if(argc > 1 && strcmp(argv[1], "-v") == 0){
          verbose ++;
          argc -= 1;
          argv += 1;
	}
	if(argc > 2 && strcmp(argv[1], "-nrels") == 0){
	    nrelsmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	if(argc > 2 && strcmp(argv[1], "-slice") == 0){
	    slice = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	if(argc > 2 && strcmp(argv[1], "-out") == 0){
	    outname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
    }

    if(nrelsmax == 0)
      {
        fprintf(stderr, "Error, missing -nrels ... option\n");
        exit(1);
      }

    if(slice > 0)
	Hsize = nrelsmax >> slice;
    else
	Hsize = nrelsmax;

    fic = argv+1;
    nfic = argc-1;
    //	fic = extractFic(&nfic, &nrelsmax, argv[3]);
    fprintf(stderr, "Number of relations is %u\n", nrelsmax);
#if AGRESSIVE_MODE == 0
    Hab.hashmod = getHashMod (Hsize + Hsize / 10, 1);
    Hab.hashtab = (uint64_t*) malloc (Hab.hashmod * sizeof(uint64_t));
    memset (Hab.hashtab, 0, Hab.hashmod * sizeof(uint64_t));
#else
    fprintf(stderr, "AGRESSIVE_MODE used\n");
    Hab.hashmod = getHashMod(((unsigned long)nrelsmax) * AGRESSIVE_LAMBDA, 1);
    Hab.len = 1 + (Hab.hashmod / NBITS);
    Hab.tab = (unsigned long *)malloc(Hab.len * sizeof(unsigned long));
    memset(Hab.tab, 0, Hab.len * sizeof(unsigned long));
#endif

    if (outname != NULL)
      outfile = gzip_open(outname, "w");
    else
      outfile = stdout;
    fprintf(stderr, "reading files of relations...\n");
    nrels = 0;
    for(slice0 = 0; slice0 < (1<<slice); slice0++){
	if(slice0)
#if AGRESSIVE_MODE == 0
          memset (Hab.hashtab, 0, Hab.hashmod * sizeof(uint64_t));
#else
          hashClear(&Hab);
#endif
	if(slice)
	    fprintf(stderr, "Performing slice0=%d\n", slice0);
	ret = remove_duplicates(outfile, fic, nfic, &nrels, &Hab, slice,
                                slice0, verbose);
    }
    if (outname != NULL)
      gzip_close(outfile, outname);
    fprintf(stderr, "Number of relations left: %u\n", nrels);

#if AGRESSIVE_MODE == 0
    free (Hab.hashtab);
#else
    hashFree (&Hab);
#endif

    return 0;
}
