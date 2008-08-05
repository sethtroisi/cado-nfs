#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "hashpair.h"
#include "../utils/utils.h"

void
hashClear(hashtable_t *H)
{
    memset(H->hashcount, 0, H->hashmod * sizeof(int));
    memset(H->hashtab_p, 0, H->hashmod * sizeof(long));
    memset(H->hashtab_r, 0, H->hashmod * sizeof(unsigned long));
}

unsigned long
getHashMod(unsigned long n, int verbose)
{
  n = ulong_nextprime (n);
  if (verbose)
    fprintf (stderr, "# hashmod = %lu\n", n);
  return n;
}

void
hashInit(hashtable_t *H, unsigned int n, int verbose)
{
    H->hashmod = getHashMod((3*n)/2, verbose);
    H->hashcount = (int *)malloc(H->hashmod * sizeof(int));
    H->hashtab_p = (long *)malloc(H->hashmod * sizeof(long));
    H->hashtab_r = (unsigned long *)malloc(H->hashmod * sizeof(unsigned long));
    if (verbose)
      {
        unsigned long alloc_hashcount = H->hashmod * sizeof(int);
        unsigned long alloc_hashtab_p = H->hashmod * sizeof(long);
        unsigned long alloc_hashtab_r = H->hashmod * sizeof(unsigned long);
        fprintf (stderr, "Allocated hash tables of total size %luMb\n",
                 (alloc_hashcount + alloc_hashtab_p + alloc_hashtab_r) / 1000000);
      }
    if(sizeof(unsigned long) == 8){
	H->HC0 = 314159265358979323UL;
	H->HC1 = 271828182845904523UL;
    }
    else{
	H->HC0 = 3141592653UL;
	H->HC1 = 2718281828UL;
    }
    hashClear(H);
}

void
hashFree(hashtable_t *H)
{
    free(H->hashcount);
    free(H->hashtab_p);
    free(H->hashtab_r);
}

// Returns a new address or the one already containing (p, r), starting from
// h.
int
getHashAddrAux(hashtable_t *H, long p, unsigned long r, unsigned int h)
{
    static unsigned int cptmax = 0;
    unsigned int cpt = 0; /* number of iterations to find a free location */

#if DEBUG >= 2
    printf("H(%ld, %lu) = %d\n", p, r, h);
#endif

    while(1){
	if(H->hashcount[h] == 0)
	    break;
	if((H->hashtab_p[h] == p) && (H->hashtab_r[h] == r))
	    // (p, r) already known
	    break;
	h++;
        cpt++;
	if(h >= H->hashmod)
	    h = 0;
    }
#if DEBUG >= 1
    if(cpt >= 50)
	printf("#W# cpt = %u >= 50\n", cpt);
#endif
    if(cpt > cptmax){
	cptmax = cpt;
	if(cptmax > 100) fprintf(stderr, "HASH: cptmax = %u\n", cptmax);
    }
    return h;
}

int
getHashAddr(hashtable_t *H, long p, unsigned long r)
{
    unsigned int h = getInitialAddress((unsigned long)p, r, H->HC0, H->HC1, H->hashmod);

    return getHashAddrAux(H, p, r, h);
}

int
hashInsert(hashtable_t *H, long p, unsigned long r)
{
    int h = getHashAddr(H, p, r);
    if(H->hashcount[h] == 0){
	// new empty place
	H->hashtab_p[h] = p;
	H->hashtab_r[h] = r;
    }
    H->hashcount[h]++;
    return h;
}

void
hashCheck (hashtable_t *H)
{
    int nb = 0;
    unsigned int i;

    for(i = 0; i < H->hashmod; i++)
	if(H->hashcount[i] > 0)
	    nb++;
    fprintf(stderr, "Hash table %1.1f%% full (%u entries vs size %lu)\n",
            100.0 * (double) nb / (double) H->hashmod, nb, H->hashmod);
}

