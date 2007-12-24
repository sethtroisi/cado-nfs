#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "hashpair.h"

#define HC0 314159265358979323UL
#define HC1 271828182845904523UL

void
hashClear(hashtable_t *H)
{
    memset(H->hashcount, 0, H->hashmod * sizeof(int));
    memset(H->hashtab_p, 0, H->hashmod * sizeof(long));
    memset(H->hashtab_r, 0, H->hashmod * sizeof(unsigned long));
}

void
hashInit(hashtable_t *H, int n)
{
    int ntab = 13, i;
    unsigned long tab[] = {500009, 750019, 1000003, 
			   5000011, 7500013, 10000019, 
			   50000017, 
			   100000007, 200000033, 300000007, 400000009,
			   500000003, 1000000007};

    n = (3*n)/2;
    for(i = 0; i < ntab; i++)
	if(n < tab[i]){
	    H->hashmod = tab[i];
	    break;
	}
    //    ASSERT(i <= ntab);
    fprintf(stderr, "# hashmod[%d] = %lu\n", i, H->hashmod);
    H->hashcount = (int *)malloc(H->hashmod * sizeof(int));
    H->hashtab_p = (long *)malloc(H->hashmod * sizeof(long));
    H->hashtab_r = (unsigned long *)malloc(H->hashmod * sizeof(unsigned long));
    hashClear(H);
}

void
hashFree(hashtable_t *H)
{
    free(H->hashcount);
    free(H->hashtab_p);
    free(H->hashtab_r);
}

// Returns a new address or the one already containing (p, r).
int
getHashAddr(hashtable_t *H, long p, unsigned long r)
{
    unsigned long tmp = HC0 * ((unsigned long)p) + HC1 * r;
    static int cptmax = 0;
    int h, h0;

    h = h0 = (int)(tmp % H->hashmod);
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
	if(h >= H->hashmod)
	    h = 0;
    }
#if DEBUG >= 1
    if((h-h0) >= 50)
	printf("#W# cpt = %d >= 50\n", h-h0);
#endif
    if((h-h0) >= cptmax){
	cptmax = h-h0;
	if(cptmax > 100) fprintf(stderr, "HASH: cptmax = %d\n", cptmax);
    }
    return h;
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

void hashCheck(hashtable_t *H, int nprimes)
{
    int i, nb = 0;

    for(i = 0; i < H->hashmod; i++)
	if(H->hashcount[i] > 0)
	    nb++;
    fprintf(stderr, "// #[hash>0]=%d vs. nprimes=%d\n", nb, nprimes);
}

