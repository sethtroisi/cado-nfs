#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#include "hashpair.h"
#include "utils.h"

void
hashClear(hashtable_t *H)
{
    memset(H->hashcount, 0, H->hashmod * sizeof(int32_t));
    if (H->need64)
      {
        memset(H->hashtab64_p, 0, H->hashmod * sizeof(int64_t));
        memset(H->hashtab64_r, 0, H->hashmod * sizeof(int64_t));
      }
    else
      {
        memset(H->hashtab32_p, 0, H->hashmod * sizeof(int32_t));
        memset(H->hashtab32_r, 0, H->hashmod * sizeof(int32_t));
      }
}

unsigned long
getHashMod(unsigned long n, int verbose)
{
  n = ulong_nextprime (n);
  if (verbose)
    fprintf (stderr, "# hashmod = %lu\n", n);
  return n;
}

/* need64 is non-zero if we need 64-bit fields for primes and/or ideals */
void
hashInit (hashtable_t *H, unsigned int n, int verbose, int need64)
{
    H->need64 = need64;
    H->size = sizeof(int32_t) /* hashcount */
      + ((need64) ? sizeof(int64_t) : sizeof(int32_t)) /* hashtab_p */
      + ((need64) ? sizeof(uint64_t) : sizeof(uint32_t)); /* hashtab_r */
    H->hashmod = getHashMod (3*(n/2), verbose);
    H->hashcount = (int *)malloc(H->hashmod * sizeof(int));
    if (H->need64)
      {
        H->hashtab64_p = (int64_t*) malloc (H->hashmod * sizeof(int64_t));
        H->hashtab64_r = (uint64_t*) malloc (H->hashmod * sizeof(uint64_t));
      }
    else
      {
        H->hashtab32_p = (int32_t*) malloc (H->hashmod * sizeof(int32_t));
        H->hashtab32_r = (uint32_t*) malloc (H->hashmod * sizeof(uint32_t));
      }
    if (verbose)
      {
        fprintf (stderr, "Using %d-bit types\n",
                 (need64) ? 64 : 32);
        fprintf (stderr, "Allocated hash table of total size %"PRIu64"Mb\n",
                 (H->hashmod * H->size) >> 20);
      }
    if(sizeof(unsigned long) == 8){
      H->HC0 = MAGIC_HC0;
      H->HC1 = MAGIC_HC1;
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
    if (H->need64)
      {
        free (H->hashtab64_p);
        free (H->hashtab64_r);
      }
    else
      {
      free (H->hashtab32_p);
      free (H->hashtab32_r);
      }
}

// Returns a new address or the one already containing (p, r), starting from
// h.
int
getHashAddrAux (hashtable_t *H, long p, unsigned long r, unsigned int h)
{
#if DEBUG >= 2
  printf("H(%ld, %lu) = %d\n", p, r, h);
#endif

  while (1)
    {
      if (H->hashcount[h] == 0)
        break;
      if ((GET_HASH_P(H,h) == p) && (GET_HASH_R(H,h) == r))
        // (p, r) already known
        break;
      h++;
      if(h >= H->hashmod)
        h = 0;
    }
  return h;
}

/* Map (p,r) into an hash address, find the corresponding entry in the hash
   table, and return the corresponding index. */
int
getHashAddr (hashtable_t *H, long p, unsigned long r)
{
    uint64_t mask = (H->need64) ? (uint64_t) (-1) : (uint32_t) (-1);
    unsigned int h;

    r = r & mask; /* ensures that -1 and -2 are mapped to -1 and -2
                     in 32-bit mode */
    h = getInitialAddressMod((unsigned long) p, r, H->HC0, H->HC1, H->hashmod);

    return getHashAddrAux (H, p, r, h);
}

int
hashInsert(hashtable_t *H, long p, unsigned long r)
{
    int h = getHashAddr(H, p, r);
    if(H->hashcount[h] == 0){
	// new empty place
        SET_HASH_P(H,h,p);
        SET_HASH_R(H,h,r);
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
    fprintf(stderr, "Hash table %1.1f%% full (%u entries vs size %"PRIu64")\n",
            100.0 * (double) nb / (double) H->hashmod, nb, H->hashmod);
}

