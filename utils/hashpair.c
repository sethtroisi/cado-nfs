#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#include "hashpair.h"
#include "portability.h"
#include "utils.h"

void
hashClear(hashtable_t *H)
{
  MEMSETZERO(H->hc, H->hm);
  MEMSETZERO(H->ht, H->hm);
}

/* need64 is non-zero if we need 64-bit fields for primes and/or ideals */
void
hashInit (hashtable_t *H, p_r_values_t n, unsigned int verbose)
{
  H->hm = (p_r_values_t) (ulong_nextprime ((unsigned long) ((n>>1)*3)));
  if (verbose) fprintf (stderr, "# hashmod = %lu\n", (unsigned long) H->hm);
  SMALLOC(H->hc, H->hm, "hashInit 1");
  SMALLOC(H->ht, H->hm, "hashInit 2");
  if (verbose) fprintf (stderr, "Hash table:\n"
			"Hash (p,r) uses %zu-bit type, hash index uses %zu-bit type,\n"
			"Allocated hash table of total size %lu MB\n",
			sizeof(index_t) * CHAR_BIT, sizeof(p_r_values_t) * CHAR_BIT,
			(unsigned long) (H->hm * (sizeof(HC_T) + sizeof(ht_t))) >> 20);
  H->hr = NULL;
  hashClear (H);
}

void
hashFree (hashtable_t *H)
{
  SFREE(H->hc);
  SFREE(H->ht);
  SFREE(H->hr);
}

void
hashCheck (hashtable_t *H)
{
    p_r_values_t nb = 0;
    HC_T *p = H->hc, *m = &(H->hc[H->hm]);

    while (p != m) if (*p++) nb++;
    fprintf(stderr, "Hash table %1.1f%% full (%lu entries vs size %lu)\n",
            100.0 * (double) nb / (double) H->hm, (unsigned long) nb, (unsigned long) H->hm);
}

/* special values of the root r:
   p   for a projective root (can happen both in normal relations when p
                              divides b, or in free relations when p divides
                              the leading coefficient of f)
   p+2 for rational ideals */
p_r_values_t
hashInsertWithKey(hashtable_t *H, index_t p, index_t r, p_r_values_t h, unsigned int *f)
{
  for ( ; ; ) {
    if (H->ht[h].p == p && H->ht[h].r == r) {
      if (H->hc[h] < UMAX(H->hc[h])) (H->hc[h])++;
      *f = 0;
      return h;
    }
    if (!(H->hc[h])) {
      H->ht[h] = (ht_t) { p, r };
      H->hc[h] = 1;
      *f = 1;
      return h;
    }
    h++;
    if (h >= H->hm) h = 0;
  }
}

p_r_values_t
hashInsert(hashtable_t *H, index_t p, index_t r, unsigned int *f)
{
  return HASHINSERT(H,p,r,f);
}

