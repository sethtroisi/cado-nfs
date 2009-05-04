#ifndef CADO_LINALG_HASH_PAIR_H_
#define CADO_LINALG_HASH_PAIR_H_

#include <stdint.h>

typedef struct {
    int need64;
    int size; /* size in bytes for one entry */
    uint64_t hashmod, HC0, HC1;
    int32_t *hashcount;
    int32_t *hashtab32_p;
    int64_t *hashtab64_p;
    uint32_t *hashtab32_r;
    uint64_t *hashtab64_r;
} hashtable_t;

#define GET_HASH_P(H,h) \
  (((H)->need64) ? (H)->hashtab64_p[h] : (H)->hashtab32_p[h])
#define SET_HASH_P(H,h,v)                                       \
  do {                                                          \
    if ((H)->need64)                                            \
      (H)->hashtab64_p[h] = (v);                                \
    else                                                        \
      (H)->hashtab32_p[h] = (v); \
  } while (0)

#define GET_HASH_R(H,h) \
  (((H)->need64) ? (H)->hashtab64_r[h] : (H)->hashtab32_r[h])
#define SET_HASH_R(H,h,v)                                       \
  do {                                                          \
    if ((H)->need64)                                            \
      (H)->hashtab64_r[h] = (v);                                \
    else                                                        \
      (H)->hashtab32_r[h] = (v); \
  } while (0)

/* The constants C0 = floor(Pi*10^17) and C1 = floor(exp(1)*10^17) are given
   in Cavallar's thesis. They guarantee that if a, b have at most 53 bits,
   then H(a,b) = C0*a+C1*b is injective, i.e., has no collision. */
#define MAGIC_HC0 (uint64_t) 314159265358979323ULL
#define MAGIC_HC1 (uint64_t) 271828182845904523ULL

#define getInitialAddress(a,b,HC0,HC1) ((HC0)*((uint64_t) (a))+(HC1)*((uint64_t) (b)))
#define getInitialAddressMod(a,b,HC0,HC1,M) (getInitialAddress(a,b,HC0,HC1) % M)

#ifdef  __cplusplus
extern "C" { 
#endif

extern unsigned long getHashMod(unsigned long n, int verbose);
extern void hashClear(hashtable_t *H);
extern void hashInit(hashtable_t *H, unsigned int n, int verbose, int);
extern void hashFree(hashtable_t *H);
extern int getHashAddrAux(hashtable_t *H, long p, unsigned long r, unsigned int h);
extern int getHashAddr(hashtable_t *H, long p, unsigned long r);
extern int hashInsert(hashtable_t *H, long p, unsigned long r);
extern void hashCheck (hashtable_t *H);

#ifdef	__cplusplus
}	/* extern "C" */
#endif

#endif	/* CADO_LINALG_HASH_PAIR_H_ */
