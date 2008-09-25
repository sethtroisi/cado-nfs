#ifndef CADO_LINALG_HASH_PAIR_H_
#define CADO_LINALG_HASH_PAIR_H_

#ifdef  __cplusplus
extern "C" { 
#endif

typedef struct {
    unsigned long hashmod, HC0, HC1;
    int *hashcount;
    long *hashtab_p;
    unsigned long *hashtab_r;
} hashtable_t;

/* The constants C0 = floor(Pi*10^17) and C1 = floor(exp(1)*10^17) are given
   in Cavallar's thesis. They guarantee that if a, b have at most 53 bits,
   then H(a,b) = C0*a+C1*b is injective, i.e., has no collision. */
#define MAGIC_HC0 314159265358979323UL
#define MAGIC_HC1 271828182845904523UL

#define getInitialAddress(a,b,HC0,HC1) ((HC0)*(a)+(HC1)*(b))
#define getInitialAddressMod(a,b,HC0,HC1,M) (getInitialAddress(a,b,HC0,HC1) % M)

extern unsigned long getHashMod(unsigned long n, int verbose);
extern void hashClear(hashtable_t *H);
extern void hashInit(hashtable_t *H, unsigned int n, int verbose);
extern void hashFree(hashtable_t *H);
extern int getHashAddrAux(hashtable_t *H, long p, unsigned long r, unsigned int h);
extern int getHashAddr(hashtable_t *H, long p, unsigned long r);
extern int hashInsert(hashtable_t *H, long p, unsigned long r);
extern void hashCheck (hashtable_t *H);

#ifdef	__cplusplus
}	/* extern "C" */
#endif

#endif	/* CADO_LINALG_HASH_PAIR_H_ */
