#ifndef CADO_LINALG_HASH_PAIR_H_
#define CADO_LINALG_HASH_PAIR_H_

#ifdef  __cplusplus
extern "C" { 
#endif

typedef struct {
    unsigned long hashmod;
    int *hashcount;
    long *hashtab_p;
    unsigned long *hashtab_r;
} hashtable_t;

#define HC0 314159265358979323UL
#define HC1 271828182845904523UL

#define getInitialAddress(a, b, M) ((unsigned int)((HC0*(a)+HC1*(b)) % M))

extern unsigned long getHashMod(unsigned long n);
extern void hashClear(hashtable_t *H);
extern void hashInit(hashtable_t *H, unsigned int n);
extern void hashFree(hashtable_t *H);
extern int getHashAddr(hashtable_t *H, long p, unsigned long r);
extern int hashInsert(hashtable_t *H, long p, unsigned long r);

#ifdef	__cplusplus
}	/* extern "C" */
#endif

#endif	/* CADO_LINALG_HASH_PAIR_H_ */
