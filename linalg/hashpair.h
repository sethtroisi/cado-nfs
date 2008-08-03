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

#define getInitialAddress(a, b, HC0, HC1, M) ((unsigned int)(((HC0)*(a)+(HC1)*(b)) % M))

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
