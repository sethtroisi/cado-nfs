#ifndef CADO_LINALG_HASH_PAIR_H_
#define CADO_LINALG_HASH_PAIR_H_

#ifdef  __cplusplus
extern "C" { 
#endif

typedef struct {
    unsigned long hashmod;
    int *hashcount;
    unsigned long *hashtab_p, *hashtab_r;
} hashtable_t;

extern void hashClear(hashtable_t *H);
extern void hashInit(hashtable_t *H, int n);
extern void hashFree(hashtable_t *H);
extern int getHashAddr(hashtable_t *H, unsigned long p, unsigned long r);
extern int hashInsert(hashtable_t *H, unsigned long p, unsigned long r);

#ifdef	__cplusplus
}	/* extern "C" */
#endif

#endif	/* CADO_LINALG_HASH_PAIR_H_ */
