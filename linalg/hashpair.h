typedef struct {
    unsigned long hashmod;
    int *hashcount;
    unsigned long *hashtab_p, *hashtab_r;
} hashtable_t;

extern void hashClear(hashtable_t *H);
extern void hashInit(hashtable_t *H);
extern void hashFree(hashtable_t *H);
extern int getHashAddr(hashtable_t *H, unsigned long p, unsigned long r);
extern int hashInsert(hashtable_t *H, unsigned long p, unsigned long r);



