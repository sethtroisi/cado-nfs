#ifndef DCLIST_H_
#define DCLIST_H_

// doubly chained lists
typedef struct dclist {
    int32_t j;
    struct dclist *prev, *next;
} *dclist;

#ifdef __cplusplus
extern "C" {
#endif

extern dclist dclistCreate(int32_t j);
extern void dclistClear (dclist dcl);
extern void dclistTex(FILE *file, dclist dcl);
extern int dclistLength(dclist dcl);
extern dclist dclistInsert(dclist dcl, int32_t j);
extern void dclistPrint(FILE *file, dclist dcl);
extern void dclistRemove (dclist dcl);
extern int dclistRemoveMessy(dclist * p_dcl);
extern dclist dclistFirst (dclist dcl);
/* consistency check. Protect this presumably with ifndef NDEBUG */
extern void dclistCheck(dclist dcl);

#ifdef __cplusplus
}
#endif

#endif	/* DCLIST_H_ */
