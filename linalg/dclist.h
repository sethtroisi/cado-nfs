#ifndef DCLIST_H_
#define DCLIST_H_

// doubly chained lists
typedef struct dclist{
    INT j;
    struct dclist *prev, *next;
} *dclist;

#ifdef __cplusplus
extern "C" {
#endif

extern dclist dclistCreate(INT j);
extern void dclistTex(FILE *file, dclist dcl);
extern int dclistLength(dclist dcl);
extern dclist dclistInsert(dclist dcl, INT j);
extern void dclistPrint(FILE *file, dclist dcl);

#ifdef __cplusplus
}
#endif

#endif	/* DCLIST_H_ */
