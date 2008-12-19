// doubly chained lists
typedef struct dclist{
    INT j;
    struct dclist *prev, *next;
} *dclist;

extern dclist dclistCreate(INT j);
extern void dclistTex(FILE *file, dclist dcl);
extern int dclistLength(dclist dcl);
