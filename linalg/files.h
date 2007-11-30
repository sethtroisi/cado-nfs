typedef struct{
    FILE *relfile, *purgedfile, *indexfile, *kerfile;
} files_t;

extern void getRelation(relation_t *rel, FILE *relfile, int ind);
extern void jumpToRelation(relation_t *rel, FILE *relfile, int old, int new);
extern void getAB(long *a, unsigned long *b, FILE *relfile, int ind);
