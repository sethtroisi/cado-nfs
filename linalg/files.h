#ifndef	CADO_LINALG_FILES_H_
#define	CADO_LINALG_FILES_H_

#include <stdio.h>

#ifdef	__cplusplus
extern "C" {
#endif

#define PURGE_INT_FORMAT "%x"

typedef struct{
    FILE *relfile, *purgedfile, *indexfile, *kerfile;
} files_t;

extern void getRelation(relation_t *rel, FILE *relfile, int ind);
extern void jumpToRelation(relation_t *rel, FILE *relfile, int old, int new);
extern void getAB(long *a, unsigned long *b, FILE *relfile, int ind);

extern char **extractFic(int *nfic, int *nrelmax, char *input);

#ifdef	__cplusplus
}	/* extern "C" */
#endif

#endif	/* CADO_LINALG_FILES_H_ */
