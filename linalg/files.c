#include <stdlib.h>
#include <stdio.h>
#include "utils/utils.h"

#include "files.h"

// TODO: make this more efficient...
void
getRelation(relation_t *rel, FILE *relfile, int ind)
{
    char str[1024];
    int i;

    rewind(relfile);
    for(i = 0; i < ind; i++)
	fgets(str, 1024, relfile);
    fread_relation(relfile, rel);
}

// We are ready to read relation old and we want new >= old.
void
jumpToRelation(relation_t *rel, FILE *relfile, int old, int new)
{
    char str[1024];
    int i;

    for(i = old; i < new; i++)
	fgets(str, 1024, relfile);
    fread_relation(relfile, rel);
}

void
getAB(long *a, unsigned long *b, FILE *relfile, int ind)
{
    relation_t rel;

    getRelation(&rel, relfile, ind);
    *a = rel.a;
    *b = rel.b;
}
