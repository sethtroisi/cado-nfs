#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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

char **
extractFic(int *nfic, int *nrelmax, char *input)
{
    FILE *ifile = fopen(input, "r");
    char str[1024];
    char **fic;
    int i, nrel;

    *nfic = 0;
    while(fscanf(ifile, "%d %s", &nrel, str) != EOF)
	*nfic += 1;
    fic = (char **)malloc(*nfic * sizeof(char *));
    i = 0;
    rewind(ifile);
    *nrelmax = 0;
    while(fscanf(ifile, "%d %s", &nrel, str) != EOF){
	fic[i] = (char *)malloc((strlen(str)+1) * sizeof(char));
	strncpy(fic[i], str, strlen(str)+1);
	i++;
	*nrelmax += nrel;
    }
    fclose(ifile);
    return fic;
}
