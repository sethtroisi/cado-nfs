/* 
 * Program: filter
 * Author : F. Morain
 * Purpose: removing duplicate relations
 * 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <gmp.h>
#include "mod_ul.c"

#define WANT_ASSERT

#include "utils/utils.h"

#include "hashpair.h"
#include "files.h"

#define STR_LEN_MAX 1024

// TODO: read line by line and do not parse/print relations
// but only print strbuf when new.

void
fprint_relation_raw(FILE *file, relation_t rel)
{
    int i, j;

    fprintf(file, "%ld,%lu:", rel.a, rel.b);
    for(i = 0; i < rel.nb_rp; ++i){
	for(j = 0; j < rel.rp[i].e; j++){
	    fprintf(file, "%lx", rel.rp[i].p);
	    if(j < rel.rp[i].e-1)
		fprintf(file, ",");
	}
	if(i < rel.nb_rp-1)
	    fprintf(file, ",");
    }
    fprintf(file, ":");
    for(i = 0; i < rel.nb_ap; ++i){
	for(j = 0; j < rel.ap[i].e; j++){
	    fprintf(file, "%lx", rel.ap[i].p);
	    if(j < rel.ap[i].e-1)
                fprintf(file, ",");
        }
        if(i < rel.nb_ap-1)
            fprintf(file, ",");
    }
    fprintf(file, "\n");
}

// stolen from relation.c
int
fread_buf(char str[STR_LEN_MAX], FILE *file)
{
  int c, i;
    
  // skip spaces and commented lines
  do {
    // skip spaces
    do {
      c = fgetc(file);
      if (c == EOF)
	return -1;
    } while (isspace(c));
    // skip commented lines
    if (c == '#') {
      do {
	c = fgetc(file);
	if (c == EOF)
	  return -1;
      } while (c != '\n');
    } else {
      ungetc(c, file);
      break;
    }
  } while (1);
  
  // copy line into str
  i = 0;
  do {
    c = fgetc(file);
    if (c == EOF)
      return -1;
    str[i++] = c;
    if (i == STR_LEN_MAX) {
      fprintf(stderr, "warning: line too long\n");
      return 0;
    }
  } while (c != '\n');

  str[i] = '\0';
  return 1;
}

// No comment.
void
get_ab(long *a, unsigned long *b, char str[STR_LEN_MAX])
{
    int ia, ib;

    for(ia = 0; str[ia] != ','; ia++);
    str[ia] = ' ';
    for(ib = ia+1; str[ib] != ':'; ib++);
    str[ib] = ' ';
    sscanf(str, "%ld %lu", a, b);
    str[ia] = ',';
    str[ib] = ':';
}

int
remove_duplicates_from_file(int *irel, unsigned int *nrels, hashtable_t *Hab, FILE *file)
{
    int ret, hab;
    long a;
    unsigned long b;
    unsigned long file_duplicates = 0;         /* duplicates in this file */
    static unsigned long total_duplicates = 0; /* duplicates in all files */
    char str[STR_LEN_MAX];

    while(1){
	ret = fread_buf(str, file);
	if(ret != 1)
	    break;
	*irel += 1;
	if(!(*irel % 100000))
	    fprintf(stderr, "nrel = %d at %2.2lf\n", *irel, seconds());
	get_ab(&a, &b, str);
	hab = hashInsert(Hab, a, b);
	if(Hab->hashcount[hab] > 1){
	    if(file_duplicates ++ < 10)
		fprintf(stderr, "(%ld, %lu) appears more than once\n", a, b);
	    continue;
	}
	else{
	    *nrels += 1;
	    printf("%s", str);
	}
    }
    total_duplicates += file_duplicates;
    fprintf(stderr, "Found %lu duplicates in this file (total %lu)\n",
	    file_duplicates, total_duplicates);
    return ret;
}

// Read all relations from file.
int
remove_duplicates(char *ficname[], int nbfic, unsigned int *nrels, hashtable_t *Hab)
{
    FILE *relfile;
    int ret = 0, irel = -1;
    int i;
    
    ASSERT(nbfic > 0);
    *nrels = 0;
    for(i = 0; i < nbfic; i++){
	relfile = fopen(ficname[i], "r");
	if(relfile == NULL){
	    fprintf(stderr, "Pb opening file %s\n", ficname[i]);
	    exit(1);
	}
	fprintf(stderr, "Adding file %s\n", ficname[i]);
	ret = remove_duplicates_from_file(&irel, nrels, Hab, relfile);
	if(ret == 0) {
	    fprintf(stderr, "Warning: error when reading file %s\n", ficname[i]);
	    break;
	}
	fclose(relfile);
    }
    fprintf(stderr, "Scanned %d relations\n", irel+1);
    
    return (ret == -1);
}

int
main(int argc, char **argv)
{
    hashtable_t Hab;
    char **fic;
    unsigned int nfic;
    int ret, k;
    unsigned int nrelsmax = 0, nrels;
    
    if(argc == 1) {
	fprintf(stderr, "usage: %s [filename]\n", argv[0]);
	fprintf(stderr, "  stdin input is not yet available, sorry.\n");
	exit(1);
    } 
    if(argc < 4) {
	fprintf(stderr, "usage: %s nrel file1 ... filen\n", argv[0]);
	fprintf(stderr, "  if no filename is given, takes input on stdin\n");
	exit(1);
    }

    fprintf (stderr, "%s.r%s", argv[0], REV);
    for (k = 1; k < argc; k++)
      fprintf (stderr, " %s", argv[k]);
    fprintf (stderr, "\n");
    
    while(argc > 1 && argv[1][0] == '-'){
	if(argc > 2 && strcmp(argv[1], "-nrels") == 0){
	    nrelsmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
    }

    if(nrelsmax == 0)
      {
        fprintf(stderr, "Error, missing -nrels ... option\n");
        exit(1);
      }

    fic = argv+1;
    nfic = argc-1;
    //	fic = extractFic(&nfic, &nrelsmax, argv[3]);
    fprintf(stderr, "Number of relations is %u\n", nrelsmax);
    hashInit(&Hab, nrelsmax);

    fprintf(stderr, "reading files of relations...\n");
    ret = remove_duplicates(fic, nfic, &nrels, &Hab);
    fprintf(stderr, "Number of relations left: %u\n", nrels);

    return 0;
}
