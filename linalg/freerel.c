/* 
 * Program: free relations
 * Author : F. Morain
 * Purpose: creating free relations in a suitable format
 * 
 */


#include "cado.h"
#include <gmp.h>
#include "mod_ul.c"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "utils/utils.h"

#include <string.h>

#include "sieve/fb.h"

int
countFreeRelations(int *deg, char *roots)
{
    FILE *ifile = fopen(roots, "r");
    char str[1024], str0[128], str1[128], str2[128], *t;
    int nfree = 0, nroots;

    *deg = -1;
    // look for the degree
    while(1){
	if(!fgets(str, 1024, ifile))
	    break;
	sscanf(str, "%s %s %s", str0, str1, str2);
	if((str0[0] == '#') && !strcmp(str1, "DEGREE:")){
	    *deg = atoi(str2);
	    break;
	}
    }
    if(*deg == -1){
	fprintf(stderr, "No degree found, sorry\n");
	exit(1);
    }
    fprintf(stderr, "DEGREE = %d\n", *deg);
    while(1){
	if(!fgets(str, 1024, ifile))
            break;
	// format is "11: 2,6,10"
	sscanf(str, "%s %s", str0, str1);
	nroots = 1;
	for(t = str1; *t != '\0'; t++)
	    if(*t == ',')
		nroots++;
	if(nroots == *deg)
	    nfree++;
    }
    fclose(ifile);
    return nfree;
}

void
addFreeRelations(char *roots, int deg)
{
    factorbase_degn_t *fb = fb_read(roots, 1.0, 0), *fbptr;
    fbprime_t p;
    int i;

    for(fbptr = fb; fbptr->p != 0; fbptr = fb_next(fbptr)){
	p = fbptr->p;
	if(fbptr->nr_roots == deg){
	    // free relation, add it!
	    // (p, 0) -> p-0*m
	    printf("%d,0:%lx:", p, (long)p);
	    for(i = 0; i < fbptr->nr_roots; i++){
		printf("%x", fbptr->roots[i]);
		if(i < fbptr->nr_roots-1)
		    printf(",");
	    }
	    printf("\n");
	}
    }
}

static void
usage (char *argv0)
{
  fprintf (stderr, "Usage: %s -poly xxx.poly -fb xxx.roots\n", argv0);
  exit (1);
}

int
main(int argc, char *argv[])
{
    int nfree, deg;
    char *fbfilename = NULL, *polyfilename = NULL;
    char *argv0 = argv[0];
    cado_poly cpoly;

    while (argc > 1 && argv[1][0] == '-')
      {
        if (argc > 2 && strcmp (argv[1], "-fb") == 0)
          {
            fbfilename = argv[2];
            argc -= 2;
            argv += 2;
          }
        else if (argc > 2 && strcmp (argv[1], "-poly") == 0)
          {
            polyfilename = argv[2];
            argc -= 2;
            argv += 2;
          }
        else
          usage (argv0);
      }

    if (argc != 1 || polyfilename == NULL || fbfilename == NULL)
      usage (argv0);

    if (!read_polynomial (cpoly, polyfilename))
      {
        fprintf (stderr, "Error reading polynomial file\n");
        exit (EXIT_FAILURE);
      }

    if (mpz_cmp_ui (cpoly->g[1], 1) != 0)
      {
        fprintf (stderr, "Error, non-monic linear polynomial not yet treated");
        fprintf (stderr, " (more theory needed)\n");
        exit (1);
      }

    nfree = countFreeRelations (&deg, fbfilename);
    fprintf (stderr, "# Free relations: %d\n", nfree);
    
    fprintf (stderr, "Handling free relations...\n");
    addFreeRelations (fbfilename, deg);
    clear_polynomial (cpoly);

    return 0;
}
