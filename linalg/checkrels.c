/* 
 * Program: checkrels
 * Author : F. Morain
 * Purpose: checking relations
 * 
 *
 */

#include <gmp.h>
#include "mod_ul.c"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#define WANT_ASSERT

#include "utils/utils.h"

#include "files.h"

void
checkRelations(char **fic, int nbfic, cado_poly pol)
{
    FILE *ifile;
    relation_t rel;
    int ret, i, j, irel = 0, nfic;
    mpz_t m1, m2, a, tmp, bm2;

    // let's write a*m1+b*m2
    mpz_init(m1);
    mpz_init(m2);
    mpz_init(a);
    mpz_init(tmp);
    mpz_init(bm2);
    mpz_set(m1, pol->g[1]);
    mpz_set(m2, pol->g[0]);
    gmp_fprintf(stdout, "m1:=%Zd; m2:=%Zd\n", m1, m2);
    for(nfic = 0; nfic < nbfic; nfic++){
	fprintf(stderr, "Checking relations from file %s\n", fic[nfic]);
	ifile = fopen(fic[nfic], "r");
	while(1){
	    irel++;
	    if(!(irel % 100000))
		fprintf(stderr, "Treating %d-th relation at %2.2lf\n",
			irel, seconds());
	    ret = fread_relation (ifile, &rel);
	    if(ret != 1)
		break;
	    for(i = 0; i < rel.nb_rp; i++)
		for(j = i+1; j < rel.nb_rp; j++)
		    if(rel.rp[j].p == rel.rp[i].p)
			fprintf(stderr,"GASP[%d]: rp[%d] = rp[%d]\n",irel,i,j);
	    for(i = 0; i < rel.nb_ap; i++)
		for(j = i+1; j < rel.nb_ap; j++)
		    if((rel.ap[j].p == rel.ap[i].p) 
		       && (rel.ap[j].r == rel.ap[i].r))
			fprintf(stderr,"GASP[%d]: ap[%d] = ap[%d]\n",irel,i,j);
	    mpz_mul_si(tmp, m1, rel.a);
	    mpz_mul_ui(bm2, m2, rel.b);
	    mpz_add(tmp, tmp, bm2);
	    for(i = 0; i < rel.nb_rp; i++){
		if(!mpz_divisible_ui_p(tmp, rel.rp[i].p)){
		    fprintf(stderr, "GASP: m1*a+m2*b not divisible by %ld\n",
			    (long)rel.rp[i].p);
		    printf("a:=%ld;b:=%lu;p:=%ld;\n",rel.a,rel.b,rel.rp[i].p);
		}
	    }
	    clear_relation(&rel);
	}
	fclose(ifile);
    }
    mpz_clear(a);
    mpz_clear(m1);
    mpz_clear(m2);
    mpz_clear(tmp);
    mpz_clear(bm2);
}

int
main(int argc, char *argv[])
{
    char *polyname = NULL;
    cado_poly pol;

    while(argc > 1 && argv[1][0] == '-'){
	if(argc > 2 && strcmp (argv[1], "-poly") == 0){
	    polyname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
    }    
    read_polynomial(pol, polyname);
    checkRelations(argv+1, argc-1, pol);
    return 0;
}
