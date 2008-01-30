#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>
#include <string.h>

#include "cado.h"
#include "utils/utils.h"
#include "files.h"
#include "hashpair.h"

#define DEBUG 0
#define MAPLE 0

#define WANT_ASSERT

int
checkVector(int *vec, int ncols)
{
    int ok, i;

    ok = 1;
    for(i = 0; i < ncols; i++)
	if(vec[i] & 1){
	    ok = 0;
	    break;
	}
#if 0
    if(ok)
	printf(" -> y\n");
    else
	printf(" -> n (%d)\n", i);
#endif
    return ok;
}

void
str2Vec(int *vec, char *str)
{
    char *t = str;
    int k = 0;

    // skip first integer which is the index
    for(; *t != ' '; t++);
    // skip second integer which is the number of primes
    for(++t; *t != ' '; t++);
    t++;
    while(1){
	if((*t == '\n') || (*t == ' ')){
	    // new integer read
	    vec[k]++;
	    k = 0;
	    if(*t == '\n')
		break;
	}
	else
	    k = k*10+(*t - '0');
	t++;
    }
}

void
treatRationalRelation(hashtable_t *H, relation_t rel)
{
    int j, h;

#if MAPLE >= 1
    fprintf(stderr, "P:=P * (%ld-m*%lu):\n", rel.a, rel.b);
#endif
    for(j = 0; j < rel.nb_rp; j++){
#if MAPLE >= 1
	fprintf(stderr, "R2:=R2*%ld^%d:\n", rel.rp[j].p, rel.rp[j].e);
#endif
	h = getHashAddr(H, rel.rp[j].p, -2);
	if(H->hashcount[h] == 0){
	    // new empty place
	    H->hashtab_p[h] = rel.rp[j].p;
	    H->hashtab_r[h] = (unsigned long)(-2);
	}
	H->hashcount[h] += rel.rp[j].e;
    }
}

void
computefree(relation_t *rel)
{
    int i;

    for(i = 0; i < rel->nb_ap; i++){
	rel->ap[i].r = rel->ap[i].p;
	rel->ap[i].p = rel->a;
	rel->ap[i].e = 1;
    }
}

void
treatAlgebraicRelation(FILE *algfile, hashtable_t *H, relation_t rel)
{
    int j, h;

    if(rel.b > 0)
	computeroots(&rel);
    else
	computefree(&rel);
    fprintf(algfile, "%ld %lu\n", rel.a, rel.b);
#if MAPLE >= 1
    fprintf(stderr, "NORM:=NORM*mynorm(f, %ld, %lu):\n", rel.a, rel.b);
    fprintf(stderr, "AX:=AX*(%ld-%lu*X):\n", rel.a, rel.b);
#endif
    for(j = 0; j < rel.nb_ap; j++){
#if MAPLE >= 2
	fprintf(stderr, "A2:=A2*[%ld, %ld]^%d:\n", 
		rel.ap[j].p, rel.ap[j].r, rel.ap[j].e);
#endif
	h = getHashAddr(H, rel.ap[j].p, rel.ap[j].r);
	if(H->hashcount[h] == 0){
	    // new empty place
	    H->hashtab_p[h] = rel.ap[j].p;
	    H->hashtab_r[h] = rel.ap[j].r;
	}
	H->hashcount[h] += rel.ap[j].e;
    }
}

void
finishRationalSqrt(FILE *ratfile, hashtable_t *H, cado_poly pol)
{
    mpz_t prod;
    int j;
    unsigned int i;

    mpz_init_set_ui(prod, 1);
#if MAPLE >= 1
    fprintf(stderr, "R:=1;\n");
#endif
    for(i = 0; i < H->hashmod; i++){
	if(H->hashcount[i] > 0){
	    if(H->hashtab_r[i] != ((unsigned long)-2))
		continue;
	    ASSERT(!(H->hashcount[i] & 1));
#if MAPLE >= 1
	    fprintf(stderr, "R:=R*%ld^%d:\n",
		    H->hashtab_p[i], H->hashcount[i]>>1);
#endif
	    // TODO: do better
	    for(j = 0; j < (H->hashcount[i]>>1); j++)
		mpz_mul_ui(prod, prod, H->hashtab_p[i]);
	    mpz_mod(prod, prod, pol->n);
	}
    }
#if DEBUG >= 1
    gmp_fprintf(stderr, "prod:=%Zd;\n", prod);
    fprintf(stderr, "# We print the squareroot of the rational side...\n");
#endif
    // TODO: humf, I should say...!
    gmp_fprintf(ratfile, "%Zd\n", prod);
    mpz_clear(prod);
}

static void
finishAlgebraicSqrt(FILE *algfile, hashtable_t *H /*, cado_poly pol*/)
{
    unsigned int i;

    fprintf(algfile, "0 0\n");
#if DEBUG >= 1
    fprintf(stderr, "# Now, we print the factors of sqrt(norm): p r e\n");
#endif
    for(i = 0; i < H->hashmod; i++)
	if(H->hashcount[i] > 0){
	    ASSERT(!(H->hashcount[i] & 1));
	    if(H->hashtab_r[i] == ((unsigned long)-2))
		continue;
	    fprintf(algfile, "%ld %ld %d\n", 
		    H->hashtab_p[i], H->hashtab_r[i], H->hashcount[i]>>1);
#if DEBUG >= 1
	    fprintf(stderr, "# H[%d] = %d\n", i, H->hashcount[i]);
#endif
	}
    fprintf(algfile, "0 0\n");
}

// returns the sign of m1*a+m2*b
int
treatSign(relation_t rel, cado_poly pol)
{
    mpz_t tmp1, tmp2;
    int s;

    mpz_init(tmp1);
    mpz_mul_si(tmp1, pol->g[1], rel.a);
    mpz_init(tmp2);
    mpz_mul_ui(tmp2, pol->g[0], rel.b);
    mpz_add(tmp1, tmp1, tmp2);
    s = (mpz_sgn(tmp1) >= 0 ? 1 : -1);
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    return s;
}

// RETURN VALUE: -1 if file exhausted.
//                0 if product is negative...
//                1 if product is positive...
int
treatDep(char *ratname, char *algname, FILE *relfile, FILE *purgedfile, FILE *indexfile, FILE *kerfile, cado_poly pol, int nrows, int ncols, char *small_row_used, int small_nrows, hashtable_t *H, int nlimbs, char *rel_used, int *vec, int rora, int verbose)
{
    FILE *ratfile, *algfile;
    relation_t rel;
    unsigned long w;
    int ret, i, j, nrel, r, irel, nr, sg, ind;
    char str[1024];

    memset(small_row_used, 0, small_nrows * sizeof(char));
    // now use this dep
    for(i = 0; i < nlimbs; ++i){
	ret = fscanf(kerfile, "%lx", &w);
	if(ret == -1)
	    return ret;
	ASSERT (ret == 1);
	if(verbose)
	    fprintf(stderr, "w=%lx\n", w);
	for(j = 0; j < GMP_NUMB_BITS; ++j){
	    if(w & 1UL){
		ind = (i * GMP_NUMB_BITS)+j;
		if(verbose)
		    fprintf(stderr, "+R_%d\n", ind);
		small_row_used[ind] = 1;
	    }
	    w >>= 1;
	}
    }
    // now map to the rels of the purged matrix
    memset(rel_used, 0, nrows * sizeof(char));
    rewind(indexfile);
    fscanf(indexfile, "%d %d", &i, &j); // skip first line
    for(i = 0; i < small_nrows; i++){
	fscanf(indexfile, "%d", &nrel);
	for(j = 0; j < nrel; j++){
	    fscanf(indexfile, "%d", &r);
	    if(small_row_used[i]){
#if DEBUG >= 1
		fprintf(stderr, "# Small[%d] -> %d\n", i, r);
#endif
#if DEBUG >= 1
		if(rel_used[r])
		    fprintf(stderr, "WARNING: flipping rel_used[%d]\n", r);
#endif
		rel_used[r] ^= 1;
	    }
	}
    }
#if MAPLE >= 1
    if((rora == 1) || (rora == 3))
	fprintf(stderr, "R2:=1; P:=1;\n");
    if((rora == 2) || (rora == 3)){
	fprintf(stderr, "A2:=1;\n");
	fprintf(stderr, "AX:=1;\n");
	fprintf(stderr, "NORM:=1;\n");
    }
#endif
    memset(vec, 0, ncols * sizeof(int));
    // FIXME: sg should be removed...
    sg = 1;
    ratfile = fopen(ratname, "w");
    algfile = fopen(algname, "w");
    // now really read the purged matrix in
    rewind(purgedfile);
    fgets(str, 1024, purgedfile); // skip first line
    // we assume purgedfile is stored in increasing order of the indices
    // of the real relations, so that one pass in the rels file is needed...!
    rewind(relfile);
    irel = 0; // we are ready to read relation irel
    for(i = 0; i < nrows; i++){
	fgets(str, 1024, purgedfile);
	if(rel_used[i]){
	    sscanf(str, "%d", &nr);
#if DEBUG >= 1
	    fprintf(stderr, "Reading in rel %d of index %d\n", i, nr);
#endif
	    str2Vec(vec, str);
	    jumpToRelation(&rel, relfile, irel, nr);
	    irel = nr+1;
	    if((rora == 1) || (rora == 3))
		treatRationalRelation(H, rel);
	    if((rora == 2) || (rora == 3))
		treatAlgebraicRelation(algfile, H, rel);
	    sg *= treatSign(rel, pol);
	    clear_relation(&rel);
	}
    }
    ASSERT(checkVector(vec, ncols));
    if(sg == -1)
	fprintf(stderr, "prod(a-b*m) < 0\n");
    else{
	if((rora == 2) || (rora == 3))
          finishAlgebraicSqrt(algfile, H /*, pol*/);
	if((rora == 1) || (rora == 3))
	    finishRationalSqrt(ratfile, H, pol);
    }
    fclose(ratfile);
    fclose(algfile);
    return (sg == -1 ? 0 : 1);
}

void
SqrtWithIndexAll(char *prefix, FILE *relfile, FILE *purgedfile, FILE *indexfile, FILE *kerfile, cado_poly pol, int rora, int ndepmin, int ndepmax, int verbose)
{
    char ratname[200], algname[200];
    hashtable_t H;
    unsigned long w;
    int i, j, ret, nlimbs, nrows, ncols, small_nrows, small_ncols;
    char *small_row_used, *rel_used;
    int *vec; // useful to check dependency relation in the purged matrix

    fscanf(purgedfile, "%d %d", &nrows, &ncols);
    fscanf(indexfile, "%d %d", &small_nrows, &small_ncols);

    nlimbs = (small_nrows / GMP_NUMB_BITS) + 1;
    // first read used rows in the small matrix
    small_row_used = (char *)malloc(small_nrows * sizeof(char));
    rel_used = (char *)malloc(nrows * sizeof(char));
    vec = (int *)malloc(ncols * sizeof(int));

    // skip first ndepmin-1 relations
    for(j = 0; j < ndepmin; j++)
	for(i = 0; i < nlimbs; ++i)
	    ret = fscanf(kerfile, "%lx", &w);

    // use a hash table to rebuild P2
    hashInit(&H, nrows);
    while(1){
	if(ndepmin >= ndepmax)
	    break;
	sprintf(ratname, "%s.rat.%03d", prefix, ndepmin);
	sprintf(algname, "%s.alg.%03d", prefix, ndepmin);
	ret = treatDep(ratname, algname, relfile, purgedfile, indexfile, kerfile, pol, nrows, ncols, small_row_used, small_nrows, &H, nlimbs, rel_used, vec, rora, verbose);
	if(ret == -1)
	    break;
	fprintf(stderr, "# Treated dependency #%d at %2.2lf\n",
		ndepmin, seconds());
	hashClear(&H);
	ndepmin++;
    }
    hashFree(&H);
    free(small_row_used);
    free(rel_used);
    free(vec);
}

int main(int argc, char *argv[])
{
    char *relname, *purgedname, *indexname, *kername, *polyname;
    FILE *relfile, *purgedfile, *indexfile, *kerfile;
    cado_poly pol;
    int verbose = 1, ndepmin, ndepmax, rora, ret;

    if(argc != 10){
	fprintf(stderr, "Usage: %s relname purgedname indexname", argv[0]);
	fprintf(stderr, " kername polyname ndepmin ndepmax r|a prefix\n");
	fprintf(stderr, "Dependency relation i will be put in files prefix.i\n");
	return 0;
    }

    relname = argv[1];
    purgedname = argv[2];
    indexname = argv[3];
    kername = argv[4];
    polyname = argv[5];
    ndepmin = atoi(argv[6]);
    ndepmax = atoi(argv[7]);
    
    relfile = fopen(relname, "r");
    purgedfile = fopen(purgedname, "r");
    indexfile = fopen(indexname, "r");
    kerfile = fopen(kername, "r");

    ret = read_polynomial(pol, polyname);
    ASSERT (ret);

    if(!strcmp(argv[8], "r"))
	rora = 1;
    else if(!strcmp(argv[8], "a"))
	rora = 2;
    else if(!strcmp(argv[8], "ar") || !strcmp(argv[8], "ra"))
	rora = 3;
    else
      {
        fprintf (stderr, "Error, 8th argument must be r, a, ar or ra\n");
        exit (1);
      }
    verbose = 0;
    SqrtWithIndexAll(argv[9], relfile, purgedfile, indexfile, kerfile, pol, rora, ndepmin, ndepmax, verbose);

    fclose(relfile);
    fclose(purgedfile);
    fclose(indexfile);
    fclose(kerfile);

    return 0;
}
