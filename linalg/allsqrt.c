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
treatAlgebraicRelation(FILE *algfile, hashtable_t *H, relation_t rel)
{
    int j, h;

    computeroots(&rel);
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
    int i, j;

    mpz_init_set_ui(prod, 1);
#if MAPLE >= 1
    fprintf(stderr, "R:=1;\n");
#endif
    for(i = 0; i < H->hashmod; i++){
	if(H->hashcount[i] > 0){
	    if(H->hashtab_r[i] != ((unsigned long)-2))
		continue;
	    assert(!(H->hashcount[i] & 1));
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

void
finishAlgebraicSqrt(FILE *algfile, hashtable_t *H, cado_poly pol)
{
    int i;

    fprintf(algfile, "0 0\n");
#if DEBUG >= 1
    fprintf(stderr, "# Now, we print the factors of sqrt(norm): p r e\n");
#endif
    for(i = 0; i < H->hashmod; i++)
	if(H->hashcount[i] > 0){
	    assert(!(H->hashcount[i] & 1));
	    if(H->hashtab_r[i] == ((unsigned long)-2))
		continue;
	    fprintf(algfile, "%lu %ld %d\n", 
		    H->hashtab_p[i], H->hashtab_r[i], H->hashcount[i]>>1);
#if DEBUG >= 1
	    fprintf(stderr, "# H[%d] = %d\n", i, H->hashcount[i]);
#endif
	}
    fprintf(algfile, "0 0\n");
}

// returns the sign of a-b*m
int
treatSign(relation_t rel, cado_poly pol)
{
    mpz_t m, a;
    int s;

    mpz_init(m);
    mpz_neg(m, pol->g[0]);
    mpz_init_set_si(a, rel.a);
    mpz_mul_ui(m, m, rel.b);
    // s = -1 iff a-b*m < 0
    s = (mpz_cmp(m, a) > 0 ? -1 : 1);
    mpz_clear(a);
    mpz_clear(m);
    return s;
}

int
treatDep(FILE *ratfile, FILE *algfile, FILE *relfile, FILE *purgedfile, FILE *indexfile, FILE *kerfile, cado_poly pol, int nrows, int ncols, char *small_row_used, int small_nrows, hashtable_t *H, int nlimbs, char *rel_used, int *vec, int rora, int verbose)
{
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
	assert (ret == 1);
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
	}
    }
    assert(checkVector(vec, ncols));
    if(sg == -1){
	fprintf(stderr, "prod(a-b*m) < 0\n");
    }
    else{
	if((rora == 2) || (rora == 3))
	    finishAlgebraicSqrt(algfile, H, pol);
	if((rora == 1) || (rora == 3))
	    finishRationalSqrt(ratfile, H, pol);
    }
    return 1;
}

void
SqrtWithIndexAll(char *prefix, FILE *relfile, FILE *purgedfile, FILE *indexfile, FILE *kerfile, cado_poly pol, int rora, int ndep, int verbose)
{
    FILE *ratfile, *algfile;
    char ratname[200], algname[200];
    hashtable_t H;
    unsigned long w;
    int i, j, ret, nlimbs, nrows, ncols, small_nrows, small_ncols;
    char *small_row_used, *rel_used;
    int *vec; // useful to check dependancy relation in the purged matrix

    fscanf(purgedfile, "%d %d", &nrows, &ncols);
    fscanf(indexfile, "%d %d", &small_nrows, &small_ncols);

    nlimbs = (small_nrows / GMP_NUMB_BITS) + 1;
    // first read used rows in the small matrix
    small_row_used = (char *)malloc(small_nrows * sizeof(char));
    rel_used = (char *)malloc(nrows * sizeof(char));
    vec = (int *)malloc(ncols * sizeof(int));

    // skip first ndep-1 relations
    for(j = 0; j < ndep; j++)
	for(i = 0; i < nlimbs; ++i)
	    ret = fscanf(kerfile, "%lx", &w);

    // use a hash table to rebuild P2
    hashInit(&H);
    while(1){
	fprintf(stderr, "# Operating on dependancy #%d\n", ndep);
	sprintf(ratname, "%s.rat.%03d", prefix, ndep);
	ratfile = fopen(ratname, "w");
	sprintf(algname, "%s.alg.%03d", prefix, ndep);
	algfile = fopen(algname, "w");
	ret = treatDep(ratfile, algfile, relfile, purgedfile, indexfile, kerfile, pol, nrows, ncols, small_row_used, small_nrows, &H, nlimbs, rel_used, vec, rora, verbose);
	fclose(ratfile);
	fclose(algfile);
	if(ret == -1)
	    break;
	hashClear(&H);
	ndep++;
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
    int verbose = 1, ndep, rora, ret;

    if(argc != 9){
	fprintf(stderr, "Usage: %s relname purgedname indexname", argv[0]);
	fprintf(stderr, " kername polyname ndep r|a prefix\n");
	fprintf(stderr, "Dependancy relation i will be put in files prefix.i\n");
	return 0;
    }

    relname = argv[1];
    purgedname = argv[2];
    indexname = argv[3];
    kername = argv[4];
    polyname = argv[5];
    ndep = atoi(argv[6]);
    
    relfile = fopen(relname, "r");
    purgedfile = fopen(purgedname, "r");
    indexfile = fopen(indexname, "r");
    kerfile = fopen(kername, "r");

    ret = read_polynomial(pol, polyname);
    assert (ret);

    if(!strcmp(argv[7], "r"))
	rora = 1;
    else if(!strcmp(argv[7], "a"))
	rora = 2;
    else if(!strcmp(argv[7], "ar") || !strcmp(argv[7], "ra"))
	rora = 3;
    verbose = 0;
    SqrtWithIndexAll(argv[8], relfile, purgedfile, indexfile, kerfile, pol, rora, ndep, verbose);

    fclose(relfile);
    fclose(purgedfile);
    fclose(indexfile);
    fclose(kerfile);

    return 0;
}
