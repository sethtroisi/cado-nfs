#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>
#include <string.h>

void printVec(int *vec, int ncols)
{
    int i;

    fprintf(stderr, "vec=");
    for(i = 0; i < ncols; i++)
	fprintf(stderr, "%d", vec[i] & 1);
    fprintf(stderr, "\n");
}

int
checkSparse(FILE *matfile, FILE *kerfile, int ncols, int nlimbs, int *vec, int compact, int verbose)
{
    long a;
    unsigned long w, b;
    int i, j, ret, nc, cc, jj;
    char c;

    memset(vec, 0, ncols * sizeof(int));
    // get next dependance relation
    rewind(matfile);
    fscanf(matfile, "%d %d", &i, &j);
    while((c = getc(matfile)) != '\n');
    for(i = 0; i < nlimbs; ++i){
	ret = fscanf(kerfile, "%lx", &w);
	if(ret == -1)
	    break;
	assert (ret == 1);
	if(verbose)
	    fprintf(stderr, "w=%lx\n", w);
	for(j = 0; j < GMP_NUMB_BITS; ++j){
	    if(w & 1UL){
		if(verbose)
		    fprintf(stderr, "+R_%d\n", (i * GMP_NUMB_BITS)+j);
		if(compact == 0)
		    fscanf(matfile, "%ld %lu %d", &a, &b, &nc);
		else
		    fscanf(matfile, "%d", &nc);
		for(jj = 0; jj < nc; jj++){
		    fscanf(matfile, "%d", &cc);
		    if(verbose >= 2)
			fprintf(stderr, "vec[%d]++\n", cc);
		    vec[cc]++;
		}
		if(verbose >= 2)
		    printVec(vec, ncols);
		while((c = getc(matfile)) != '\n')
		    if(feof(matfile))
			break;
	    }
	    else{
		while((c = getc(matfile)) != '\n')
		    if(feof(matfile))
			break;
	    }
	    w >>= 1;
	}
    }
    return ret;
}

void
checkVector(int *vec, int ncols)
{
    int ok, i;

    ok = 1;
    for(i = 0; i < ncols; i++)
	if(vec[i] & 1){
	    ok = 0;
	    break;
	}
    if(ok)
	fprintf(stderr, "y");
    else
	fprintf(stderr, "[n (%d)]", i);
}

void
checkSparseAll(char *matname, char *kername, int compact, int verbose)
{
    FILE *matfile, *kerfile;
    int nrows, ncols, nlimbs, *vec, ndep, ret;

    matfile = fopen(matname, "r");
    kerfile = fopen(kername, "r");
    fscanf(matfile, "%d %d", &nrows, &ncols);
    nlimbs = (nrows / GMP_NUMB_BITS) + 1;
    vec = (int *)malloc(ncols * sizeof(int));
    ndep = 0;
    while(1){
	if(verbose)
	    fprintf(stderr, "Checking dep %d\n", ndep++);
	ret = checkSparse(matfile,kerfile,ncols,nlimbs,vec,compact,verbose);
	if(ret == -1)
	    break;
	if(verbose >= 2)
	    printVec(vec, ncols);
	checkVector(vec, ncols);
    }
    fprintf(stderr, "\n");
    free(vec);
    fclose(matfile);
    fclose(kerfile);
}

int
checkWithIndex(FILE *purgedfile, FILE *indexfile, FILE *kerfile, int verbose, int nrows, int ncols, int small_nrows, int small_ncols, int *vec)
{
    unsigned long w;
    int i, j, k, ret, nlimbs, ind, nrel, r, nc;
    char *small_row_used, *rel_used;

    nlimbs = (small_nrows / GMP_NUMB_BITS) + 1;
    // first read used rows in the small matrix
    small_row_used = (char *)malloc(small_nrows * sizeof(char));
    memset(small_row_used, 0, small_nrows * sizeof(char));
    for(i = 0; i < nlimbs; ++i){
	ret = fscanf(kerfile, "%lx", &w);
	if(ret == -1){
	    free(small_row_used);
	    return ret;
	}
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
    // now map to the rels of the big matrix
    rel_used = (char *)malloc(nrows * sizeof(char));
    memset(rel_used, 0, nrows * sizeof(char));
    rewind(indexfile);
    fscanf(indexfile, "%d %d", &i, &j); // skip first line
    for(i = 0; i < small_nrows; i++){
	fscanf(indexfile, "%d", &nrel);
	for(j = 0; j < nrel; j++){
	    fscanf(indexfile, "%d", &r);
	    if(small_row_used[i])
		rel_used[r] ^= 1;
	}
    }
    // now really read the big matrix in
    rewind(purgedfile);
    fscanf(purgedfile, "%d %d", &i, &j);
    memset(vec, 0, ncols * sizeof(int));
    for(i = 0; i < nrows; i++){
	fscanf(purgedfile, "%d %d", &j, &nc);
	for(j = 0; j < nc; j++){
	    fscanf(purgedfile, "%d", &k);
	    if(rel_used[i])
		vec[k]++;
	}
    }
    return 1;
}

void
checkWithIndexAll(char *purgedname, char *indexname, char *kername, int verbose)
{
    FILE *purgedfile = fopen(purgedname, "r");
    FILE *indexfile = fopen(indexname, "r");
    FILE *kerfile = fopen(kername, "r");
    int nrows, ncols, *vec, small_nrows, small_ncols, ret, ndep;

    fscanf(purgedfile, "%d %d", &nrows, &ncols);
    fscanf(indexfile, "%d %d", &small_nrows, &small_ncols);
    vec = (int *)malloc(ncols * sizeof(int));
    ndep = 0;
    while(1){
        if(verbose)
            fprintf(stderr, "Checking dep %d\n", ndep++);
        ret = checkWithIndex(purgedfile,indexfile,kerfile,verbose,nrows,ncols,small_nrows, small_ncols, vec);
	if(ret == -1)
	    break;
	checkVector(vec, ncols);
    }
    fprintf(stderr, "\n");
    free(vec);
    fclose(kerfile);
    fclose(purgedfile);
    fclose(indexfile);
}

int main(int argc, char *argv[])
{
    int verbose = 1;

#if 0
    int compact = 0;
    if(argc > 3){
	compact = atoi(argv[3]);
	if(argc > 4)
	    verbose = atoi(argv[4]);
    }
    checkSparseAll(argv[1], argv[2], compact, verbose);
#else
    verbose = 0;
    checkWithIndexAll(argv[1], argv[2], argv[3], verbose);
#endif
    return 0;
}
