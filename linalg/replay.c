/* 
 * Program: history
 * Author : F. Morain
 * Purpose: managing history of merges
 * 
 * Algorithm:
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define DEBUG 0

// newrows[i] contains a new row formed of old rows that correspond to
// true original relations (and not multirelations).
//
// After computing newrows, we deduce for all old rows the list of newrows
// containing it.

void
printRow(FILE *file, int *row)
{
    int i;

    fprintf(file, "[%d]", row[0]);
    for(i = 1; i <= row[0]; i++)
	fprintf(file, " %d", row[i]);
}

void
printOldRows(int **oldrows, int nrows)
{
    int i;
    
    // where oldrows intervene
    for(i = 0; i < nrows; i++)
	if(oldrows[i][0] != 0){
	    fprintf(stderr, "old row_%d takes part to new row(s)", i);
	    printRow(stderr, oldrows[i]);
	    fprintf(stderr, "\n");
	}
}

void
printBuf(FILE *file, int *buf, int ibuf)
{
    int i;

    fprintf(file, "buf[%d] =", ibuf);
    for(i = 0; i < ibuf; i++)
	fprintf(file, " %d", buf[i]);
}

void
removeWeight(int **sparsemat, int *colweight, int i)
{
    int k;

    for(k = 1; k <= sparsemat[i][0]; k++)
	colweight[sparsemat[i][k]]--;
}

// add buf[0..ibuf[ to row i of sparsemat
void
addrel(int **sparsemat, int *colweight, int *buf, int ibuf, int i)
{
    int *rowi = sparsemat[i], *tmp, tmp_len, k1, k2, k;

    if(rowi == NULL){
	// first time, dump buf
	rowi = (int *)malloc((ibuf+1) * sizeof(int));
	rowi[0] = ibuf;
	memcpy(rowi+1, buf, ibuf * sizeof(int));
	sparsemat[i] = rowi;
	for(k = 1; k <= rowi[0]; k++)
	    colweight[rowi[k]]++;
    }
    else{
#if DEBUG >= 1
	printRow(stderr, rowi); fprintf(stderr, "\n");
	printBuf(stderr, buf, ibuf); fprintf(stderr, "\n");
#endif
	// rowi is bound to disappear
	removeWeight(sparsemat, colweight, i);
	tmp_len = rowi[0] + ibuf + 1;
	tmp = (int *)malloc(tmp_len * sizeof(int));
	k1 = 1;
	k2 = 0;
	k = 1;
	while((k1 <= rowi[0]) && (k2 < ibuf)){
	    if(rowi[k1] < buf[k2]){
		tmp[k++] = rowi[k1++];
		colweight[rowi[k1-1]]++;
	    }
	    else if(rowi[k1] > buf[k2]){
		tmp[k++] = buf[k2++];
		colweight[buf[k2-1]]++;
	    }
	    else{
		k1++; k2++;
	    }
	}
	for(; k1 <= rowi[0]; k1++){
	    tmp[k++] = rowi[k1];
	    colweight[rowi[k1]]++;
	}
	for(; k2 < ibuf; k2++){
	    tmp[k++] = buf[k2];
	    colweight[buf[k2]]++;
	}
	tmp[0] = k-1;
	free(rowi);
	sparsemat[i] = tmp;
#if DEBUG >= 1
	printRow(stderr, tmp); fprintf(stderr, "\n");
	if(tmp[0] != (tmp_len-1))
	    printf("#W# shorter length\n"); // who cares, really?
#endif
    }
}

static int 
cmp(const void *p, const void *q) {
    int x = *((int *)p);
    int y = *((int *)q);
    return (x <= y ? -1 : 1);
}

void
makeSparse(int **sparsemat, int *colweight, FILE *purgedfile, int nrows, int ncols, int **oldrows)
{
    int i, j, nj, *buf, buf_len, ibuf, ind, k;

    buf_len = 100;
    buf = (int *)malloc(buf_len * sizeof(int));
    fprintf(stderr, "Reading and treating relations from purged file\n");
    rewind(purgedfile);
    fscanf(purgedfile, "%d %d", &i, &j); // skip first line; check?
    ind = 0;
    while(fscanf(purgedfile, "%d %d", &i, &nj) != EOF){
	if(!(ind % 1000))
	    fprintf(stderr, "Treating old rel #%d\n", ind);
	// store primes in rel
	if(nj > buf_len){
	    fprintf(stderr, "WARNING: realloc for buf [nj=%d]\n", nj);
	    buf = (int *)realloc(buf, nj * sizeof(int));
	    buf_len = nj;
	}
	// take everybody, clean later...!
	ibuf = 0;
	for(i = 0; i < nj; i++){
	    fscanf(purgedfile, "%d", buf+ibuf);
	    ibuf++;
	}
	qsort(buf, ibuf, sizeof(int), cmp);
	// now, we "add" relation [a, b] in all new relations in which it
	// participates
	for(k = 1; k <= oldrows[ind][0]; k++)
	    addrel(sparsemat, colweight, buf, ibuf, oldrows[ind][k]);
	ind++;
    }
}

void
flushSparse(char *sparsename, int **sparsemat, int small_nrows, int small_ncols, int *code)
{
    FILE *ofile = fopen(sparsename, "w");
    unsigned long W = 0;
    int i, j;

    fprintf(ofile, "%d %d\n", small_nrows, small_ncols);
    for(i = 0; i < small_nrows; i++){
	W += sparsemat[i][0];
	fprintf(ofile, "%d", sparsemat[i][0]);
	for(j = 1; j <= sparsemat[i][0]; j++){
#if DEBUG >= 1
	    assert(code[sparsemat[i][j]] > 0);
#endif
	    fprintf(ofile, " %d", code[sparsemat[i][j]]-1); // FIXME
	}
	fprintf(ofile, "\n");
    }
    fclose(ofile);
    fprintf(stderr, "# Weight(M_small) = %lu\n", W);
}

// Do we have enough space for this approach?
void
makeabFile(char *abname, char *purgedname, int nrows, int **newrows, int small_nrows, int small_ncols)
{
    FILE *purgedfile = fopen(purgedname, "r"), *abfile;
    long *ta = (long *)malloc(nrows * sizeof(long));
    unsigned long *tb = (unsigned long *)malloc(nrows * sizeof(unsigned long));
    int i, j, ind;
    char c;

    fprintf(stderr, "Fetching all ab pairs from file %s\n", purgedname);
    fscanf(purgedfile, "%d %d", &i, &j); // skip first line; check?
    ind = 0;
    while(fscanf(purgedfile, "%ld %lu", ta+ind, tb+ind) != EOF){
	if(!(ind % 1000))
	    fprintf(stderr, "Treating rel #%d\n", ind);
	// skip till end of line
	while((c = getc(purgedfile)) != '\n');
	ind++;
    }
    fclose(purgedfile);
    fprintf(stderr, "Building ab file\n");
    abfile = fopen(abname, "w");
    fprintf(abfile, "%d %d\nmultiab\n", small_nrows, small_ncols);
    for(i = 0; i < nrows; i++)
	if(newrows[i] != NULL){
	    fprintf(abfile, "%d\n", newrows[i][0]);
	    for(j = 1; j <= newrows[i][0]; j++)
		fprintf(abfile, "%ld %lu\n", ta[newrows[i][j]], tb[newrows[i][j]]);
	}
    fclose(abfile);
}

// Do we have enough space for this approach?
void
makeIndexFile(char *indexname, int nrows, int **newrows, int small_nrows, int small_ncols)
{
    FILE *indexfile;
    int i, j;

    indexfile = fopen(indexname, "w");
    fprintf(indexfile, "%d %d\n", small_nrows, small_ncols);
    for(i = 0; i < nrows; i++)
	if(newrows[i] != NULL){
	    fprintf(indexfile, "%d", newrows[i][0]);
	    for(j = 1; j <= newrows[i][0]; j++)
		fprintf(indexfile, " %d", newrows[i][j]);
	    fprintf(indexfile, "\n");
	}
    fclose(indexfile);
}

// i1 += i2
// A row is row[0..max] where row[0] = max and the real components are
// row[1..max].
void
addrows(int **sparsemat, int i1, int i2)
{
    int k1, k2, k, len, *tmp, *tmp2;

    assert(sparsemat[i1] != NULL);
    assert(sparsemat[i2] != NULL);
#if DEBUG >= 1
    fprintf(stderr, "R[%d] =", i1); printRow(stderr, sparsemat[i1]); 
    fprintf(stderr, "\n");
    fprintf(stderr, "R[%d] =", i2); printRow(stderr, sparsemat[i2]);
    fprintf(stderr, "\n");
#endif
    len = sparsemat[i1][0] + sparsemat[i2][0] + 1;
    tmp = (int *)malloc(len * sizeof(tmp));
    k = k1 = k2 = 1;

    // loop while everybody is here
    while((k1 <= sparsemat[i1][0]) && (k2 <= sparsemat[i2][0])){
	if(sparsemat[i1][k1] < sparsemat[i2][k2])
	    tmp[k++] = sparsemat[i1][k1++];
	else if(sparsemat[i1][k1] > sparsemat[i2][k2])
            tmp[k++] = sparsemat[i2][k2++];
	else{
	    fprintf(stderr, "WARNING: j1=j2=%d in addrows\n", k1);
	    k1++;
	    k2++;
	}
    }
    // finish with k1
    for( ; k1 <= sparsemat[i1][0]; k1++)
	tmp[k++] = sparsemat[i1][k1];
    // finish with k2
    for( ; k2 <= sparsemat[i2][0]; k2++)
	tmp[k++] = sparsemat[i2][k2];
    assert(k <= len);
    // copy back
    free(sparsemat[i1]);
    tmp2 = (int *)malloc(k * sizeof(int));
    memcpy(tmp2, tmp, k * sizeof(int));
    tmp2[0] = k-1;
    sparsemat[i1] = tmp2;
    free(tmp);
#if DEBUG >= 1
    fprintf(stderr, "row[%d]+row[%d] =", i1, i2);
    printRow(stderr, sparsemat[i1]); fprintf(stderr, "\n");
#endif
}

// on input, colweight[j] contains the weight; on exit, colweight[j]
// contains the new index for j.
void
renumber(int *small_ncols, int *colweight, int ncols)
{
    int j, k, nb, *tmp;

    tmp = (int *)malloc((ncols<<1) * sizeof(int));
    memset(tmp, 0, (ncols<<1) * sizeof(int));
    for(j = 0, nb = 0; j < ncols; j++)
	if(colweight[j] > 0){
	    tmp[nb++] = colweight[j];
	    tmp[nb++] = j;
	}
    *small_ncols = nb>>1;
    qsort(tmp, nb>>1, 2*sizeof(int), cmp);
    memset(colweight, 0, ncols * sizeof(int));
    for(j = 0, k = 1; j < nb; j += 2)
	colweight[tmp[j+1]] = k++; // always this +1 trick
    free(tmp);
}

// a line is "i i1 ... ik", row[i] is to be added to row i1 and destroyed
// at the end of the process...
void
doAllAdds(int **newrows, char *str)
{
    char *t = str;
    int i, ii;

    for(t = str, i = 0; *t != ' '; t++)
	i = 10 * i + (*t - '0');
#if DEBUG >= 1
    fprintf(stderr, "first i is %d\n", i);
#endif
    ++t;
    ii = 0;
    while(1){
	if((*t == '\n') || (*t == ' ')){
#if DEBUG >= 1
	    fprintf(stderr, "next ii is %d\n", ii);
#endif
	    addrows(newrows, ii, i);
	    ii = 0;
	}
	else
	    ii = 10 * ii + (*t - '0');
	if(*t == '\n')
	    break;
	t++;
    }
    // destroy initial row!
    free(newrows[i]);
    newrows[i] = NULL;
}

// We start from M_purged which is nrows x ncols;
// we build M_small which is small_nrows x small_ncols.
// newrows[i] if != NULL, contains a list of the indices of the rows in
// M_purged that were added together to form this new row in M_small.
// TODO: replace this index by the index to rels directly to skip one
// indirection???
int
main(int argc, char *argv[])
{
    FILE *hisfile, *purgedfile;
    int nrows, ncols;
    int **newrows, i, j, nb, *nbrels, **oldrows, *colweight;
    int ind, small_nrows, small_ncols, **sparsemat;
    char str[1024];

    if(argc != 5){
	fprintf(stderr, "Usage: %s purgedfile hisfile sparsefile abfile\n", argv[0]);
	return 0;
    }
    
    purgedfile = fopen(argv[1], "r");
    assert(purgedfile != NULL);
    // read parameters that should be the same as in purgedfile!
    hisfile = fopen(argv[2], "r");
    assert(hisfile != NULL);
    fgets(str, 1024, hisfile);
    sscanf(str, "%d %d", &nrows, &ncols);
    newrows = (int **)malloc(nrows * sizeof(int *));
    for(i = 0; i < nrows; i++){
	newrows[i] = (int *)malloc(2 * sizeof(int));
	newrows[i][0] = 1;
	newrows[i][1] = i;
    }
    fprintf(stderr, "Reading row additions\n");
    while(fgets(str, 1024, hisfile))
	doAllAdds(newrows, str);
    fclose(hisfile);
    nbrels = (int *)malloc(nrows * sizeof(int));
    memset(nbrels, 0, nrows * sizeof(int));
    // nbrels[oldi] contains the number of new relations in which old
    // relation of index oldi takes part
    small_nrows = 0;
    for(i = 0; i < nrows; i++)
	if(newrows[i] != NULL){
	    small_nrows++;
#if DEBUG >= 1
	    fprintf(stderr, "New row %d:", small_nrows-1);
#endif
	    for(j = 1; j <= newrows[i][0]; j++){
#if DEBUG >= 1
		fprintf(stderr, " %d", newrows[i][j]);
#endif
		nbrels[newrows[i][j]] += 1;
	    }
#if DEBUG >= 1
	    fprintf(stderr, "\n");
#endif
	}
    fprintf(stderr, "Allocating oldrows\n");
    oldrows = (int **)malloc(nrows * sizeof(int *));
    for(i = 0; i < nrows; i++){
	oldrows[i] = (int *)malloc((nbrels[i]+1) * sizeof(int));
	oldrows[i][0] = 0;
    }
    fprintf(stderr, "Filling oldrows\n");
    for(i = 0, nb = 0; i < nrows; i++)
        if(newrows[i] != NULL){
	    // this is row of index nb in the new matrix
	    for(j = 1; j <= newrows[i][0]; j++){
		ind = newrows[i][j];
		oldrows[ind][0]++;
		oldrows[ind][oldrows[ind][0]] = nb;
	    }
	    nb++;
	}
#if DEBUG >= 1
    printOldRows(oldrows, nrows);
#endif
    colweight = (int *)malloc(ncols * sizeof(int *));
    memset(colweight, 0, ncols * sizeof(int *));
    fprintf(stderr, "Building sparse representation\n");
    sparsemat = (int **)malloc(small_nrows * sizeof(int *));
    for(i = 0; i < small_nrows; i++)
	sparsemat[i] = NULL;
    makeSparse(sparsemat, colweight, purgedfile, nrows, ncols, oldrows);
    fclose(purgedfile);

    fprintf(stderr, "Renumbering columns (including sorting w.r.t. weight)\n");
    renumber(&small_ncols, colweight, ncols);

    fprintf(stderr, "small_nrows=%d small_ncols=%d\n",small_nrows,small_ncols);

    fprintf(stderr, "Writing sparse representation to file\n");
    flushSparse(argv[3], sparsemat, small_nrows, small_ncols, colweight);

#if 0
    makeabFile(argv[4], argv[1], nrows, newrows, small_nrows, small_ncols);
#else
    makeIndexFile(argv[4], nrows, newrows, small_nrows, small_ncols);
#endif

    for(i = 0; i < small_nrows; i++)
	if(sparsemat[i] != NULL)
	    free(sparsemat[i]);
    free(sparsemat);

    for(i = 0; i < nrows; i++){
	if(newrows[i] != NULL)
	    free(newrows[i]);
	if(oldrows[i] != NULL)
	    free(oldrows[i]);
    }
    free(newrows);
    free(oldrows);
    free(colweight);
    return 0;
}
