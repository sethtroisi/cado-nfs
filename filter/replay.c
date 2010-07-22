#define _BSD_SOURCE     /* strdup */
/* 
 * Program: replay
 * Author : F. Morain
 * Purpose: replaying history of merges to build the small matrix
 * 
 * Algorithm:
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "files.h"

#include "sparse.h"
#include "gzip.h"

#define DEBUG 0

#define TRACE_COL -1 // 231 // put to -1 if not...!

// newrows[i] contains a new row formed of old rows that correspond to
// true original relations (and not multirelations).
//
// After computing newrows, we deduce for all old rows the list of newrows
// containing it.

#if DEBUG >= 1
static void
printOldRows(int **oldrows, int nrows)
{
    int i;
    
    // where oldrows intervene
    for(i = 0; i < nrows; i++)
	if(oldrows[i][0] != 0){
	    fprintf(stderr, "old row_%d takes part to new row(s)", i);
	    fprintRow(stderr, oldrows[i]);
	    fprintf(stderr, "\n");
	}
}

static void
printBuf(FILE *file, int *buf, int ibuf)
{
    int i;

    fprintf(file, "buf[%d] =", ibuf);
    for(i = 0; i < ibuf; i++)
	fprintf(file, " %d", buf[i]);
}
#endif

// add buf[0..ibuf[ to row i of sparsemat
static void
addrel(int **sparsemat, int *colweight, int *buf, int ibuf, int i)
{
    int *rowi = sparsemat[i], *tmp, tmp_len, k1, k2, k;

    if(rowi == NULL){
	// first time, dump buf
	rowi = (int *)malloc((ibuf+1) * sizeof(int));
	rowi[0] = ibuf;
	memcpy(rowi+1, buf, ibuf * sizeof(int));
	sparsemat[i] = rowi;
	for(k = 1; k <= rowi[0]; k++){
#if TRACE_COL >= 0
	    if(rowi[k] == TRACE_COL)
		fprintf(stderr, "addrel: j=%d appears 1st in R_%d\n",TRACE_COL,i);
#endif
	    colweight[rowi[k]]++;
	}
    }
    else{
#if DEBUG >= 1
	fprintRow(stderr, rowi); fprintf(stderr, "\n");
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
#if TRACE_COL >= 0
	    if(rowi[k1] == TRACE_COL)
		fprintf(stderr, "addrel: j=%d appears in R_%d\n",TRACE_COL,i);
	    if(buf[k2] == TRACE_COL)
		fprintf(stderr, "addrel: j=%d appears in buf_%d\n",TRACE_COL,i);
#endif
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
	fprintRow(stderr, tmp); fprintf(stderr, "\n");
	//	if(tmp[0] != (tmp_len-1))
	//	    printf("#W# shorter length\n"); // who cares, really?
#endif
	//	fprintf(stderr, "W: 231 -> %d\n", colweight[231]);
    }
}

static void
makeSparse(int **sparsemat, int *colweight, FILE *purgedfile,
           int jmin, int jmax, int **oldrows, int verbose, int nslices)
{
    int i, j, nj, *buf, buf_len, ibuf, ind, k;
    int report = (verbose == 0) ? 100000 : 10000;
    int rc;

    buf_len = 100;
    buf = (int *) malloc (buf_len * sizeof(int));
    fprintf(stderr, "Reading and treating relations from purged file\n");
    rc = fscanf(purgedfile, "%d %d", &i, &j); // skip first line; check?
    ASSERT_ALWAYS(rc == 2);

    ind = 0;
    while(fscanf(purgedfile, "%d %d", &i, &nj) != EOF){
	if(!(ind % report))
	    fprintf(stderr, "Treating old rel #%d at %2.2lf\n",ind,seconds());
	// store primes in rel
	if(nj > buf_len){
	    fprintf(stderr, "WARNING: realloc for buf [nj=%d]\n", nj);
	    buf = (int *) realloc (buf, nj * sizeof(int));
	    buf_len = nj;
	}
	// take everybody and clean
	ibuf = 0;
	for(i = 0; i < nj; i++){
	    rc = fscanf(purgedfile, PURGE_INT_FORMAT, buf+ibuf);
            ASSERT_ALWAYS(rc == 1);

	    // what a trick!!!!
	    if((buf[ibuf] >= jmin) && (buf[ibuf] < jmax))
		ibuf++;
	}
	if(ibuf > 0){
	    qsort(buf, ibuf, sizeof(int), cmp);
	    // now, we "add" relation [a, b] in all new relations in which it
	    // participates
	    for(k = 1; k <= oldrows[ind][0]; k++)
		addrel(sparsemat, colweight, buf, ibuf, oldrows[ind][k]);
	}
	if(nslices == 0){
	    // free space just in case
	    free(oldrows[ind]);
	    oldrows[ind] = NULL;
	}
	ind++;
    }
    free (buf);
}

static unsigned long
flushSparse_ascii(const char *sparsename, int **sparsemat, int small_nrows, int small_ncols, int *code)
{
    FILE *ofile;
    unsigned long W = 0;
    int i, j;

    ofile = gzip_open(sparsename, "w");
    fprintf(ofile, "%d %d\n", small_nrows, small_ncols);
    for(i = 0; i < small_nrows; i++){
	if(sparsemat[i] == NULL)
	    fprintf(ofile, "0");
	else{
	    W += sparsemat[i][0];
	    fprintf(ofile, "%d", sparsemat[i][0]);
	    for(j = 1; j <= sparsemat[i][0]; j++){
#if DEBUG >= 1
		ASSERT(code[sparsemat[i][j]] > 0);
#endif
		fprintf(ofile, " %d", code[sparsemat[i][j]]-1); // FIXME
	    }
	}
	fprintf(ofile, "\n");
    }
    gzip_close (ofile, sparsename);
    return W;
}

static unsigned long
flushSparse_binary(const char *sparsename, int **sparsemat, int small_nrows, int small_ncols, int *code)
{
    FILE *ofile;
    unsigned long W = 0;
    char * zip = has_suffix(sparsename, ".gz") ? ".gz" : NULL;
    uint32_t x;
    uint32_t * weights = malloc(small_ncols * sizeof(uint32_t));
    memset(weights, 0, small_ncols * sizeof(uint32_t));

    char * zname = strdup(sparsename);
    if (zip) { zname[strlen(zname)-3]='\0'; }
    char * bname = strdup(zname);
    if (has_suffix(bname, ".bin")) { bname[strlen(bname)-4]='\0'; }

    char * name = derived_filename(bname, "bin", zip);
    ofile = gzip_open(name, "w");
    for(int i = 0; i < small_nrows; i++){
	if(sparsemat[i] == NULL) {
	    x = 0;
            fwrite(&x, sizeof(uint32_t), 1, ofile);
        } else {
	    W += sparsemat[i][0];
	    x = sparsemat[i][0];
            fwrite(&x, sizeof(uint32_t), 1, ofile);
	    for(int j = 1; j <= sparsemat[i][0]; j++){
#if DEBUG >= 1
		ASSERT(code[sparsemat[i][j]] > 0);
#endif
		x=code[sparsemat[i][j]]-1;
                fwrite(&x, sizeof(uint32_t), 1, ofile);
                weights[x]++;
	    }
	}
    }
    gzip_close (ofile, name);
    free(name);

    name = derived_filename(bname, "rw.bin", zip);
    ofile = gzip_open(name, "w");
    for(int i = 0; i < small_nrows; i++){
	x = sparsemat[i] == NULL ? 0 : sparsemat[i][0];
        fwrite(&x, sizeof(uint32_t), 1, ofile);
    }
    gzip_close (ofile, name);
    free(name);

    name = derived_filename(bname, "cw.bin", zip);
    ofile = gzip_open(name, "w");
    for(int j = 0; j < small_ncols; j++){
	x = weights[j];
        fwrite(&x, sizeof(uint32_t), 1, ofile);
    }
    gzip_close (ofile, name);
    free(name);

    free(bname);
    free(zname);

    free(weights);

    return W;
}

static unsigned long
flushSparse(const char *sparsename, int **sparsemat, int small_nrows, int small_ncols, int *code, int bin)
{
    if (bin)
        return flushSparse_binary(sparsename, sparsemat, small_nrows, small_ncols, code);
    else
        return flushSparse_ascii(sparsename, sparsemat, small_nrows, small_ncols, code);
}

/*
static void
infos4Manu(char *name, char *sparsename, int small_nrows, int small_ncols, int nvslices, int *tabnc, unsigned long *Wslice)
{
    FILE *ofile;
    int j;
    int jtotal;

    sprintf(name, "%s.info", sparsename);
    fprintf(stderr, "Creating file %s\n", name);
    ofile = fopen(name, "w");
    fprintf(ofile, "%d %d\n", small_nrows, small_ncols);
    fprintf(ofile, "%d %d\n", 1, nvslices);
    const char * sparse_basename;
    sparse_basename = strrchr(sparsename, '/');
    if (sparse_basename == NULL) {
        sparse_basename = sparsename;
    } else {
        sparse_basename++;
    }

    jtotal=0;
    for(j = 0; j < nvslices; j++){
	sprintf(name, "%s.%02d", sparse_basename, j);
	fprintf(ofile, "%d %d %d %d", 0, j, 0, jtotal);
	fprintf(ofile, " %d %d %lu %s\n",small_nrows,tabnc[j],Wslice[j],name);
        jtotal += tabnc[j];
    }
    fclose(ofile);
}
*/

// dump of newrows in indexname.
static void
makeIndexFile(const char *indexname, int nrows, int **newrows, int small_nrows, int small_ncols)
{
    FILE *indexfile;
    int i, j;

    indexfile = gzip_open(indexname, "w");
    fprintf(indexfile, "%d %d\n", small_nrows, small_ncols);
    for(i = 0; i < nrows; i++)
	if(newrows[i] != NULL){
	    fprintf(indexfile, "%d", newrows[i][0]);
	    for(j = 1; j <= newrows[i][0]; j++){
		fprintf(indexfile, " ");
		fprintf(indexfile, PURGE_INT_FORMAT, newrows[i][j]);
	    }
	    fprintf(indexfile, "\n");
	}
    gzip_close(indexfile, indexname);
}

// on input, colweight[j] contains the weight; on exit, colweight[j]
// contains the new index for j.
static void
renumber(int *small_ncols, int *colweight, int ncols)
{
    int j, k, nb, *tmp;

    tmp = (int *)malloc((ncols<<1) * sizeof(int));
    memset(tmp, 0, (ncols<<1) * sizeof(int));
    for(j = 0, nb = 0; j < ncols; j++)
	if(colweight[j] > 0){
#if DEBUG >= 1
	    fprintf(stderr, "J %d %d\n", j, colweight[j]);
#endif
	    tmp[nb++] = colweight[j];
	    tmp[nb++] = j;
	}
    *small_ncols = nb>>1;
    qsort(tmp, nb>>1, 2*sizeof(int), cmp);
    memset(colweight, 0, ncols * sizeof(int));
#if 0
    // useful for Gauss only...
    fprintf(stderr, "Sorting in INcreasing weight order of j\n");
    for(j = 0, k = 1; j < nb; j += 2)
	colweight[tmp[j+1]] = k++; // always this +1 trick
#else
    // useful for BW + skipping heavy part only...
    fprintf (stderr, "Sorting columns by decreasing weight\n");
    for(j = nb-1, k = 1; j >= 0; j -= 2)
	colweight[tmp[j]] = k++; // always this +1 trick
#endif
    free(tmp);
}

// A line is "i i1 ... ik".
// If i >= 0 then
//     row[i] is to be added to rows i1...ik and destroyed at the end of
//     the process.
//     Works also is i is alone (hence: destroyed row).
// If i < 0 then
//     row[-i-1] is to be added to rows i1...ik and NOT destroyed.
//
static void
doAllAdds(int **newrows, char *str)
{
    char *t = str;
    int i, ii, destroy = 1;

    if(*t == '-'){
	destroy = 0;
	t++;
    }
    for(i = 0; (*t != ' ') && (*t != '\n'); t++)
	i = 10 * i + (*t - '0');
    if(!destroy)
	i--; // what a trick, man!
#if DEBUG >= 1
    fprintf(stderr, "first i is %d\n", i);
#endif
    if(*t != '\n'){
	++t;
	ii = 0;
	while(1){
	    if((*t == '\n') || (*t == ' ')){
#if DEBUG >= 1
		fprintf(stderr, "next ii is %d\n", ii);
#endif
		addRows(newrows, ii, i, -1);
		ii = 0;
	    }
	    else
		ii = 10 * ii + (*t - '0');
	    if(*t == '\n')
		break;
	    t++;
	}
    }
    if(destroy){
	// destroy initial row!
	free(newrows[i]);
	newrows[i] = NULL;
    }
}

#define STRLENMAX 2048

static int
oneFile(const char *sparsename, int **sparsemat, int *colweight, const char *purgedname, FILE *purgedfile, int **oldrows, int ncols, int small_nrows, int verbose, int bin)
{
    unsigned long W;
    int small_ncols;

    makeSparse(sparsemat, colweight, purgedfile, 0, ncols, oldrows, verbose,0);
    gzip_close(purgedfile, purgedname);

    fprintf(stderr, "Renumbering columns (including sorting w.r.t. weight)\n");
    renumber(&small_ncols, colweight, ncols);

    fprintf(stderr, "small_nrows=%d small_ncols=%d\n",small_nrows,small_ncols);

    double tt = seconds();
    fprintf(stderr, "Writing sparse representation to file\n");
    W = flushSparse(sparsename, sparsemat, small_nrows, small_ncols, colweight, bin);
    fprintf(stderr, "#T# writing sparse: %2.2lf\n", seconds()-tt);
    fprintf(stderr, "# Weight(M_small) = %lu\n", W);

    return small_ncols;
}

static int
manyFiles(const char *sparsename, int **sparsemat, int *colweight, const char *purgedname, FILE *purgedfile, int **oldrows, int ncols, int small_nrows, int verbose, int bin, int nslices)
{
    char *name;
    int small_ncols = 1; // always the +1 trick
    int slice, i, jmin, jmax, jstep, j, *tabnc;
    unsigned long *Wslice;

    // to add ".00" or to add ".infos"
    name = (char *)malloc((strlen(sparsename)+7) * sizeof(char));
    jmin = jmax = 0;
    // tabnc[slice] = number of columns in slice
    tabnc = (int *)malloc(nslices * sizeof(int));
    Wslice = (unsigned long *)malloc(nslices * sizeof(unsigned long));
    // we want nslices always
    if((ncols % nslices) == 0)
	jstep = ncols/nslices;
    else
	jstep = ncols/(nslices-1);
    for(slice = 0; slice < nslices; slice++){
	// we operate on [jmin..jmax[
	jmin = jmax;
	jmax += jstep;
	if(slice == (nslices-1))
	    jmax = ncols;
	sprintf(name, "%s.%02d", sparsename, slice);
	fprintf(stderr, "Dealing with M[%d..%d[ -> %s\n", jmin, jmax, name);
	makeSparse(sparsemat,colweight,purgedfile,jmin,jmax,oldrows,verbose,
		   nslices);
	gzip_close(purgedfile, purgedname);
	// colweight[jmin..jmax[ contains the new weights of the
	// corresponding columns
	tabnc[slice] = 0;
	for(j = jmin; j < jmax; j++)
	    if(colweight[j] > 0){
		colweight[j] = small_ncols++;
		tabnc[slice]++;
	    }
	// warning: take care to the +1 trick
	Wslice[slice] = flushSparse(name, sparsemat, small_nrows, small_ncols,
				    colweight, bin);
	// do not forget to clean the rows in sparsemat to be ready for
	// next time if any
	if(slice < nslices-1){
	    for(i = 0; i < small_nrows; i++)
		if(sparsemat[i] != NULL){
		    free(sparsemat[i]);
		    sparsemat[i] = NULL;
		}
	    purgedfile = gzip_open(purgedname, "r");
	    ASSERT(purgedfile != NULL);
	}
    }
    small_ncols--; // undoing the +1 trick
    // infos4Manu(name,sparsename,small_nrows,small_ncols,nslices,tabnc,Wslice);
    free(name);
    free(tabnc);
    free(Wslice);
    return small_ncols;
}

static void
build_newrows_from_file(int **newrows, int nrows, FILE *hisfile, uint64_t bwcostmin)
{
    uint64_t bwcost;
    unsigned long addread = 0;
    int i;
    char str[STRLENMAX];

    for(i = 0; i < nrows; i++){
	newrows[i] = (int *)malloc(2 * sizeof(int));
	newrows[i][0] = 1;
	newrows[i][1] = i;
    }
    fprintf(stderr, "Reading row additions\n");
    while(fgets(str, STRLENMAX, hisfile)){
	addread++;
	if((addread % 100000) == 0)
	    fprintf(stderr, "%lu lines read at %2.2lf\n", addread, seconds());
	if(str[strlen(str)-1] != '\n'){
	    fprintf(stderr, "Gasp: not a complete a line!");
	    fprintf(stderr, " I stop reading and go to the next phase\n");
	    break;
	}
	if(strncmp(str, "BWCOST", 6) != 0)
	    doAllAdds(newrows, str);
	else{
	    if(strncmp(str, "BWCOSTMIN", 9) != 0){
		sscanf(str+8, "%"PRIu64"", &bwcost);
		//	fprintf(stderr, "Read bwcost=%"PRIu64"\n", bwcost);
		if((bwcostmin != 0) && (bwcost == bwcostmin)){
		    // what a damn tricky feature!!!!!!!
		    fprintf(stderr, "Activating tricky stopping feature");
		    fprintf(stderr, " since I reached %"PRIu64"\n", bwcostmin);
		    break;
		}
	    }
	}
    }
}

static void
read_newrows_from_file(int **newrows, int nrows, FILE *file)
{
    int small_nrows, small_ncols, i, nc, k, *tmp;
    int rc;

    rc = fscanf(file, "%d %d", &small_nrows, &small_ncols);
    ASSERT_ALWAYS(rc == 2);

    for(i = 0; i < small_nrows; i++){
	rc = fscanf(file, "%d", &nc);
        ASSERT_ALWAYS(rc == 1);

	tmp = (int *)malloc((1+nc) * sizeof(int));
	tmp[0] = nc;
	for(k = 0; k < nc; k++) {
	    rc = fscanf(file, PURGE_INT_FORMAT, tmp+k+1);
            ASSERT_ALWAYS(rc == 1);
        }

	newrows[i] = tmp;
    }
    for( ; i < nrows; i++)
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
    FILE *hisfile, *purgedfile, *fromfile;
    uint64_t bwcostmin = 0;
    int nrows, ncols, nslices = 0;
    int **newrows, i, j, nb, *nbrels, **oldrows, *colweight;
    int ind, small_nrows, small_ncols, **sparsemat;
    int verbose = 0, writeindex;
    char str[STRLENMAX];
    int bin=0;
    char * rp;

    // printing the arguments as everybody does these days
    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

    param_list(pl);

    param_list_init(pl);
    argv++,argc--;
    param_list_configure_knob(pl, "--verbose", &verbose);
    param_list_configure_knob(pl, "--binary", &bin);
    param_list_configure_alias(pl, "--verbose", "-v");

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        abort();
    }
    const char * purgedname = param_list_lookup_string(pl, "purged");
    const char * hisname = param_list_lookup_string(pl, "his");
    const char * sparsename = param_list_lookup_string(pl, "out");
    const char * indexname = param_list_lookup_string(pl, "index");
    const char * fromname = param_list_lookup_string(pl, "from");
    param_list_parse_int(pl, "binary", &bin);
    param_list_parse_int(pl, "nslices", &nslices);
    param_list_parse_uint64(pl, "bwcostmin", &bwcostmin);
    if (has_suffix(sparsename, ".bin") || has_suffix(sparsename, ".bin.gz"))
        bin=1;

    purgedfile = gzip_open(purgedname, "r");
    ASSERT(purgedfile != NULL);

    hisfile = fopen(hisname, "r");
    ASSERT(hisfile != NULL);
    rp = fgets(str, STRLENMAX, hisfile);
    ASSERT_ALWAYS(rp);

    // read parameters that should be the same as in purgedfile!
    sscanf(str, "%d %d", &nrows, &ncols);
    fprintf(stderr, "Original matrix has size %d x %d\n", nrows, ncols);
    newrows = (int **)malloc(nrows * sizeof(int *));

    if(fromname == NULL){
	writeindex = 1;
	build_newrows_from_file(newrows, nrows, hisfile, bwcostmin);
    }
    else{
	writeindex = 0;
	fromfile = fopen(fromname, "r");
	ASSERT(fromfile != NULL);
	read_newrows_from_file(newrows, nrows, fromfile);
    }
    fclose(hisfile);

    nbrels = (int *) malloc(nrows * sizeof(int));
    memset (nbrels, 0, nrows * sizeof(int));
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
    free (nbrels);
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
    if(nslices == 0) {
	small_ncols = oneFile(sparsename, sparsemat, colweight, 
			      purgedname, purgedfile,
			      oldrows, ncols, small_nrows, verbose, bin);
    } else {
	small_ncols = manyFiles(sparsename, sparsemat, colweight, 
				purgedname, purgedfile,
				oldrows, ncols, small_nrows, verbose, bin,
				nslices);
	exit(0);
    }
    if(writeindex){
	// this part depends on newrows only
	double tt = seconds();
	fprintf(stderr, "Writing index file\n");
	makeIndexFile(indexname, nrows, newrows, small_nrows, small_ncols);
	fprintf(stderr, "#T# writing index file: %2.2lf\n", seconds()-tt);
    }

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
    param_list_clear(pl);
    free(newrows);
    free(oldrows);
    free(colweight);
    return 0;
}
