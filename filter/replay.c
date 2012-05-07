/* replay --- replaying history of merges to build the small matrix

Copyright 2008, 2009, 2010, 2011, 2012 Francois Morain, Emmanuel Thome,
                                       Paul Zimmermann

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

#include "sparse.h"
#include "gzip.h"

#define DEBUG 0

#define TRACE_COL -1 // 231 // put to -1 if not...!

#if DEBUG >= 1
static unsigned long row_additions = 0;
#endif

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

static unsigned long
flushSparse(const char *sparsename, int **sparsemat, int small_nrows, int small_ncols, int *code, int skip, int bin)
{
    const struct {
        const char * ext;
        const char * smat;
        const char * srw;
        const char * scw;
        const char * dmat;
        const char * drw;
        const char * dcw;
    } suffixes[2] = {
        {
          .ext = ".txt",
          .smat = "txt",
          .srw = "rw.txt",
          .scw = "cw.txt",
          .dmat = "dense.txt",
          .drw = "dense.rw.txt",
          .dcw = "dense.cw.txt",
        },
        {
          .ext = ".bin",
          .smat = "bin",
          .srw = "rw.bin",
          .scw = "cw.bin",
          .dmat = "dense.bin",
          .drw = "dense.rw.bin",
          .dcw = "dense.cw.bin",
        },
    }, * suf = &(suffixes[bin]);

    unsigned long W = 0;
    unsigned long DW = 0;
    char * zip = has_suffix(sparsename, ".gz") ? ".gz" : NULL;
    uint32_t * weights = malloc(small_ncols * sizeof(uint32_t));
    ASSERT_ALWAYS(weights != NULL);
    memset(weights, 0, small_ncols * sizeof(uint32_t));

    char * base = strdup(sparsename);
    if (zip) { base[strlen(base)-3]='\0'; }
    if (has_suffix(base, suf->ext)) { base[strlen(base)-4]='\0'; }

    char * smatname = NULL;
    FILE * smatfile = NULL;
    char * srwname  = NULL;
    FILE * srwfile  = NULL;
    char * scwname  = NULL;
    FILE * scwfile  = NULL;

    {
        smatname = derived_filename(base, suf->smat, zip);
        srwname  = derived_filename(base, suf->srw, zip);
        smatfile = gzip_open(smatname, "w");
        srwfile  = gzip_open(srwname, "w");
        if (!bin)
            fprintf(smatfile, "%d %d\n", small_nrows, small_ncols - skip);
    }

    char * dmatname = NULL;
    FILE * dmatfile = NULL;
    char * drwname  = NULL;
    FILE * drwfile  = NULL;
    char * dcwname  = NULL;
    FILE * dcwfile  = NULL;

    if (skip) {
        dmatname = derived_filename(base, suf->dmat, zip);
        drwname  = derived_filename(base, suf->drw, zip);
        dmatfile = gzip_open(dmatname, "w");
        drwfile  = gzip_open(drwname, "w");
        if (!bin)
            fprintf(dmatfile, "%d %d\n", small_nrows, skip);
    }

    for(int i = 0; i < small_nrows; i++){
	if(sparsemat[i] == NULL) {
          /* this should not happen, unless the corresponding combination of
             relations yields a square, thus we have found a dependency in the
             merge process */
            if (bin) {
                const uint32_t x = 0;
                fwrite32_little(&x, 1, smatfile);
                if (skip) fwrite32_little(&x, 1, dmatfile);
            } else {
                fprintf(smatfile, "0");
                if (skip) fprintf(dmatfile, "0");
            }
        } else {
            uint32_t dw = 0;
            uint32_t sw = 0;
	    for(int j = 1; j <= sparsemat[i][0]; j++){
		if (code[sparsemat[i][j]]-1 < skip) {
                    dw++;
                    DW++;
                } else {
                    sw++;
                    W++;
                }
            }
            if (bin) {
                fwrite32_little(&sw, 1, smatfile);
                fwrite32_little(&sw, 1, srwfile);
                if (skip) fwrite32_little(&dw, 1, dmatfile);
                if (skip) fwrite32_little(&dw, 1, drwfile);
            } else {
                fprintf(smatfile, "%"PRIu32"", sw);
                fprintf(srwfile, "%"PRIu32"\n", sw);
                if (skip) fprintf(dmatfile, "%"PRIu32"", dw);
                if (skip) fprintf(drwfile, "%"PRIu32"\n", dw);
            }
	    for(int j = 1; j <= sparsemat[i][0]; j++){
#if DEBUG >= 1
		ASSERT(code[sparsemat[i][j]] > 0);
#endif
		uint32_t x = code[sparsemat[i][j]]-1;
                weights[x]++;
                if ((int) x < skip) {
                    ASSERT_ALWAYS(skip);
                    if (bin) {
                        fwrite32_little(&x, 1, dmatfile);
                    } else {
                        fprintf(dmatfile, " %"PRIu32"", x);
                    }
                } else {
                    x-=skip;
                    if (bin) {
                        fwrite32_little(&x, 1, smatfile);
                    } else {
                        fprintf(smatfile, " %"PRIu32"", x);
                    }
                }
	    }
	}
	if (!bin) {
            fprintf(smatfile, "\n");
            if (skip) fprintf(dmatfile, "\n");
        }
    }

    {
        gzip_close (smatfile, smatname);
        gzip_close (srwfile, srwname);
        free(smatname);
        free(srwname);
    }
    if (skip) {
        fprintf(stderr, "%lu coeffs (out of %lu total) put into %s (%.1f%%)\n",
                DW, DW+W, dmatname,
                100.0 * (double) DW / (DW+W));

        gzip_close (dmatfile, dmatname);
        gzip_close (drwfile, drwname);
        free(dmatname);
        free(drwname);
    }

    if (skip) {
        dcwname = derived_filename(base, suf->dcw, zip);
        dcwfile = gzip_open(dcwname, "w");
        for(int j = 0; j < skip; j++){
            uint32_t x = weights[j];
            if (bin) {
                fwrite32_little(&x, 1, dcwfile);
            } else {
                fprintf(dcwfile, "%"PRIu32"\n", x);
            }
        }
        gzip_close(dcwfile, dcwname);
        free(dcwname);
    }

    {
        scwname = derived_filename(base, suf->scw, zip);
        scwfile = gzip_open(scwname, "w");
        for(int j = skip; j < small_ncols; j++){
            uint32_t x = weights[j];
            if (bin) {
                fwrite32_little(&x, 1, scwfile);
            } else {
                fprintf(scwfile, "%"PRIu32"\n", x);
            }
        }
        gzip_close(scwfile, scwname);
        free(scwname);
    }

    free(base);
    free(weights);

    return W;
}

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
// contains the new index for j. Heavier columns are in front of the new
// matrix.
static void
renumber(const char * sosname, int *small_ncols, int *colweight, int ncols)
{
    int j, k, nb, *tmp;
    uint64_t w;

    tmp = (int *)malloc((ncols<<1) * sizeof(int));
    ASSERT_ALWAYS(tmp != NULL);
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
    fprintf (stderr, "Sorting %d columns by decreasing weight\n",
             *small_ncols);
    qsort(tmp, nb>>1, 2*sizeof(int), cmp);
    for (w = 0, j = 2; j <= 4000 && j <= nb; j += 2)
      {
        w += (uint64_t) tmp[nb - j];
        if ((j & 63) == 0)
          fprintf (stderr, "Total weight of heaviest %d columns is %"
                   PRIu64 "\n", j >> 1, w);
      }
    memset(colweight, 0, ncols * sizeof(int));
#if 0
    // useful for Gauss only...
    fprintf(stderr, "Sorting in Increasing weight order of j\n");
    for(j = 0, k = 1; j < nb; j += 2)
	colweight[tmp[j+1]] = k++; // always this +1 trick
#else
    // useful for BW + skipping heavy part only...
    for(j = nb-1, k = 1; j >= 0; j -= 2)
	colweight[tmp[j]] = k++; // always this +1 trick
    if (sosname) {
        fprintf(stderr, "Saving ideal renumbering info to %s\n", sosname);
        FILE * f = fopen(sosname, "w");
        ASSERT_ALWAYS(f);
        for(j = nb-1, k = 0; j >= 0; j -= 2) {
            // column number k in the output matrix corresponds to ideal
            // number tmp[j] (whose meaning can be fetched from the sos
            // file produced by purge)
            uint32_t z = tmp[j];
            int r = fwrite32_little(&z, 1, f);
            ASSERT_ALWAYS(r == 1);
        }
        fclose(f);
    }

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
    fprintf(stderr, "first i is %d in %s", i, str);
#endif
    if(*t != '\n'){
	++t;
	ii = 0;
	while(1){
	    if((*t == '\n') || (*t == ' ')){
#if DEBUG >= 1
		fprintf(stderr, "next ii is %d\n", ii);
                row_additions ++;
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

// sparsemat is small_nrows x small_ncols, after small_ncols is found using
// renumbering.
static int
toFlush(const char *sparsename, const char * sosname, int **sparsemat, int *colweight, int ncols, int small_nrows, int skip, int bin)
{
    unsigned long W;
    int small_ncols;

    fprintf(stderr, "Renumbering columns (including sorting w.r.t. weight)\n");
    renumber(sosname, &small_ncols, colweight, ncols);

    fprintf(stderr, "small_nrows=%d small_ncols=%d\n",small_nrows,small_ncols);

    double tt = seconds();
    fprintf(stderr, "Writing sparse representation to file\n");
    W = flushSparse(sparsename, sparsemat, small_nrows, small_ncols, colweight, skip, bin);
    fprintf(stderr, "#T# writing sparse: %2.2lf\n", seconds()-tt);
    fprintf(stderr, "# Weight(M_small) = %lu\n", W);

    return small_ncols;
}

static void
build_newrows_from_file(int **newrows, FILE *hisfile, uint64_t bwcostmin)
{
    uint64_t bwcost;
    unsigned long addread = 0;
    char str[STRLENMAX];

    fprintf(stderr, "Reading row additions\n");
    double tt = wct_seconds();
    while(fgets(str, STRLENMAX, hisfile)){
	addread++;
	if((addread % 100000) == 0)
	    fprintf(stderr, "%lu lines read at %2.2lf\n", addread, wct_seconds() - tt);
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
	ASSERT_ALWAYS(tmp != NULL);
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

// Feed sparsemat with M_purged
static void
readPurged(int **sparsemat, purgedfile_stream ps, int verbose)
{
    fprintf(stderr, "Reading sparse matrix from purged file\n");
    for(int i = 0 ; purgedfile_stream_get(ps, NULL) >= 0 ; i++) {
	if (verbose && purgedfile_stream_disp_progress_now_p(ps))
	    fprintf(stderr, "Treating old rel #%d at %2.2lf\n",
                    ps->rrows,ps->dt);
	if(ps->nc == 0)
	    fprintf(stderr, "Hard to believe: row[%d] is NULL\n", i);
	qsort(ps->cols, ps->nc, sizeof(int), cmp);
	sparsemat[i] = (int *)malloc((ps->nc+1) * sizeof(int));
	ASSERT_ALWAYS(sparsemat[i] != NULL);
	sparsemat[i][0] = ps->nc;
        for(int k = 0; k < ps->nc; k++)
	    sparsemat[i][k+1] = ps->cols[k];
    }
}

static void
toIndex(int **newrows, const char *indexname, FILE *hisfile,
	uint64_t bwcostmin, int nrows, int small_nrows, int small_ncols)
{
    char *rp, str[STRLENMAX];
    int small_nrows2;

    // rewind
    rewind(hisfile);
    rp = fgets(str, STRLENMAX, hisfile);
    ASSERT_ALWAYS(rp);
    // reallocate
    for(int i = 0; i < nrows; i++){
      /* we only need to free the "crunched" part */
        if (i < small_nrows)
          free(newrows[i]);
	newrows[i] = (int *)malloc(2 * sizeof(int));
	ASSERT_ALWAYS(newrows[i] != NULL);
	newrows[i][0] = 1;
	newrows[i][1] = i;
    }
    // replay hisfile
    build_newrows_from_file(newrows, hisfile, bwcostmin);
    // re-determining small_nrows to check
    small_nrows2 = 0;
    for(int i = 0; i < nrows; i++)
	if(newrows[i] != NULL)
	    small_nrows2++;
    ASSERT_ALWAYS(small_nrows2 == small_nrows);
    // now, make index
    double tt = wct_seconds();
    fprintf(stderr, "Writing index file\n");
    // WARNING: small_ncols is not used, but...
    makeIndexFile(indexname, nrows, newrows, small_nrows, small_ncols);
    fprintf(stderr, "#T# writing index file: %2.2lf\n", wct_seconds()-tt);
}

static void
fasterVersion(int **newrows, 
	      const char *sparsename, const char *sosname, 
	      const char *indexname, 
	      FILE *hisfile, purgedfile_stream ps,
	      uint64_t bwcostmin, int nrows, int ncols,
	      int skip, int bin, int writeindex)
{
    int *colweight;
    int small_nrows, small_ncols;

    // 1st pass

    // read M_purged
    readPurged(newrows, ps, 1);
#if DEBUG >= 1
    for(int i = 0; i < nrows; i++){
	printf("row[%d]=", i);
	for(int k = 1; k <= newrows[i][0]; k++)
	    printf(" %d", newrows[i][k]);
	printf("\n");
    }
#endif
    // replay additions
    build_newrows_from_file(newrows, hisfile, bwcostmin);
    // compute weights of columns
    colweight = (int *)malloc(ncols * sizeof(int *));
    ASSERT_ALWAYS(colweight != NULL);
    memset(colweight, 0, ncols * sizeof(int *));

    /* compute column weights and crunch empty rows */
    for (int i = small_nrows = 0; i < nrows; i++)
      {
	if (newrows[i] != NULL)
          {
	    for(int k = 1; k <= newrows[i][0]; k++)
              colweight[newrows[i][k]] += 1;
            newrows[small_nrows++] = newrows[i];
          }
      }
    small_ncols = toFlush(sparsename, sosname, newrows, colweight, ncols,
			  small_nrows, skip, bin);
    free(colweight);
    if(writeindex)
        toIndex(newrows, indexname, hisfile, bwcostmin, nrows,
                small_nrows, small_ncols);
    for(int i = 0; i < nrows; i++)
       if(newrows[i] != NULL)
        free (newrows[i]);
    free(newrows);
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
    FILE *hisfile, *fromfile;
    uint64_t bwcostmin = 0;
    int nrows, ncols, nslices = 0;
    int **newrows;
    int verbose = 0, writeindex;
    int bin=0;
    int skip=0;
    char *rp, str[STRLENMAX];
    double wct0 = wct_seconds ();

    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    // printing the arguments as everybody does these days
    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (int i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

    param_list pl;

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
    const char * sosname = param_list_lookup_string(pl, "sos");
    param_list_parse_int(pl, "binary", &bin);
    param_list_parse_int(pl, "nslices", &nslices);
    param_list_parse_int(pl, "skip", &skip);
    param_list_parse_uint64(pl, "bwcostmin", &bwcostmin);
    if (has_suffix(sparsename, ".bin") || has_suffix(sparsename, ".bin.gz"))
        bin=1;

    purgedfile_stream ps;
    purgedfile_stream_init(ps);
    purgedfile_stream_openfile(ps, purgedname);

    hisfile = fopen(hisname, "r");
    ASSERT(hisfile != NULL);
    rp = fgets(str, STRLENMAX, hisfile);
    ASSERT_ALWAYS(rp);

    // read parameters that should be the same as in purgedfile!
    sscanf(str, "%d %d", &nrows, &ncols);
    ASSERT_ALWAYS(nrows == ps->nrows);
    ASSERT_ALWAYS(ncols == ps->ncols);

    fprintf(stderr, "Original matrix has size %d x %d\n", nrows, ncols);

    newrows = (int **)malloc(nrows * sizeof(int *));
    ASSERT_ALWAYS(newrows != NULL);

    // at the end of the following operations, newrows[i] is either
    // NULL
    // or k i_1 ... i_k which means that M_small will contain a row formed
    // of the addition of the rows of indices i_1 ... i_k in the original
    // matrix
    if(fromname == NULL){
	// generic case
	writeindex = 1;
    } else {
	// rare case, probably very very rare
	abort(); // to be clarified before use...!
	writeindex = 0;
	fromfile = fopen(fromname, "r");
	ASSERT(fromfile != NULL);
	read_newrows_from_file(newrows, nrows, fromfile);
    }

    fasterVersion(newrows, sparsename, sosname, indexname, hisfile, ps,
		  bwcostmin, nrows, ncols, skip, bin, writeindex);

    fclose(hisfile);

    purgedfile_stream_closefile(ps);
    purgedfile_stream_clear(ps);

    param_list_clear(pl);

    print_timing_and_memory (wct0);

    return 0;
}
