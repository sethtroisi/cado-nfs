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
renumber (int *small_ncols, int *colweight, int ncols)
{
    int j, k, nb, *tmp;

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
    memset(colweight, 0, ncols * sizeof(int));
    // useful for BW + skipping heavy part only...
    for(j = nb-1, k = 1; j >= 0; j -= 2)
	colweight[tmp[j]] = k++; // always this +1 trick
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
toFlush (const char *sparsename, int **sparsemat, int *colweight, int ncols,
         int small_nrows, int skip, int bin)
{
    unsigned long W;
    int small_ncols;

    fprintf(stderr, "Renumbering columns (including sorting w.r.t. weight)\n");
    renumber (&small_ncols, colweight, ncols);

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

#if DEBUG >= 1
/* weight of merge of row i and row k */
static uint32_t
weight_merge (int **newrows, int i, int k, int skip)
{
  uint32_t w, li, lk, ix, kx;

  ASSERT(newrows[i] != NULL);
  ASSERT(newrows[k] != NULL);
  li = newrows[i][0];
  lk = newrows[k][0];
  for (w = 0, ix = 1, kx = 1; ix <= li && kx <= lk; )
    {
      if (newrows[i][ix] < newrows[k][kx])
        w += newrows[i][ix++] >= skip;
      else if (newrows[i][ix] == newrows[k][kx])
        ix++, kx++;
      else
        w += newrows[k][kx++] >= skip;
    }
  while (ix <= li)
    w += newrows[i][ix++] >= skip;
  while (kx <= lk)
    w += newrows[k][kx++] >= skip;
  return w;
}
#endif

void
add_relation (int **cols, int *len_col, int j, int i)
{
  len_col[j] ++;
  cols[j] = (int*) realloc (cols[j], len_col[j] * sizeof(int));
  cols[j][len_col[j]-1] = i;
}

void
sub_relation (int **cols, int *len_col, int j, int i)
{
  int k;

  for (k = 0; k < len_col[j]; k++)
    if (cols[j][k] == i)
      break;
  if (k >= len_col[j])
    {
      fprintf (stderr, "Error, row i=%d not found in column j=%d\n", i, j);
      exit (1);
    }
  len_col[j] --;
  cols[j][k] = cols[j][len_col[j]];
  cols[j] = (int*) realloc (cols[j], len_col[j] * sizeof(int));
}

#if DEBUG >= 1
static void
print_row (int **newrows, int i)
{
  int j;

  printf ("%d:", newrows[i][0]);
  for (j = 1; j <= newrows[i][0]; j++)
    printf (" %d", newrows[i][j]);
  printf ("\n");
}

static void
print_col (int **cols, int *len_col, int j)
{
  int i;

  for (i = 0; i < len_col[j]; i++)
    printf ("%d ", cols[j][i]);
  printf ("\n");
}

static void
check_row (int **newrows, int i, int **cols, int *len_col)
{
  int j, k, l;

  for (k = 1; k <= newrows[i][0]; k++)
    {
      /* check that newrows[i][k] is in the tranpose matrix */
      j = newrows[i][k];
      for (l = 0; l < len_col[j]; l++)
        if (cols[j][l] == i)
          break;
      if (l >= len_col[j])
        {
          fprintf (stderr, "Error, row %d contains column %d but column %d does not contain row %d\n", i, j, j, i);
          exit (1);
        }
    }
}

static void
check_col (int **cols, int *len_col, int j, int **newrows)
{
  int i, k, l;

  for (l = 0; l < len_col[j]; l++)
    {
      i = cols[j][l];
      for (k = 1; k <= newrows[i][0]; k++)
        if (newrows[i][k] == j)
          break;
      ASSERT_ALWAYS (k <= newrows[i][0]);
    }
}
#endif

static int
weight_row (int **newrows, int i, int skip)
{
  int w = 0, k;

  ASSERT(newrows[i] != NULL);
  for (k = 1; k <= newrows[i][0]; k++)
    w += newrows[i][k] >= skip;
  return w;
}

/* row cols[j][i] += cols[j][k] */
static void
do_merge (int **newrows, int **cols, int *len_col, int j, int i, int k,
          int skip, FILE *hisfile)
{
  int ii, kk, li, lk, *tmp, ni, nk, ltmp;

  ii = cols[j][i];
  kk = cols[j][k];
  fprintf (hisfile, "-%u %u\n", kk + 1, ii);
  fflush (hisfile);
  ASSERT(newrows[ii] != NULL);
  ASSERT(newrows[kk] != NULL);
  li = newrows[ii][0];
  lk = newrows[kk][0];
  ASSERT(weight_row (newrows, ii, skip) >= weight_row (newrows, kk, skip));
  tmp = (int*) malloc ((1 + li + lk) * sizeof(int));
  for (ni = 1, nk = 1, ltmp = 0; ni <= li && nk <= lk;)
    {
      if (newrows[ii][ni] < newrows[kk][nk])
        {
          /* newrows[ii][ni] is already in row ii */
          tmp[++ltmp] = newrows[ii][ni++];
        }
      else if (newrows[ii][ni] > newrows[kk][nk])
        {
          /* newrows[kk][nk] is new in row ii */
          if (newrows[kk][nk] >=  skip)
            add_relation (cols, len_col, newrows[kk][nk], ii);
          tmp[++ltmp] = newrows[kk][nk++];
        }
      else
        {
          /* newrows[ii][ni] disappears in row ii */
          if (newrows[ii][ni] >= skip)
            sub_relation (cols, len_col, newrows[ii][ni], ii);
          ni++, nk++;
        }
    }
  /* only one of the following loops is non-empty */
  while (ni <= li)
    tmp[++ltmp] = newrows[ii][ni++];
  while (nk <= lk)
    {
      if (newrows[kk][nk] >= skip)
        add_relation (cols, len_col, newrows[kk][nk], ii);
      tmp[++ltmp] = newrows[kk][nk++];
    }
  free (newrows[ii]);
  tmp = (int*) realloc (tmp, (1 + ltmp) * sizeof(int));
  tmp[0] = ltmp;
  newrows[ii] = tmp;
}

static int
try_merge (int **newrows, int **cols, int *len_col, int j, int skip,
           FILE *hisfile)
{
  int i, k, gain, gain_max, imax, kmax, *W, ii, kk, *J, s, t;
  mpz_t *M, and;

  /* compute weights of rows containing the ideal j */
  W = (int*) malloc (len_col[j] * sizeof(int));
  for (i = 0, s = 0; i < len_col[j]; i++)
    {
      W[i] = weight_row (newrows, cols[j][i], skip);
      s += W[i];
    }

  /* compute all ideals appearing in the rows containing j */
  if (j >= skip)
    s -= len_col[j];
  J = (int*) malloc (s * sizeof(int));
  for (ii = 0, t = 0; ii < len_col[j]; ii++)
    {
      i = cols[j][ii];
      for (k = 1; k <= newrows[i][0]; k++)
        if ((newrows[i][k] >= skip) && (newrows[i][k] != j))
          J[t++] = newrows[i][k];
    }
  ASSERT(t == s);
  qsort (J, t, sizeof(int), cmp);

  /* extract ideals appearing at least twice */
  for (i = 1, k = 0; i < t; i++)
    if ((J[i-1] == J[i]) && (k == 0 || (J[k-1] != J[i])))
      J[k++] = J[i];
  J[k++] = INT_MAX;
  t = k;

  /* compute bit vectors for ideals appearing at least twice */
  M = (mpz_t*) malloc (len_col[j] * sizeof(mpz_t));
  mpz_init (and);
  for (ii = 0; ii < len_col[j]; ii++)
    {
      mpz_init (M[ii]);
      mpz_realloc2 (M[ii], t - 1);
      i = cols[j][ii];
      for (s = 0, k = 1; k <= newrows[i][0]; k++)
        {
          while (J[s] < newrows[i][k])
            s++;
          /* now J[s] >= newrows[i][k] */
          if (J[s] == newrows[i][k])
            mpz_setbit (M[ii], s++);
        }
    }
  for (gain_max = 0, i = 0; i < len_col[j]; i++)
    {
      for (k = i + 1; k < len_col[j]; k++)
        {
          if (W[i] >= W[k])
            ii = i, kk = k;
          else
            ii = k, kk = i;
          /* we need wmin > gain_max since the maximum gain is wmin */
          if (W[kk] <= gain_max)
            continue;
          mpz_and (and, M[ii], M[kk]);
          /* the +2 accounts for the ideal j */
          gain = 2 * mpz_popcount (and) + 2 - W[kk];
#if DEBUG >= 1
          {
            int wik;
            wik = weight_merge (newrows, cols[j][ii], cols[j][kk], skip);
            ASSERT_ALWAYS (gain == W[ii] - wik);
          }
#endif
          if (gain > gain_max)
            {
              gain_max = gain;
              imax = ii;
              kmax = kk;
            }
        }
    }
  mpz_clear (and);
  for (i = 0; i < len_col[j]; i++)
    mpz_clear (M[i]);
  free (M);
  free (J);
  free (W);
  if (gain_max > 0)
    do_merge (newrows, cols, len_col, j, imax, kmax, skip, hisfile);
  return gain_max;
}

static int
cmp_ge (const void *p, const void *q)
{
  int x = *((int *)p);
  int y = *((int *)q);
  return (x >= y ? -1 : 1);
}

/* we append the new merges in file 'hisname', so that they are considered
   when writing the index file afterwards */
static void
optimize (int **newrows, int nrows, int *colweight, int ncols, int skip,
          const char *hisname)
{
  int **cols, *len_col, i, j, k, small_ncols, pass = 0, *perm_cols,
    small_nrows;
  uint64_t total_weight;
  FILE *hisfile;

  hisfile = fopen (hisname, "a");
  ASSERT_ALWAYS(hisfile != NULL);

  /* first permute columns to get heavier columns first */
  perm_cols = (int*) malloc (2 * ncols * sizeof(int));
  for (j = 0; j < ncols; j++)
    {
      perm_cols[2*j] = colweight[j];
      perm_cols[2*j+1] = j;
    }
  qsort (perm_cols, ncols, 2*sizeof(int), cmp_ge);
#if DEBUG >= 1
  for (j = 1; j < ncols; j++)
    ASSERT_ALWAYS(perm_cols[2*(j-1)] >= perm_cols[2*j]);
#endif
  printf ("Skip %d heaviest columns\n", skip);
  for (j = 0; j < ncols && perm_cols[2*j] != 0; j++);
  small_ncols = j;
  printf ("New number of columns: %d\n", small_ncols);

  /* compute the inverse permutation in colweight[] */
  for (j = 0; j < ncols; j++)
    colweight[perm_cols[2*j+1]] = j;

  ncols = j;

  /* renumber the direct matrix */
  for (i = 0, small_nrows = 0; i < nrows; i++)
    if (newrows[i] != NULL)
      {
        small_nrows ++;
        for (k = 1; k <= newrows[i][0]; k++)
          {
            j = newrows[i][k];
            newrows[i][k] = colweight[j];
          }
        qsort (newrows[i] + 1, newrows[i][0], sizeof(int), cmp);
#if DEBUG >= 1
        for (k = 2; k <= newrows[i][0]; k++)
          ASSERT_ALWAYS(newrows[i][k-1] <= newrows[i][k]);
#endif
      }
  /* update colweight */
  for (j = 0; j < ncols; j++)
    colweight[j] = perm_cols[2*j];
  free (perm_cols);

  cols = (int**) malloc (ncols * sizeof(int*));
  len_col = (int*) malloc (ncols * sizeof(int));
  for (j = 0, total_weight = 0; j < ncols; j++)
    {
      if (j < skip)
        cols[j] = NULL;
      else
        {
          cols[j] = (int*) malloc (colweight[j] * sizeof(int));
          total_weight += colweight[j];
        }
      len_col[j] = 0;
    }

  printf ("Optimize: small_nrows=%d small_ncols=%d weight=%lu (av. %1.2f)\n",
          small_nrows, small_ncols, total_weight,
          (double) total_weight / (double) small_nrows);

  /* compute transpose matrix */
  for (i = 0; i < nrows; i++)
    if (newrows[i] != NULL)
      {
        for (k = 1; k <= newrows[i][0]; k++)
          {
            j = newrows[i][k];
            /* we assume the row elements are sorted by increasing order */
            if (k > 1 && newrows[i][k-1] > j)
              printf ("i=%d k=%d newrows[i][k-1]=%d j=%d\n", i, k,
                      newrows[i][k-1], j);
            ASSERT(k == 1 || newrows[i][k-1] <= j);
            if (j >= skip)
              {
                cols[j][len_col[j]] = i;
                len_col[j]++;
              }
          }
      }

  printf ("Computed transpose matrix\n");
  fflush (stdout);

#if DEBUG >= 1
  for (j = skip; j < ncols; j++)
    ASSERT_ALWAYS(len_col[j] == colweight[j]);
#endif

  bit_vector active;
  bit_vector_init_set (active, ncols, 1);
  do
    {
      int gain;

      printf ("Pass %d: decrease total weight to ", ++pass);
      fflush (stdout);
      for (j = skip; j < ncols; j++)
        if (len_col[j] <= pass + 1 && bit_vector_getbit (active, j))
          {
            gain = try_merge (newrows, cols, len_col, j, skip, hisfile);
            if (gain == 0) /* we assume a column which does not give any
                              gain will never give a gain in the future */
              bit_vector_clearbit (active, j);
            total_weight -= gain;
          }
      printf ("%lu (%1.2f per row)\n", total_weight,
              (double) total_weight / (double) small_nrows);
      fflush (stdout);
    }
  while (pass < 20);
  bit_vector_clear (active);

#if DEBUG >= 1
  for (j = 0; j < skip; j++)
    ASSERT_ALWAYS(len_col[j] == 0);
#endif

  /* recompute colweight[] since it is wrong for heavy columns */
  memset (colweight, 0, skip * sizeof(int));
  for (i = 0; i < nrows; i++)
    if (newrows[i] != NULL)
      for (k = 1; k <= newrows[i][0]; k++)
        {
          j = newrows[i][k];
          if (j >= skip)
            break;
          colweight[j] ++;
        }

  for (j = 0; j < ncols; j++)
    free (cols[j]);
  free (len_col);
  free (cols);
  fclose (hisfile);
}

static void
fasterVersion(int **newrows, const char *sparsename, const char *indexname,
	      const char *hisname, purgedfile_stream ps, uint64_t bwcostmin,
              int nrows, int ncols, int skip, int bin, int writeindex)
{
    FILE *hisfile;
    int *colweight;
    int small_nrows, small_ncols;
    char str[STRLENMAX], *rp MAYBE_UNUSED;

    hisfile = fopen (hisname, "r");
    ASSERT_ALWAYS(hisfile != NULL);
    /* read first line */
    rp = fgets (str, STRLENMAX, hisfile);

    // 1st pass: read the relations in *.purged and put them in newrows[][]
    readPurged(newrows, ps, 1);
#if DEBUG >= 1
    for(int i = 0; i < nrows; i++){
	printf("row[%d]=", i);
	for(int k = 1; k <= newrows[i][0]; k++)
	    printf(" %d", newrows[i][k]);
	printf("\n");
    }
#endif

    // read merges in the *.merge.his file and replay them
    build_newrows_from_file(newrows, hisfile, bwcostmin);

    /* compute column weights */
    colweight = (int*) malloc (ncols * sizeof(int *));
    ASSERT_ALWAYS(colweight != NULL);
    memset (colweight, 0, ncols * sizeof(int *));
    for (int i = 0; i < nrows; i++)
      if (newrows[i] != NULL)
        for(int k = 1; k <= newrows[i][0]; k++)
          colweight[newrows[i][k]] += 1;

    /* comment out the following line to disable cycle optimization */
    optimize (newrows, nrows, colweight, ncols, skip, hisname);

    /* crunch empty rows */
    for (int i = small_nrows = 0; i < nrows; i++)
      if (newrows[i] != NULL)
        newrows[small_nrows++] = newrows[i];

    /* renumber columns after sorting them by decreasing weight */
    small_ncols = toFlush(sparsename, newrows, colweight, ncols,
			  small_nrows, skip, bin);
    free (colweight);
    if(writeindex)
        toIndex(newrows, indexname, hisfile, bwcostmin, nrows,
                small_nrows, small_ncols);
    for(int i = 0; i < nrows; i++)
       if(newrows[i] != NULL)
        free (newrows[i]);
    free(newrows);

    fclose (hisfile);
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
    FILE *hisfile;
    uint64_t bwcostmin = 0;
    int nrows, ncols, nslices = 0;
    int **newrows;
    int verbose = 0;
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
    param_list_parse_int(pl, "binary", &bin);
    param_list_parse_int(pl, "nslices", &nslices);
    param_list_parse_int(pl, "skip", &skip);
    param_list_parse_uint64(pl, "bwcostmin", &bwcostmin);
    if (has_suffix(sparsename, ".bin") || has_suffix(sparsename, ".bin.gz"))
        bin=1;

    purgedfile_stream ps;
    purgedfile_stream_init(ps);
    purgedfile_stream_openfile(ps, purgedname);

    hisfile = fopen (hisname, "r");
    ASSERT_ALWAYS(hisfile != NULL);
    rp = fgets(str, STRLENMAX, hisfile);
    ASSERT_ALWAYS(rp);
    fclose (hisfile);

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

    fasterVersion (newrows, sparsename, indexname, hisname, ps,
                   bwcostmin, nrows, ncols, skip, bin, 1);


    purgedfile_stream_closefile(ps);
    purgedfile_stream_clear(ps);

    param_list_clear(pl);

    print_timing_and_memory (wct0);

    return 0;
}
