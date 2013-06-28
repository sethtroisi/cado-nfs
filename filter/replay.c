/* replay --- replaying history of merges to build the small matrix

Copyright 2008, 2009, 2010, 2011, 2012, 2013
          Francois Morain, Emmanuel Thome, Paul Zimmermann,
          Cyril Bouvier, Pierrick Gaudry

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

#include "portability.h"
#include "utils.h"

#include "merge_opts.h"
#include "filter_matrix.h"
#include "sparse.h"
#include "gzip.h"

#define DEBUG 0
#define STAT_FFS

// newrows[i] contains a new row formed of old rows that correspond to
// true original relations (and not multirelations).
//
// After computing newrows, we deduce for all old rows the list of newrows
// containing it.

static unsigned long
flushSparse(const char *sparsename, typerow_t **sparsemat, int small_nrows,
            int small_ncols, int *code, int skip, int bin)
{
#ifdef FOR_DL
  ASSERT_ALWAYS (skip == 0);
#endif
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

    char wmode[3] = "w";
#ifdef HAVE_MINGW
    if (bin)
      strcpy (wmode, "wb");
#endif

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
        smatfile = fopen_maybe_compressed(smatname, wmode);
        srwfile  = fopen_maybe_compressed(srwname, wmode);
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
        dmatfile = fopen_maybe_compressed(dmatname, wmode);
        drwfile  = fopen_maybe_compressed(drwname, wmode);
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
	    for(int j = 1; j <= rowLength(sparsemat, i); j++){
		if (code[rowCell(sparsemat, i, j)]-1 < skip) {
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


#if 0 //#ifdef FOR_DL //for DL sort the index
      for (int k = 1; k < rowLength(sparsemat, i); k++)
        {
          for (int l = k+1; l <= rowLength(sparsemat, i); l++)
            {
              uint32_t x = code[rowCell(sparsemat, i, k)]-1;
              uint32_t y = code[rowCell(sparsemat, i, l)]-1;
              if (x > y)
                {
                  typerow_t tmp = sparsemat[i][k];
                  sparsemat[i][k] = sparsemat[i][l];
                  sparsemat[i][l] = tmp;
                }
            }
        }
#endif

	    for(int j = 1; j <= rowLength(sparsemat, i); j++){
#if DEBUG >= 1
		ASSERT(code[rowCell(sparsemat, i, j)] > 0);
#endif
		uint32_t x = code[rowCell(sparsemat, i, j)]-1;
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
#ifdef FOR_DL
                        uint32_t e = (uint32_t) sparsemat[i][j].e;
                        fwrite32_little(&e, 1, smatfile);
#endif
                    } else {
                        fprintf(smatfile, " %"PRIu32"", x);
#ifdef FOR_DL
                        fprintf(smatfile, ":%d", sparsemat[i][j].e);
#endif
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
        fclose_maybe_compressed (smatfile, smatname);
        fclose_maybe_compressed (srwfile, srwname);
        free(smatname);
        free(srwname);
    }
    if (skip) {
        fprintf(stderr, "%lu coeffs (out of %lu total) put into %s (%.1f%%)\n",
                DW, DW+W, dmatname,
                100.0 * (double) DW / (DW+W));

        fclose_maybe_compressed (dmatfile, dmatname);
        fclose_maybe_compressed (drwfile, drwname);
        free(dmatname);
        free(drwname);
    }

    if (skip) {
        dcwname = derived_filename(base, suf->dcw, zip);
        dcwfile = fopen_maybe_compressed(dcwname, wmode);
        for(int j = 0; j < skip; j++){
            uint32_t x = weights[j];
            if (bin) {
                fwrite32_little(&x, 1, dcwfile);
            } else {
                fprintf(dcwfile, "%"PRIu32"\n", x);
            }
        }
        fclose_maybe_compressed(dcwfile, dcwname);
        free(dcwname);
    }

    {
        scwname = derived_filename(base, suf->scw, zip);
        scwfile = fopen_maybe_compressed(scwname, wmode);
        for(int j = skip; j < small_ncols; j++){
            uint32_t x = weights[j];
            if (bin) {
                fwrite32_little(&x, 1, scwfile);
            } else {
                fprintf(scwfile, "%"PRIu32"\n", x);
            }
        }
        fclose_maybe_compressed(scwfile, scwname);
        free(scwname);
    }

    free(base);
    free(weights);

    return W;
}


/* we also compare x[1] and y[1] to make the code deterministic
   since in case x[0] = y[0] qsort() may give different results on
   different machines */
static int
cmp2 (const void *p, const void *q)
{
  int *x = (int*) p;
  int *y = (int*) q;

  if (x[0] < y[0])
    return -1;
  else if (x[0] > y[0])
    return 1;
  else
    return (x[1] < y[1]) ? 1 : -1;
}

// on input, colweight[j] contains the weight; on exit, colweight[j]
// contains the new index for j. Heavier columns are in front of the new
// matrix.
static void
renumber (int *small_ncols, int *colweight, int ncols,
          MAYBE_UNUSED const char *idealsfilename)
{
    int j, k, nb, *tmp;

#ifdef FOR_DL
    FILE *renumberfile = fopen (idealsfilename, "w+");
    if (renumberfile == NULL)
      {
        fprintf (stderr, "Error while opening file to save permutation of"
                         "ideals\n");
        exit(1);
      }
#endif

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
    qsort(tmp, nb>>1, 2*sizeof(int), cmp2);
    memset(colweight, 0, ncols * sizeof(int));
    // useful for BW + skipping heavy part only...
    for(j = nb-1, k = 1; j >= 0; j -= 2)
      {
        colweight[tmp[j]] = k++; // always this +1 trick
#ifdef FOR_DL
        fprintf (renumberfile, "%d %x\n", colweight[tmp[j]]-1, tmp[j]);
#endif
      }

#ifdef FOR_DL
    fclose(renumberfile);
#endif
    free(tmp);
}

// A line is "i i1 ... ik [#j]".
// If i >= 0 then
//     row[i] is to be added to rows i1...ik and destroyed at the end of
//     the process.
//     Works also is i is alone (hence: destroyed row).
// If i < 0 then
//     row[-i-1] is to be added to rows i1...ik and NOT destroyed.
//
// If given, j is the index of the column used for pivoting (used in DL).
static void
doAllAdds(typerow_t **newrows, char *str, MAYBE_UNUSED FILE *outdelfile,
        index_data_t index_data)
{
  int32_t j;
  int32_t ind[MERGE_LEVEL_MAX], i0;
  int ni, destroy;

  ni = parse_hisfile_line (ind, str, &j);

  if (ind[0] < 0)
    {
      destroy = 0;
      i0 = -ind[0]-1;
    }
  else
    {
      destroy = 1;
      i0 = ind[0];
    }
#if DEBUG >= 1
    fprintf(stderr, "first i is %d in %s", i0, str);
#endif

  for (int k = 1; k < ni; k++)
      addRowsUpdateIndex(newrows, index_data, ind[k], i0, j);

  if(destroy)
    {
      //destroy initial row!
#ifdef FOR_DL
      fprintf (outdelfile, "%x %d", j, rowLength(newrows, i0));
      for (int k = 1; k <= rowLength(newrows, i0); k++)
          fprintf (outdelfile, " %x:%d", newrows[i0][k].id, newrows[i0][k].e);
      fprintf (outdelfile, "\n");
#endif
      free(newrows[i0]);
      newrows[i0] = NULL;

      if (index_data != NULL) // ie we want an index
      {
        index_data[i0].n = 0;
        free(index_data[i0].rels);
        index_data[i0].rels = NULL;
      }
    }
}


#define STRLENMAX 2048

// sparsemat is small_nrows x small_ncols, after small_ncols is found using
// renumbering.
static int
toFlush (const char *sparsename, typerow_t **sparsemat, int *colweight,
         int ncols, int small_nrows, int skip, int bin,
         const char *idealsfilename)
{
    unsigned long W;
    int small_ncols;

    fprintf(stderr, "Renumbering columns (including sorting w.r.t. weight)\n");
    renumber (&small_ncols, colweight, ncols, idealsfilename);

    fprintf(stderr, "small_nrows=%d small_ncols=%d\n",small_nrows,small_ncols);

    double tt = seconds();
    fprintf(stderr, "Writing sparse representation to file\n");
    W = flushSparse(sparsename, sparsemat, small_nrows, small_ncols, colweight, skip, bin);
    fprintf(stderr, "#T# writing sparse: %2.2lf\n", seconds()-tt);
    fprintf(stderr, "# Weight(M_small) = %lu\n", W);

    return small_ncols;
}

static void
build_newrows_from_file(typerow_t **newrows, FILE *hisfile, uint64_t bwcostmin,
                        const char* outdelfilename, index_data_t index_data)
{
    uint64_t bwcost;
    unsigned long addread = 0;
    char str[STRLENMAX];

    FILE *outdelfile = NULL;
    if (outdelfilename != NULL)
      {
        outdelfile = fopen (outdelfilename, "w+");
        if (outdelfile == NULL)
          {
            fprintf (stderr, "Error, cannot open file %s.\n", outdelfilename);
            exit(1);
          }
      }

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
	    doAllAdds(newrows, str, outdelfile, index_data);
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

/* if for_msieve=1, generate the *.cyc file needed by msieve to construct
   its matrix, which is of the following (binary) format:
      small_nrows
      n1 i1 i2 ... in1
      n2 j1 j2 ... jn2
      ...
      nk ...
   where each value is stored as a 32-bit integer (no linebreak),
   small_nrows is the number of relation-sets of the matrix,
   n1 is the number of relations in the first relation-set,
   i1 is the index of the first relation in the first relation-set
   (should correspond to line i1+2 in *.purged.gz), and so on */

// Feed sparsemat with M_purged
static void
readPurged (typerow_t **sparsemat, purgedfile_stream ps, int verbose,
            int for_msieve)
{
  if (for_msieve == 0)
    {
      fprintf(stderr, "Reading sparse matrix from purged file\n");
      for(int i = 0 ; purgedfile_stream_get(ps, NULL) >= 0 ; i++) {
        if (verbose && purgedfile_stream_disp_progress_now_p(ps))
          fprintf(stderr, "Treating old rel #%d at %2.2lf\n",
                  ps->rrows,ps->dt);
        if(ps->nc == 0)
          fprintf(stderr, "Hard to believe: row[%d] is NULL\n", i);
      qsort(ps->cols, ps->nc, sizeof(int), cmp);
      sparsemat[i] = (typerow_t *)malloc((ps->nc+1) * sizeof(typerow_t));
      ASSERT_ALWAYS(sparsemat[i] != NULL);
      int j = 1;
      for(int k = 0; k < ps->nc; k++)
        {
#ifdef FOR_DL
          if (j > 1 && sparsemat[i][j-1].id == ps->cols[k])
            sparsemat[i][j-1].e++;
          else
#endif
          {
            setCell(sparsemat[i][j], ps->cols[k], 1);
            j++;
          }
        }
      rowLength(sparsemat, i) = j-1;
      }
    }
  else /* for_msieve */
    {
      /* to generate the .cyc file for msieve, we only need to start from the
         identity matrix with newnrows relation-sets, where relation-set i
         contains only relation i. Thus we only need to read the first line of
         the purged file, to get the number of relations-sets. */

      long i;

      for (i = 0; i < ps->nrows; i++)
        {
          sparsemat[i] = (typerow_t *) malloc(2 * sizeof(typerow_t));
          rowCell(sparsemat, i, 1) = i;
          rowLength(sparsemat, i) = 1;
        }
    }
}

static void 
writeIndex(const char *indexname, index_data_t index_data,
        int small_nrows, int small_ncols)
{
    FILE *indexfile;
    indexfile = fopen_maybe_compressed(indexname, "w");
    fprintf(indexfile, "%d %d\n", small_nrows, small_ncols);

    for (int i = 0; i < small_nrows; ++i) {
        ASSERT (index_data[i].n > 0);
        fprintf(indexfile, "%d", index_data[i].n);
        for (int j = 0; j < index_data[i].n; ++j) {
#ifdef FOR_DL
            fprintf(indexfile, " " PURGE_INT_FORMAT ":%d",
                    index_data[i].rels[j].ind_row,
                    index_data[i].rels[j].e);
#else
            ASSERT (index_data[i].rels[j].e == 1);
            fprintf(indexfile, " " PURGE_INT_FORMAT,
                    index_data[i].rels[j].ind_row);
#endif
        }
        fprintf(indexfile, "\n");
    }
    fclose_maybe_compressed(indexfile, indexname);
}

void
generate_cyc (const char *outname, typerow_t **rows, uint32_t nrows)
{
  FILE *outfile;
  uint32_t t, u, i, k;

  outfile = fopen (outname, "w");
  ASSERT_ALWAYS(outfile != NULL);

  /* first write the number of relations as a 32-bit integer */
  t = nrows;
  fwrite (&t, sizeof(uint32_t), 1, outfile);

  /* then for each relation-set write a 32-bit integer giving the number k
     of its element, followed by those k elements */
  for (i = 0; i < nrows; i++)
    {
      t = rowLength(rows, i);
      ASSERT_ALWAYS(t != 0);
      fwrite (&t, sizeof(uint32_t), 1, outfile);
      for (k = 1; k <= t; k++)
        {
          u = rowCell(rows, i, k);
          fwrite (&u, sizeof(uint32_t), 1, outfile);
        }
    }

  fclose (outfile);
}

static void
fasterVersion(typerow_t **newrows, const char *sparsename,
              const char *indexname, const char *hisname, purgedfile_stream ps,
              uint64_t bwcostmin, int nrows, int ncols, int skip, int bin,
              const char *idealsfilename, const char *outdelfilename,
              int for_msieve)
{
    FILE *hisfile;
    int *colweight;
    int small_nrows, small_ncols;
    char str[STRLENMAX], *rp MAYBE_UNUSED;
    index_data_t index_data = NULL;

    hisfile = fopen (hisname, "r");
    ASSERT_ALWAYS(hisfile != NULL);
    /* read first line */
    rp = fgets (str, STRLENMAX, hisfile);

    // 1st pass: read the relations in *.purged and put them in newrows[][]
    readPurged(newrows, ps, 1, for_msieve);
#if DEBUG >= 1
    // [PG]: FIXME: seems to be broken for FFS ?
    for(int i = 0; i < nrows; i++){
	printf("row[%d]=", i);
	for(int k = 1; k <= newrows[i][0]; k++)
	    printf(" %d", newrows[i][k]);
	printf("\n");
    }
#endif

    if (indexname != NULL) {
        // At the beginning, the index_data consists of relsets that
        // are just single relations.
        index_data = (relset_t *) malloc (nrows*sizeof(relset_t));
        for (int i = 0; i < nrows; ++i) {
            index_data[i].n = 1;
            index_data[i].rels = (multirel_t *)malloc(sizeof(multirel_t));
            index_data[i].rels[0].ind_row = i;
            index_data[i].rels[0].e = 1;
        }
    }
    else
      index_data = NULL;

    // read merges in the *.merge.his file and replay them
    build_newrows_from_file(newrows, hisfile, bwcostmin, outdelfilename,
            index_data);

    /* compute column weights */
    colweight = (int*) malloc (ncols * sizeof(int *));
    ASSERT_ALWAYS(colweight != NULL);
    memset (colweight, 0, ncols * sizeof(int *));
    for (int i = 0; i < nrows; i++)
      if (newrows[i] != NULL)
        for(int k = 1; k <= rowLength(newrows, i); k++)
          colweight[rowCell(newrows, i, k)] += 1;

    /* crunch empty rows */
    for (int i = small_nrows = 0; i < nrows; i++)
      if (newrows[i] != NULL)
        newrows[small_nrows++] = newrows[i];
    if (indexname != NULL) {
        int ii = 0;
        for (int i = 0; i < nrows; i++) 
            if (index_data[i].n > 0) 
                index_data[ii++] = index_data[i];
            else
                free(index_data[i].rels);
        ASSERT (ii == small_nrows);
    }

#if defined FOR_DL && defined STAT_DL
    uint64_t count[11] = {0,0,0,0,0,0,0,0,0,0,0};
    uint64_t nonzero = 0;
    for (int i = 0; i < small_nrows ; i++)
      {
        for(int k = 1; k <= rowLength(newrows, i); k++)
          {
            if (abs(newrows[i][k].e) > 10)
              count[0]++;
            else
              count[abs(newrows[i][k].e)]++;
            nonzero++;
          }
      }
    fprintf (stderr, "# of non zero coeff: %lu\n", nonzero);
    for (int i = 1; i <= 10 ; i++)
      fprintf (stderr, "# of %d: %lu(%.2f%%)\n", i, count[i],
                       100 * (double) count[i]/nonzero);
    fprintf (stderr, "# of > 10: %lu(%.2f%%)\n", count[0],
                     100 * (double) count[0]/nonzero);
#endif

    if (for_msieve)
      {
        /* generate the <dat_file_name>.cyc file in "indexname" */
        if (skip != 0)
          {
            fprintf (stderr, "Error, skip should be 0 with --for_msieve\n");
            exit (1);
          }
        if (bin != 0)
          {
            fprintf (stderr, "Error, --binary incompatible with --for_msieve\n");
            exit (1);
          }
        generate_cyc (sparsename, newrows, small_nrows);
      }
    else
      {
        /* renumber columns after sorting them by decreasing weight */
        small_ncols = toFlush(sparsename, newrows, colweight, ncols,
                              small_nrows, skip, bin, idealsfilename);
        // Create the index
        if (indexname != NULL) 
          writeIndex(indexname, index_data, small_nrows, small_ncols);
      }

    // Free.
    free (colweight);
    for(int i = 0; i < small_nrows; i++)
        free (newrows[i]);
    free (newrows);
    if (index_data != NULL) {
        for (int i = 0; i < small_nrows; ++i) 
            free (index_data[i].rels);
        free (index_data);
    }

    fclose (hisfile);
}

static void
usage (const char *argv0)
{
  fprintf (stderr, "Usage: %s [options]\n", argv0);
  fprintf (stderr, "\nMandatory command line options: \n");
  fprintf (stderr, "   -purged  xxx   - input (purged) file is xxx\n");
  fprintf (stderr, "   -his  xxx      - history file (produced by merge)\n");
  fprintf (stderr, "   -out  xxx      - basename for output matrices\n");
  fprintf (stderr, "\nOther command line options: \n");
  fprintf (stderr, "   -skip nnn      - dense matrice contains the nnn heaviest "
                                       "columns (default %u)\n", SKIP_DEFAULT);
  fprintf (stderr, "   -v or --verbose\n");
  fprintf (stderr, "   --binary\n");
  fprintf (stderr, "   --noindex\n");
  fprintf (stderr, "   -index <file>  - if and only if there is no --noindex\n");
  fprintf (stderr, "   -ideals\n");
  fprintf (stderr, "   -outdel\n");
  exit (1);
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
  char * argv0 = argv[0];
    FILE *hisfile;
    uint64_t bwcostmin = 0;
    int nrows, ncols;
    typerow_t **newrows;
    int verbose = 0;
    int bin=0;
    int skip = SKIP_DEFAULT;
    int noindex = 0;
    int for_msieve = 0;
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
    param_list_configure_switch(pl, "--verbose", &verbose);
    param_list_configure_switch(pl, "--binary", &bin);
    param_list_configure_switch(pl, "--noindex", &noindex);
    param_list_configure_switch(pl, "--for_msieve", &for_msieve);
    param_list_configure_alias(pl, "--verbose", "-v");

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage(argv0);
    }
    const char * purgedname = param_list_lookup_string(pl, "purged");
    const char * hisname = param_list_lookup_string(pl, "his");
    const char * sparsename = param_list_lookup_string(pl, "out");
    const char * indexname = param_list_lookup_string(pl, "index");
    const char * idealsfilename = param_list_lookup_string(pl, "ideals");
    const char * outdelfilename = param_list_lookup_string(pl, "outdel");
    param_list_parse_int(pl, "skip", &skip);
    param_list_parse_uint64(pl, "bwcostmin", &bwcostmin);

    if (purgedname == NULL)
    {
      fprintf (stderr, "Error, missing mandatory -purged option.\n");
      usage(argv0);
    }

    if (sparsename == NULL)
    {
      fprintf (stderr, "Error, missing mandatory -out option.\n");
      usage(argv0);
    }

    if (has_suffix(sparsename, ".bin") || has_suffix(sparsename, ".bin.gz"))
        bin=1;

    if (noindex && indexname != NULL) {
        fprintf (stderr, "Error: --noindex was switched on, but a "
                "name for the index file was given.\n");
        exit (1);
    }

    if (noindex == 0 && indexname == NULL) {
        fprintf (stderr, "Error: --noindex was not given, but no "
                "index file was given.\n");
        exit (1);
    }

#ifdef FOR_DL
    if (skip != 0)
      {
        fprintf (stderr, "Error, for DL -skip should be 0\n");
        exit (1);
      }
    if (idealsfilename == NULL)
      {
        fprintf (stderr, "Error, for DL -ideals should be non null\n");
        exit (1);
      }
    if (outdelfilename == NULL)
      {
        fprintf (stderr, "Error, for DL -outdel should be non null\n");
        exit (1);
      }
#endif

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

    newrows = (typerow_t **)malloc(nrows * sizeof(typerow_t *));
    ASSERT_ALWAYS(newrows != NULL);

    // at the end of the following operations, newrows[i] is either
    // NULL
    // or k i_1 ... i_k which means that M_small will contain a row formed
    // of the addition of the rows of indices i_1 ... i_k in the original
    // matrix

    fasterVersion (newrows, sparsename, indexname, hisname, ps,
                   bwcostmin, nrows, ncols, skip, bin, 
                   idealsfilename, outdelfilename, for_msieve);


    purgedfile_stream_closefile(ps);
    purgedfile_stream_clear(ps);

    param_list_clear(pl);

    print_timing_and_memory (wct0);

    return 0;
}
