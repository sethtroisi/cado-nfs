/*
 * Program: main program for merge
 * Author : F. Morain
 * Purpose: postprocessing relations
 * 
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h> /* for INT_MAX */

#include "utils.h"
#include "files.h"
#include "gzip.h"

#include "merge_opts.h"

#include "sparse.h"
#include "dclist.h"
#include "sparse_mat.h"
#include "mst.h"
#include "report.h"
#ifndef USE_MARKOWITZ
# include "swar.h"
#else
# include "markowitz.h"
#endif
#include "merge_mono.h"

#define DEBUG 0

void
flushRow(FILE *outfile, sparse_mat_t *mat, int32_t i)
{
    int k;

    fprintf(outfile, "%d", lengthRow(mat, i));
    for(k = 1; k <= lengthRow(mat, i); k++){
	fprintf(outfile, " ");
	fprintf(outfile, PURGE_INT_FORMAT, cell(mat, i, k));
    }
    fprintf(outfile, "\n");
}

void
flushMatrix(FILE *outfile, sparse_mat_t *mat)
{
    int32_t i;

    for(i = 0; i < mat->nrows; i++)
	flushRow(outfile, mat, i);
}

void
replaceWithMST(sparse_mat_t *mat, int m, int32_t *ind)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX];
    int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];
    int father[MERGE_LEVEL_MAX], height[MERGE_LEVEL_MAX], hmax;
    int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];
    int k, wold, wnew;

    wold = 0;
    for(k = 0; k < m; k++)
	wold += lengthRow(mat, ind[k]);
    fillRowAddMatrix(A, mat, m, ind);
#if DEBUG >= 1
    fprintf(stderr, "Initial rows\n");
    for(k = 0; k < m; k++){
	fprintRow(stderr, mat->rows[ind[k]]);
	fprintf(stderr, "\n");
    }
#endif
    if(m == 2){
	int w01 = weightSum(mat, ind[0], ind[1]);

	wnew = w01 + lengthRow(mat, ind[1]);
#if DEBUG >= 1
	fprintf(stderr, "wold=%d wnew=%d\n", wold, wnew);
#endif
	if(wnew < wold){
	    addRowsAndUpdate(mat, ind[0], ind[1], w01);
#if DEBUG >= 1
	    fprintRow(stderr, mat->rows[ind[0]]);
	    fprintf(stderr, "\n");
#endif
	}
    }
    else{
	hmax = minimalSpanningTree(&wnew, father, height, sons, m, A);
	wnew += lengthRow(mat, ind[0]);
#if DEBUG >= 1
	fprintf(stderr, "wold=%d wnew=%d\n", wold, wnew);
#endif
	if(wnew < wold){
	    hmax=addFatherToSons(history,mat,m,ind,A,father,height,hmax,sons);
#if DEBUG >= 1
	    hmax = 0;
	    fprintf(stderr, "New rows\n");
	    for(k = 0; k < m; k++){
		fprintRow(stderr, mat->rows[ind[k]]);
		fprintf(stderr, "\n");
		hmax += lengthRow(mat, ind[k]);
	    }
	    fprintf(stderr, "wold=%d wnewer=%d\n", wold, hmax);
#endif
	}
    }
} 

void
postProcess(FILE *outfile, sparse_mat_t *mat)
{
    int32_t dj, j, mkz;
    int32_t ind[MERGE_LEVEL_MAX];
    unsigned long oldw = mat->weight;

    while(!MkzIsQueueEmpty(mat->MKZQ)){
	MkzPopQueue(&dj, &mkz, mat->MKZQ, mat->MKZA);
	j = dj + mat->jmin;
#if DEBUG >= 1
	fprintf(stderr, "I popped j=%d mkz=%d wt=%d\n", j, mkz, mat->wt[dj]);
#endif
	if(mat->wt[dj] == 0){
	    fprintf(stderr, "Strange: wt(%d)=0\n", j);
	} // wt = 1 useless
	else if(mat->wt[dj] >= 2){
#if DEBUG >= 1
	    fprintf(stderr, "wt(%d)=%d\n", j, mat->wt[dj]);
#endif
	    fillTabWithRowsForGivenj(ind, mat, j);
	    replaceWithMST(mat, mat->wt[dj], ind);
	}
    }
    flushMatrix(outfile, mat);
    fprintf(stderr, "old_index_weight//new_index_weight: %lu // %lu\n",
	    oldw, mat->weight);
}

int
main(int argc, char *argv[])
{
    FILE *indexfile, *outfile;
    char *indexname = NULL, *outname = NULL;
    sparse_mat_t mat;
    int small_nrows, small_ncols, nrows = 0;
    int i;
    double tt;

    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");
    while(argc > 1 && argv[1][0] == '-'){
        if (argc > 2 && strcmp (argv[1], "-index") == 0){
	    indexname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-out") == 0){
	    outname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-nrows") == 0){
	    nrows = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else
	  {
	    fprintf (stderr, "Error, unknown option %s\n", argv[1]);
	    exit (1);
	  }
    }
    ASSERT_ALWAYS(nrows > 0);
    indexfile = gzip_open(indexname, "r");
    ASSERT_ALWAYS(indexfile != NULL);
    int rc;
    rc = fscanf(indexfile, "%d %d", &small_nrows, &small_ncols);
    ASSERT_ALWAYS(rc == 2);

    // a typical line is 
    // nj j1 ... j_nj where 0 <= j_r < nrows and not small_ncols...! 
    mat.nrows = small_nrows;
    mat.ncols = nrows; // ohhhhhhhhhhhhh!!!!
    mat.wmstmax = 7;
    mat.cwmax = mat.nrows;
    mat.mkztype = 1;

    tt = seconds();
    initMat(&mat, 0, nrows);
    fprintf(stderr, "Time for initMat: %2.2lf\n", seconds()-tt);

    tt = seconds();
    initWeightFromFile(&mat, indexfile, 0);
    fprintf(stderr, "Time for initWeightFromFile: %2.2lf\n", seconds()-tt);
    gzip_close(indexfile, indexname);

    fprintf(stderr, "Markowitz version (mkztype=%d)\n", mat.mkztype);
    tt = seconds();
    fillmat(&mat);
    fprintf(stderr, "Time for fillmat: %2.2lf\n", seconds()-tt);
    
    indexfile = gzip_open(indexname, "r");
    ASSERT_ALWAYS(indexfile != NULL);
    readmat(&mat, indexfile, 0, 0, 0);
    gzip_close(indexfile, indexname);
#if DEBUG >= 1
    for(i = 0; i < mat.nrows; i++){
	fprintRow(stderr, mat.rows[i]);
	fprintf(stderr, "\n");
    }
#endif

    outfile = gzip_open(outname, "w");
    ASSERT_ALWAYS(outfile != NULL);
    fprintf(outfile, "%d %d\n", small_nrows, small_ncols);

    tt = seconds();
    MkzInit(&mat);
    fprintf(stderr, "Time for MkzInit: %2.2lf\n", seconds()-tt);

    postProcess(outfile, &mat);

    MkzClose(&mat);
    gzip_close(outfile, outname);
    return 0;
}
