/*
 * Program: main program for merge
 * Author : F. Morain
 * Purpose: merging relations
 * 
 * Algorithm: digest and interpolation from Cavallar.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h> /* for INT_MAX */

#include "utils/utils.h"
#include "files.h"
#include "gzip.h"
#include "sparse.h"
#include "merge_mono.h"
#include "prune.h"

#ifdef USE_MPI
#include "mpi.h"
#include "merge_mpi.h"
#endif

int
main(int argc, char *argv[])
{
    FILE *purgedfile;
    sparse_mat_t mat;
    report_t rep;
    char *purgedname = NULL, *hisname = NULL, *outname = NULL;
    int nrows, ncols;
    int cwmax = 20, rwmax = 1000000, maxlevel = 2, iprune = 0, keep = 128;
    int verbose = 0; /* default verbose level */
    double tt;
    double kprune = 1.0; /* prune keeps kprune * (initial excess) */
    double ratio = 1.1; /* bound on cN_new/cN to stop the computation */
    int i, forbw = 0, coverNmax = 0;

    fprintf (stderr, "%s.r%s", argv[0], REV);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif

#if TEX
    fprintf(stderr, "\\begin{verbatim}\n");
#endif    
    while(argc > 1 && argv[1][0] == '-'){
        if (argc > 2 && strcmp (argv[1], "-mat") == 0){
	    purgedname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-rebuild") == 0){
	    hisname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-cwmax") == 0){
	    cwmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-rwmax") == 0){
	    rwmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-maxlevel") == 0){
	    maxlevel = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-prune") == 0){
	    kprune = strtod(argv[2], NULL);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-keep") == 0){
	    keep = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 1 && strcmp (argv[1], "-v") == 0){
            verbose ++;
	    argc -= 1;
	    argv += 1;
	}
	else if (argc > 2 && strcmp (argv[1], "-forbw") == 0){
	    forbw = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-ratio") == 0){
	    ratio = strtod(argv[2], NULL);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-coverNmax") == 0){
	    coverNmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-out") == 0){
	    outname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else 
	    break;
    }
    purgedfile = gzip_open(purgedname, "r");
    ASSERT_ALWAYS(purgedfile != NULL);
    fscanf(purgedfile, "%d %d", &nrows, &ncols);

    mat.nrows = nrows;
    mat.ncols = ncols;
    mat.delta = keep; // FIXME: change the name of delta
    mat.cwmax = cwmax;
    mat.rwmax = rwmax;
    mat.mergelevelmax = maxlevel;
    
#ifdef USE_MPI
    mpi_start_proc(outname,&mat,purgedfile,purgedname,forbw,ratio,coverNmax);
    // TODO: clean the mat data structure (?)
    MPI_Finalize();
    return 0;
#endif    
    tt = seconds();
    initMat(&mat, 0, ncols);
    fprintf(stderr, "Time for initMat: %2.2lf\n", seconds()-tt);

    tt = seconds();
    initWeightFromFile(&mat, purgedfile);
    fprintf(stderr, "Time for initWeightFromFile: %2.2lf\n", seconds()-tt);
    gzip_close(purgedfile, purgedname);
    
    tt = seconds();
    fillSWAR(&mat);
    fprintf(stderr, "Time for fillSWAR: %2.2lf\n", seconds()-tt);
    
    purgedfile = gzip_open(purgedname, "r");
    ASSERT_ALWAYS(purgedfile != NULL);
    readmat(&mat, purgedfile);
    gzip_close(purgedfile, purgedname);
#if DEBUG >= 3
    checkmat(&mat);
#endif

    init_rep(&rep, outname, &mat, 0);
    report2(&rep, mat.nrows, mat.ncols);

    /* iprune is the excess we want at the end of prune */
    iprune = (mat.nrows-mat.ncols) * kprune;
    if (iprune < mat.delta) /* ensures iprune >= DELTA */
      iprune = mat.delta;
    /* only call prune if the current excess is larger than iprune */
    if (iprune < mat.nrows - mat.ncols){
	double tt = seconds();
	prune(&rep, &mat, iprune);
	fprintf(stderr, "Pruning: nrows=%d ncols=%d %2.2lf\n",
		mat.rem_nrows, mat.rem_ncols, seconds()-tt);
    }

    // ouhhhhhhhhhh
    // we do not have a clear idea of which function to minimize...!
#if 0
    forbw = 0;
    fprintf(stderr, "WARNING: forcing forbw=0...!!!!\n");
#endif

#if M_STRATEGY <= 2
    merge(&mat, maxlevel, verbose, forbw);
#else
    mergeOneByOne(&rep, &mat, maxlevel, verbose, forbw, ratio, coverNmax);
#endif
    gzip_close(rep.outfile, outname);
    fprintf(stderr, "Final matrix has w(M)=%lu, ncols*w(M)=%lu\n",
	    mat.weight, ((unsigned long)mat.rem_ncols) * mat.weight);
#if TEX
    fprintf(stderr, "\\end{verbatim}\n");
#endif
#if DEBUG >= 1
    {
	FILE *ofile = fopen("toto", "w");
	dumpSparse(ofile, &mat);
	fclose(ofile);
    }
#endif
    closeSWAR(/*&mat*/);
    return 0;
}
