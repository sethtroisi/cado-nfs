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

#include "merge_opts.h"
#include "sparse.h"
#include "dclist.h"
#include "sparse_mat.h"
#include "report.h"
#ifndef USE_MARKOWITZ
# include "swar.h"
#else
# include "markowitz.h"
#endif
#include "merge_mono.h"

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
    char *resumename = NULL;
    int nrows, ncols;
    int cwmax = 20, rwmax = 1000000, maxlevel = 2, keep = 128;
    int verbose = 0; /* default verbose level */
    double tt;
    double ratio = 1.1; /* bound on cN_new/cN to stop the computation */
    int i, forbw = 0, coverNmax = 0;
#ifdef USE_MARKOWITZ
    int wmstmax = 7; /* use real MST minimum for wt[j] <= wmstmax */
    int mkzrnd = 0;
#endif
    int itermax = 0;

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
	else if (argc > 2 && strcmp (argv[1], "-resume") == 0){
	    resumename = argv[2];
	    argc -= 2;
	    argv += 2;
	}
#ifdef USE_MARKOWITZ
	else if (argc > 2 && strcmp (argv[1], "-wmstmax") == 0){
	    wmstmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-mkzrnd") == 0){
	    mkzrnd = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
#endif
	else if (argc > 2 && strcmp (argv[1], "-itermax") == 0){
	    itermax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else
	  {
	    fprintf (stderr, "Error, unknown option %s\n", argv[1]);
	    exit (1);
	  }
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
    mpi_start_proc(outname,&mat,purgedfile,purgedname,forbw,ratio,coverNmax,
		   resumename);
    // TODO: clean the mat data structure (?)
    MPI_Finalize();
    return 0;
#endif    
    tt = seconds();
    initMat(&mat, 0, ncols);
    fprintf(stderr, "Time for initMat: %2.2lf\n", seconds()-tt);

    tt = seconds();
    initWeightFromFile(&mat, purgedfile, 1);
    fprintf(stderr, "Time for initWeightFromFile: %2.2lf\n", seconds()-tt);
    gzip_close(purgedfile, purgedname);

#ifndef USE_MARKOWITZ
    fprintf(stderr, "SWAR version\n");
    initSWAR(&mat);
#else
    fprintf(stderr, "Markowitz version\n");
#endif
    tt = seconds();
    fillmat(&mat);
    fprintf(stderr, "Time for fillmat: %2.2lf\n", seconds()-tt);
    
    purgedfile = gzip_open(purgedname, "r");
    ASSERT_ALWAYS(purgedfile != NULL);
    readmat(&mat, purgedfile, 1);
    gzip_close(purgedfile, purgedname);
#if DEBUG >= 3
    checkmat(&mat);
#endif

    init_rep(&rep, outname, &mat, 0, MERGE_LEVEL_MAX);
    report2(&rep, mat.nrows, mat.ncols);

    if(resumename != NULL)
	resume(&rep, &mat, resumename);

    // ouhhhhhhhhhh
    // we do not have a clear idea of which function to minimize...!
#if 0
    forbw = 0;
    fprintf(stderr, "WARNING: forcing forbw=0...!!!!\n");
#endif

    mat.itermax = itermax;
#ifdef USE_MARKOWITZ
    mat.wmstmax = wmstmax;
    mat.mkzrnd = mkzrnd;
    tt = seconds();
    MkzInit(&mat);
    fprintf(stderr, "Time for MkzInit: %2.2lf\n", seconds()-tt);
#endif

#if M_STRATEGY <= 2
# ifdef USE_MARKOWITZ
    fprintf(stderr, "merge NYI for Markowitz\n");
    return 1;
# else
    merge(&mat, maxlevel, verbose, forbw);
# endif
#else
    mergeOneByOne(&rep, &mat, maxlevel, verbose, forbw, ratio, coverNmax);
#endif

    gzip_close(rep.outfile, outname);
    fprintf(stderr, "Final matrix has N=%d nc=%d (%d) w(M)=%lu N*w(M)=%"PRIu64"\n",
	    mat.rem_nrows, mat.rem_ncols, mat.rem_nrows-mat.rem_ncols,
	    mat.weight,
	    (uint64_t) mat.rem_nrows * (uint64_t) mat.weight);
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
#ifndef USE_MARKOWITZ
    closeSWAR(/*&mat*/);
#else
    MkzClose(&mat);
#endif
    return 0;
}
