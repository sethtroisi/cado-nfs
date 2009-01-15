#include "utils/utils.h"
#include "files.h"
#include "gzip.h"
#include "sparse.h"
#include "dclist.h"
#include "sparse_mat.h"
#include "report.h"

void
init_rep(report_t *rep, char *outname, sparse_mat_t *mat, int type, int bufsize)
{
    INT** tmp, i;

    rep->type = type;
    if(type == 2)
	// do nothing...!
	return;
    rep->outfile = gzip_open(outname, "w");
    switch(type){
    case 0:
	// classical one: output to a file
	break;
    case 1:
	// mostly for MPI
	tmp = (INT **)malloc(mat->nrows * sizeof(INT *));
	for(i = 0; i < mat->nrows; i++)
	    tmp[i] = NULL;
	rep->history = tmp;
	rep->mark = -1;
	rep->bufsize = bufsize;
	break;
    default:
	fprintf(stderr, "Unknown type: %d\n", type);
    }
}

// terrific hack: everybody on the same line
// the first is to be destroyed in replay!!!
void
reportn(report_t *rep, INT *ind, int n)
{
    int i;

#if DEBUG >= 1
    fprintf(stderr, "Reporting for n=%d\n", n);
#endif
    if(rep->type == 2)
	return;
    if(rep->type == 0){
	for(i = 0; i < n; i++){
	    fprintf(rep->outfile, "%ld", (long int) ind[i]);
	    if(i < n-1)
		fprintf(rep->outfile, " ");
	}
	fprintf(rep->outfile, "\n");
    }
    else if((rep->type == 1) && (rep->mark != -2)){
	// mark == -2 => we are probably resuming and we don't care
	// to consume a lot of memory that will not be used, anyway.
	rep->mark += 1;
	if(rep->history[rep->mark] == NULL)
	    rep->history[rep->mark] = 
		(INT *)malloc((rep->bufsize+1)*sizeof(INT));
	rep->history[rep->mark][0] = n;
	for(i = 0; i < n; i++)
	    rep->history[rep->mark][i+1] = ind[i];
    }
}

#if 0 // probably obsolete
// Same as reportn, but with another input type.
// TODO: destroy this, since a workaround is found using pointers.
void
reporthis(report_t *rep, INT tab[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1], int i0)
{
    int k;

#if DEBUG >= 1
    fprintf(stderr, "Reporting for tab[%d]\n", i0);
#endif
    if(rep->type == 0){
	for(k = 1; k <= tab[i0][0]; k++){
	    fprintf(rep->outfile, "%ld", (long int) tab[i0][k]);
	    if(k <= tab[i0][0]-1)
		fprintf(rep->outfile, " ");
	}
	fprintf(rep->outfile, "\n");
    }
    else{
	fprintf(stderr, "I told you this is useless\n");
	exit(0);
    }
}
#endif

void
report1(report_t *rep, INT i)
{
    reportn(rep, &i, 1);
}

void
report2(report_t *rep, INT i1, INT i2)
{
    INT tmp[2];
    tmp[0] = i1;
    tmp[1] = i2;
    reportn(rep, tmp, 2);
}

