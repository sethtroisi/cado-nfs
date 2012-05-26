/* debug functions for merge and replay */

void
dumpSparse(FILE *ofile, filter_matrix_t *mat)
{
    long w = 0;
    int i, j, k, buf[1000], ibuf, new_nrows = 0;

    fprintf(ofile, "%d %d\n", mat->rem_nrows, mat->rem_ncols);
    for(i = 0; i < mat->nrows; i++){
	if(isRowNull(mat, i))
	    continue;
	w += lengthRow(mat, i);
	new_nrows++;
	ibuf = 0;
	for(k = 1; k <= lengthRow(mat, i); k++)
	    buf[ibuf++] = cell(mat, i, k);
	fprintf(ofile, "%d", ibuf);
	for(k = 0; k < ibuf; k++)
	    fprintf(ofile, " %d", buf[k]);
	fprintf(ofile, "\n");
    }
    fprintf(stderr, "Weight = %ld\n", w);
    ASSERT(new_nrows == mat->rem_nrows);
    // check j used
    for(j = 0; j < mat->ncols; j++)
	if(mat->wt[GETJ(mat, j)] != 0){
	    fprintf(ofile, "J %d %d\n", j, abs(mat->wt[GETJ(mat, j)]));
	    if(abs(mat->wt[GETJ(mat, j)]) <= mat->mergelevelmax)
		fprintf(stderr, "Gasp: %d %d\n", j, mat->wt[GETJ(mat, j)]);
	}
}

