/* swar.c --- routines dealing with the SWAR data structure

Copyright 2008-2009 Francois Morain.
Reviewed by Paul Zimmermann, February 2009.

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

#include "utils.h"

#include "merge_opts.h"
#include "sparse.h"
#include "dclist.h"
#include "sparse_mat.h"
#include "report.h"

#ifndef USE_MARKOWITZ

#include "swar.h"
#include "merge_mono.h"

/* Initializes the SWAR data structure.
   Inputs:
      mat->ncols - number of columns of matrix mat
      mat->wt[j] - weight of column j, for 0 <= j < mat->ncols
      mat->cwmax - weight bound
      mat->S[w]  - lists containing only one element (-1), 0 <= w <= mat->cwmax
   Outputs:
      mat->S[w]  - doubly-chained list with columns j of weight w
 */
void
initSWAR (sparse_mat_t *mat)
{
    /* mat->cwmax+2 to prevent bangs */
    /* FIXME: what is the purpose of S[cwmax+1]? */
    dclist *S = (dclist *) malloc ((mat->cwmax + 2) * sizeof(dclist));
    dclist *A = (dclist *) malloc ((mat->jmax - mat->jmin) * sizeof(dclist));
    int k;

    ASSERT_ALWAYS(S != NULL);
    ASSERT_ALWAYS(A != NULL);
    /* S[0] has a meaning at least for temporary reasons */
    for(k = 0; k <= mat->cwmax + 1; k++)
	S[k] = dclistCreate (-1);
    mat->S = S;
    mat->A = A;
}

/* free the memory allocated in mat */
void
closeSWAR (sparse_mat_t *mat)
{
  int k;

  for(k = 0; k <= mat->cwmax + 1; k++)
    dclistClear (mat->S[k]);
  free (mat->S);
  free (mat->A);
}

// dump for debugging reasons
void
printSWAR(sparse_mat_t *mat, int ncols)
{
    int j, w;
    int32_t k;

    fprintf(stderr, "===== S is\n");
    for(w = 0; w <= mat->cwmax; w++){
	fprintf(stderr, "  S[%d] ", w);
	dclistPrint(stderr, mat->S[w]);
	fprintf(stderr, "\n");
    }
    fprintf(stderr, "===== Wj is\n");
    for(j = 0; j < ncols; j++)
	fprintf(stderr, "  Wj[%d]=%d\n", j, mat->wt[GETJ(mat, j)]);
    fprintf(stderr, "===== R is\n");
#ifndef USE_COMPACT_R
    for(j = 0; j < ncols; j++){
	if(mat->R[GETJ(mat, j)] == NULL)
	    continue;
	fprintf(stderr, "  R[%d]=", j);
	for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
          fprintf(stderr, " %ld", (long int) mat->R[GETJ(mat, j)][k]);
	fprintf(stderr, "\n");
    }
#else
	fprintf(stderr, "R: NYI in printSWAR\n");
	exit(1);
#endif
}

void
printStatsSWAR (sparse_mat_t *mat)
{
    int w;

    for (w = 0; w <= mat->cwmax; w++)
      fprintf (stderr, "I found %d primes of weight %d\n",
	       dclistLength (mat->S[w]), w);
}

// remove j from the stack it belongs to.
static void
remove_j_from_S(sparse_mat_t *mat, int j)
{
    dclist dcl = mat->A[GETJ(mat, j)];

    if(dcl == NULL){
	fprintf(stderr, "Column %d already removed?\n", j);
	return;
    }
#if DEBUG >= 2
    int ind = mat->wt[GETJ(mat, j)];
    if(ind > mat->cwmax)
        ind = mat->cwmax+1;
    fprintf(stderr, "S[%d]_b=", ind);
    dclistPrint(stderr, mat->S[ind]); fprintf(stderr, "\n");
    fprintf(stderr, "dcl="); dclistPrint(stderr, dcl); fprintf(stderr, "\n");
#endif
    /* remove current cell from dcl */
    dclistRemove (dcl);
#if DEBUG >= 2
    fprintf(stderr, "S[%d]_a=", ind);
    dclistPrint(stderr, mat->S[ind]); fprintf(stderr, "\n");
#endif
}

void
remove_j_from_SWAR(sparse_mat_t *mat, int j)
{
    remove_j_from_S(mat, j);
    mat->A[GETJ(mat, j)] = NULL;
}

/* remove the cell (i,j), and updates matrix correspondingly.
   Note: A[j] contains the address of the cell in S[w] where j is stored.
   
   Updates:
   - mat->wt[j] (weight of column j)
   - mat->S[w] : cell j is removed, with w = weight(j)
   - mat->S[w-1] : cell j is added
   - A[j] : points to S[w-1] instead of S[w]

   We arrive here when mat->wt[j] > 0.

*/
void
decreaseColWeightSWAR(sparse_mat_t *mat, int32_t j)
{
    int ind;

    // update weight
#if DEBUG >= 1
    fprintf(stderr, "removeCellSWAR: moving j=%d from S[%d] to S[%d]\n",
	    j, mat->wt[GETJ(mat, j)], decrS(mat->wt[GETJ(mat, j)]));
#endif
    ind = mat->wt[GETJ(mat, j)] = decrS(mat->wt[GETJ(mat, j)]);
    remove_j_from_S(mat, j);
    if(mat->wt[GETJ(mat, j)] > mat->cwmax)
	ind = mat->cwmax+1;
    // update A[j]
#if DEBUG >= 2
    fprintf(stderr, "S[%d]_b=", ind);
    dclistPrint(stderr, mat->S[ind]); fprintf(stderr, "\n");
#endif
    // TODO: replace this with a move of pointers...!
    mat->A[GETJ(mat, j)] = dclistInsert(mat->S[ind], j);
#if DEBUG >= 2
    fprintf(stderr, "S[%d]_a=", ind);
    dclistPrint(stderr, mat->S[ind]); fprintf(stderr, "\n");
#endif
}

// Returns a value < 0 if nothing was done, since the weight was too heavy.
int
addColSWAR(sparse_mat_t *mat, int32_t j)
{
    int ind;

    // update weight
#if DEBUG >= 1
    fprintf(stderr, "addColSWAR: moving j=%d from S[%d] to S[%d]?\n",
	    j, mat->wt[GETJ(mat, j)], incrS(mat->wt[GETJ(mat, j)]));
#endif
    ind = mat->wt[GETJ(mat, j)] = incrS(mat->wt[GETJ(mat, j)]);
    if (ind < 0)
	return ind;
    remove_j_from_S(mat, j);
    if(mat->wt[GETJ(mat, j)] > mat->cwmax){
#if DEBUG >= 1
	fprintf(stderr, "WARNING: column %d is too heavy (%d)\n", j,
		mat->wt[GETJ(mat, j)]);
#endif
	ind = mat->cwmax+1; // trick
    }
    // update A[j]
    mat->A[GETJ(mat, j)] = dclistInsert(mat->S[ind], j);
    return ind;
}

int
deleteEmptyColumnsSWAR (sparse_mat_t *mat)
{
    dclist dcl = mat->S[0];
    int njrem = 0;
    int32_t j;

    ASSERT(dcl != NULL);
    while (dcl->next != NULL)
      {
        j = dclistRemoveNext (dcl->next);
	njrem++;
	mat->A[GETJ(mat, j)] = NULL;
	mat->wt[GETJ(mat, j)] = 0;
	freeRj (mat, j);
      }
    mat->rem_ncols -= njrem;
    return njrem;
}

#endif /* ifndef USE_MARKOWITZ */
