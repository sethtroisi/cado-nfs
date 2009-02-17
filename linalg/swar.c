/* 
 * Program: history
 * Author : F. Morain
 * Purpose: data structure for elimination
 * 
 * Algorithm:
 *
 */

#include "utils/utils.h"

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
	dclistPrint(stderr, mat->S[w]->next);
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
	       dclistLength (mat->S[w]->next), w);
}

// remove j from the stack it belongs to.
static void
remove_j_from_S(sparse_mat_t *mat, int j)
{
    dclist dcl = mat->A[GETJ(mat, j)], foo;

    if(dcl == NULL){
	fprintf(stderr, "Column %d already removed?\n", j);
	return;
    }
#if DEBUG >= 2
    int ind = mat->wt[GETJ(mat, j)];
    if(ind > mat->cwmax)
        ind = mat->cwmax+1;
    fprintf(stderr, "S[%d]_b=", ind);
    dclistPrint(stderr, mat->S[ind]->next); fprintf(stderr, "\n");
    fprintf(stderr, "dcl="); dclistPrint(stderr, dcl); fprintf(stderr, "\n");
#endif
    foo = dcl->prev;
    foo->next = dcl->next;
    if(dcl->next != NULL)
	dcl->next->prev = foo;
#if DEBUG >= 2
    fprintf(stderr, "S[%d]_a=", ind);
    dclistPrint(stderr, mat->S[ind]->next); fprintf(stderr, "\n");
#endif
#if USE_CONNECT == 0
    free(dcl);
#endif
}

void
remove_j_from_SWAR(sparse_mat_t *mat, int j)
{
    remove_j_from_S(mat, j);
#if USE_CONNECT
    free(mat->A[GETJ(mat, j)]);
#endif
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
    dclistPrint(stderr, mat->S[ind]->next); fprintf(stderr, "\n");
#endif
    // TODO: replace this with a move of pointers...!
#if USE_CONNECT == 0
    mat->A[GETJ(mat, j)] = dclistInsert(mat->S[ind], j);
#else
    dclistConnect(mat->S[ind], mat->A[GETJ(mat, j)]);
#endif
#if DEBUG >= 2
    fprintf(stderr, "S[%d]_a=", ind);
    dclistPrint(stderr, mat->S[ind]->next); fprintf(stderr, "\n");
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
#if USE_CONNECT == 0
    mat->A[GETJ(mat, j)] = dclistInsert(mat->S[ind], j);
#else
    dclistConnect(mat->S[ind], mat->A[GETJ(mat, j)]);
#endif
    return ind;
}

int
deleteEmptyColumnsSWAR(sparse_mat_t *mat)
{
    dclist dcl = mat->S[0], foo;
    int njrem = 0;
    int32_t j;

    while(dcl->next != NULL){
	foo = dcl->next;
	j = foo->j;
	dcl->next = foo->next;
	free(foo);
	njrem++;
	mat->A[GETJ(mat, j)] = NULL;
	mat->wt[GETJ(mat, j)] = 0;
	freeRj(mat, j);
    }
    mat->rem_ncols -= njrem;
    return njrem;
}

#endif /* ifndef USE_MARKOWITZ */
