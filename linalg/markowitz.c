#include "utils/utils.h"
#include "sparse.h"
#include "merge_mono.h"
#include "markowitz.h"

#ifdef USE_MARKOWITZ

#define MKZ_DEBUG 1

// Again, a priority queue as a heap...!
// Q[0] contains the number of items in Q[], so that useful part of Q
// is Q[1..Q[0]]

#define MkZIsQueueEmpty(Q) (Q[0] == 0)

// get j-th component of Q[i] 
#define MkzGet(Q, i, j) (Q[((i)<<1)+(j)])
#define MkzSet(Q, i, j, val) (Q[((i)<<1)+(j)] = (val))

void
MkzPrintQueue(INT *Q)
{
    int level = 0, imax = 1, i;

    fprintf(stderr, "L0:");
    for(i = 1; i <= Q[0]; i++){
	fprintf(stderr, " %d", MkzGet(Q, i, 1));
	if(i == imax){
	    imax = (imax<<1)+1;
	    fprintf(stderr, "\nL%d:", ++level);
	}
    }
    fprintf(stderr, "\n");
}

void
MkzUpQueue(INT *Q, INT k)
{
    INT x = MkzGet(Q, k, 0), v = MkzGet(Q, k, 1);

    while((k > 1) && (MkzGet(Q, k/2, 1) >= v)){
	// we are at level > 0 and the father is >= son
	// the father replaces the son
	MkzSet(Q, k, 0, MkzGet(Q, k/2, 0));
	MkzSet(Q, k, 1, MkzGet(Q, k/2, 1)); // could be simplified...!
	k /= 2;
    }
    // we found the place of (x, v)
    MkzSet(Q, k, 0, x);
    MkzSet(Q, k, 1, v);
}

void
MkzInsert(INT *Q, INT j, INT c)
{
    Q[0]++;
    MkzSet(Q, Q[0], 0, j); 
    MkzSet(Q, Q[0], 1, c); 
    MkzUpQueue(Q, Q[0]);
}

// Move Q[1] down.
void
MkzDownQueue(INT *Q, INT k)
{
    INT x = MkzGet(Q, k, 0), v = MkzGet(Q, k, 1), j;

    while(k <= Q[0]/2){
	// k has at least a left son
	j = 2*k;
	if(j < Q[0])
	    // k has a right son
	    if(MkzGet(Q, j, 1) > MkzGet(Q, j+1, 1))
		j++;
	// at this point, Q[j] is the largest son
	if(v <= MkzGet(Q, j, 1))
	    break;
	else{
	    // the father takes the place of the son
	    MkzSet(Q, k, 0, MkzGet(Q, j, 0));
	    MkzSet(Q, k, 1, MkzGet(Q, j, 1));
	    k = j;
	}
    }
    // we found the place of v
    MkzSet(Q, k, 0, x);
    MkzSet(Q, k, 1, v);
}

void
MkzPopQueue(INT *j, INT *mkz, INT *Q)
{
    *j = MkzGet(Q, 1, 0);
    *mkz = MkzGet(Q, 1, 1);
    MkzSet(Q, 1, 0, MkzGet(Q, Q[0], 0));
    MkzSet(Q, 1, 1, MkzGet(Q, Q[0], 1));
    Q[0]--;
    MkzDownQueue(Q, 1);
}

void
MkzInit(sparse_mat_t *mat)
{
    INT i, j, k, mkz;

    printf("Entering initMarkowitz\n");
    mat->MKZQ = (INT *)malloc((mat->ncols+1) * 2 * sizeof(INT));
    // just to understand
    for(j = mat->jmin; j < mat->jmax; j++)
	if(mat->wt[GETJ(mat, j)] > 0){
	    mkz = mat->nrows;
#if MKZ_DEBUG >= 2
	    printf("j=%d wt=%d", j, mat->wt[GETJ(mat, j)]);
#endif
	    for(k = 1; k <= mat->R[j][0]; k++)
		if((i = mat->R[j][k]) != -1){
		    // this should be the weight of row i
#if MKZ_DEBUG >= 2
		    printf(" %d", mat->rows[i][0]);
#endif
		    if(mat->rows[i][0] < mkz)
			mkz = mat->rows[i][0];
		}
	    mkz = (mkz-1) * mat->wt[GETJ(mat, j)];
#if MKZ_DEBUG >= 2
	    printf(" => mkz=%d\n", mkz);
#endif
	    MkzInsert(mat->MKZQ, j, mkz);
	}
    // TODO: underflow in mkz???
}

void
MkzClose(sparse_mat_t *mat)
{
    free(mat->MKZQ);
}

#endif
