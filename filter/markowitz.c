#include "cado.h"
#include "utils.h"

#include "merge_opts.h"

#ifdef USE_MARKOWITZ

#include "sparse.h"
#include "dclist.h"
#include "filter_matrix.h"
#include "mst.h"
#include "report.h"
#include "markowitz.h"

#define MKZ_DEBUG 0
#define MKZ_TIMINGS 1

#if MKZ_TIMINGS
double tmkzup, tmkzdown, tmkzupdown, tmkzcount;
#endif

#define MKZ_INF -1

// Again, a priority queue as a heap...!
// Q[0] contains the number of items in Q[], so that useful part of Q
// is Q[1..Q[0]]

// Q[2*i] contains j-jmin = dj
// Q[2*i+1] contains the Markowitz count for j

// get r-th component of Q[i]
#define MkzGet(Q, i, r) (Q[((i)<<1)+(r)])
#define MkzSet(Q, i, r, val) (Q[((i)<<1)+(r)] = (val))

static int
MkzGetCount(int32_t *Q, int32_t *A, int32_t dj) MAYBE_UNUSED;
static int
MkzGetCount(int32_t *Q, int32_t *A, int32_t dj)
{
    return MkzGet(Q, A[dj], 1);
}

int
MkzIsAlive(int32_t *A, int32_t dj)
{
    return A[dj] != MKZ_INF;
}

// (Q, A)[k1] <- (Q, A)[k2]
static void
MkzAssign(int32_t *Q, int32_t *A, int32_t k1, int32_t k2)
{
    int32_t dj = MkzGet(Q, k2, 0);

    MkzSet(Q, k1, 0, dj);
    MkzSet(Q, k1, 1, MkzGet(Q, k2, 1)); // could be simplified...!
    A[dj] = k1;
}

static void
MkzPrintQueue(int32_t *Q) MAYBE_UNUSED;
static void
MkzPrintQueue(int32_t *Q)
{
    int level = 0, imax = 1, i;

    fprintf(stderr, "L0:");
    for(i = 1; i <= Q[0]; i++){
	fprintf(stderr, " [%d, %d]", MkzGet(Q, i, 1), MkzGet(Q, i, 0));
	if(i == imax){
	    imax = (imax<<1)+1;
	    fprintf(stderr, "\nL%d:", ++level);
	}
    }
    fprintf(stderr, "\n");
}

static void
MkzUpQueue(int32_t *Q, int32_t *A, int32_t k)
{
    int32_t dj = MkzGet(Q, k, 0), count = MkzGet(Q, k, 1);
#if MKZ_TIMINGS
    double tt = seconds();
#endif

    while((k > 1) && (MkzGet(Q, k/2, 1) >= count)){
	// we are at level > 0 and the father is >= son
	// the father replaces the son
	MkzAssign(Q, A, k, k/2);
	k /= 2;
    }
    // we found the place of (dj, count)
    MkzSet(Q, k, 0, dj);
    MkzSet(Q, k, 1, count);
    A[dj] = k;
#if MKZ_TIMINGS
    tmkzup += (seconds()-tt);
#endif
}

static void
MkzInsert(int32_t *Q, int32_t *A, int32_t dj, int32_t count)
{
    Q[0]++;
    MkzSet(Q, Q[0], 0, dj);
    MkzSet(Q, Q[0], 1, count);
    A[dj] = Q[0];
    MkzUpQueue(Q, A, Q[0]);
}

// Move Q[k] down.
static void
MkzDownQueue(int32_t *Q, int32_t *A, int32_t k)
{
    int32_t dj = MkzGet(Q, k, 0), count = MkzGet(Q, k, 1), j;
#if MKZ_TIMINGS
    double tt = seconds();
#endif

    while(k <= Q[0]/2){
	// k has at least a left son
	j = 2*k;
	if(j < Q[0])
	    // k has a right son
	    if(MkzGet(Q, j, 1) > MkzGet(Q, j+1, 1))
		j++;
	// at this point, Q[j] is the largest son
	if(count <= MkzGet(Q, j, 1))
	    break;
	else{
	    // the father takes the place of the son
	    MkzAssign(Q, A, k, j);
	    k = j;
	}
    }
    // we found the place of (dj, count)
    MkzSet(Q, k, 0, dj);
    MkzSet(Q, k, 1, count);
    A[dj] = k;
#if MKZ_TIMINGS
    tmkzdown += (seconds()-tt);
#endif
}

// (Q, A)[k] has just arrived, but we have to move it in the heap, so that
// it finds its place.
static void
MkzMoveUpOrDown(int32_t *Q, int32_t *A, int32_t k)
{
#if MKZ_TIMINGS
    double tt = seconds();
#endif

    // move new node up or down
    if(k == 1)
	// rare event!
	MkzDownQueue(Q, A, 1);
    else{
	// k has a father
	if(MkzGet(Q, k/2, 1) > MkzGet(Q, k, 1))
	    // we have to move up
	    MkzUpQueue(Q, A, k);
	else
	    MkzDownQueue(Q, A, k);
    }
#if MKZ_TIMINGS
    tmkzupdown += (seconds()-tt);
#endif
}

// Remove (Q, A)[k].
static void
MkzDelete(int32_t *Q, int32_t *A, int32_t k)
{
#if MKZ_DEBUG >= 1
    fprintf(stderr, "MKZ: deleting (Q, A)[%d]=[%d, %d]\n", k,
	    MkzGet(Q, k, 0), MkzGet(Q, k, 1));
#endif
    // we put Q[Q[0]] in Q[k]
    MkzAssign(Q, A, k, Q[0]);
    Q[0]--;
    MkzMoveUpOrDown(Q, A, k);
}

// Remove (Q, A)[k].
void
MkzRemove(int32_t *dj, int32_t *mkz, int32_t *Q, int32_t *A, int32_t k)
{
    *dj = MkzGet(Q, k, 0);
    *mkz = MkzGet(Q, k, 1);
    A[*dj] = MKZ_INF;
    MkzAssign(Q, A, k, Q[0]);
    Q[0]--;
    MkzMoveUpOrDown(Q, A, k);
}

static int
MkzIsHeap(int32_t *Q) MAYBE_UNUSED;
static int
MkzIsHeap(int32_t *Q)
{
    int k;

    for(k = 1; k <= Q[0]/2; k++){
	// k has a left son
	if(MkzGet(Q, k, 1) > MkzGet(Q, 2*k, 1)){
	    fprintf(stderr, "Pb: father=%d > lson=%d\n",
		    MkzGet(Q, k, 1), MkzGet(Q, 2*k, 1));
	    return 0;
	}
	if(k < Q[0]/2){
	    // k has a right son
	    if(MkzGet(Q, k, 1) > MkzGet(Q, 2*k+1, 1)){
		fprintf(stderr, "Pb: father=%d > rson=%d\n",
			MkzGet(Q, k, 1), MkzGet(Q, 2*k+1, 1));
		return 0;
	    }
	}
    }
    return 1;
}

static void
MkzCheck(filter_matrix_t *mat)
{
    int32_t dj;

    for(dj = 0; dj < mat->jmax - mat->jmin; dj++)
	if(mat->wt[dj] > 0)
	    if(MkzGet(mat->MKZQ, mat->MKZA[dj], 0) != dj)
		fprintf(stderr, "GASP: %d in MkzCheck\n", dj);
}

static int
Cavallar(filter_matrix_t *mat, int32_t j)
{
    return abs(mat->wt[GETJ(mat, j)]);
}

static int
pureMkz(filter_matrix_t *mat, int32_t j)
{
    int mkz, k, i;

#if MKZ_TIMINGS
    double tt = seconds();
#endif
    // trick to be sure that columns with wt <= 2 are treated asap
    if(mat->wt[GETJ(mat, j)] == 1)
      mkz = 1 - 2 * mat->ncols;
#if 0 /* do not treat columns of weight 2 specially for the moment */
    else if(mat->wt[GETJ(mat, j)] == 2)
      {
        int32_t ind[MERGE_LEVEL_MAX];
	fillTabWithRowsForGivenj(ind, mat, j);
	// the more this is < 0, the less the weight is
	mkz = weightSum(mat, ind[0], ind[1])-2*mat->ncols;
      }
#endif
    else
      {
        // real traditional Markowitz count
        mkz = mat->nrows;
        for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
          if((i = mat->R[GETJ(mat, j)][k]) != -1){
	    // this should be the weight of row i
	    if(mat->rows[i][0] < mkz)
              mkz = mat->rows[i][0];
          }
        /* the lightest row has weight mkz, we add it (wt-1) times,
           remove it once, and remove wt entries in the jth column */
        mkz = (mkz - 1) * (mat->wt[GETJ(mat, j)] - 2);
      }
#if MKZ_TIMINGS
    tmkzcount += (seconds()-tt);
#endif
    return mkz;
}

// forcing lighter columns first.
static int
lightColAndMkz(filter_matrix_t *mat, int32_t j)
{
    int mkz, k, i, wj, cte;
    int32_t ind[MERGE_LEVEL_MAX];
    double tfill, tMST;

#if MKZ_TIMINGS
    double tt = seconds();
#endif
    // trick to be sure that columns with wt <= 2 are treated asap
    wj = mat->wt[GETJ(mat, j)];
    cte = mat->ncols; // could be smaller
    if(wj == 1)
	return cte;
    else if(wj == 2){
	fillTabWithRowsForGivenj(ind, mat, j);
	// the more this is < 0, the less the weight is
	return 2 * cte + weightSum(mat, ind[0], ind[1]);
    }
    else if(wj <= mat->wmstmax){
	fillTabWithRowsForGivenj(ind, mat, j);
	mkz = minCostUsingMST(mat, mat->wt[GETJ(mat, j)], ind, &tfill, &tMST);
	return wj * cte + mkz;
    }
    // real traditional Markowitz count
    mkz = mat->nrows;
    for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
	if((i = mat->R[GETJ(mat, j)][k]) != -1){
	    // this should be the weight of row i
	    if(mat->rows[i][0] < mkz)
		mkz = mat->rows[i][0];
	}
    mkz = wj * cte + (mkz-1) * (wj-1);
#if MKZ_TIMINGS
    tmkzcount += (seconds()-tt);
#endif
    return mkz;
}

static int
MkzCount(filter_matrix_t *mat, int32_t j)
{
    switch(mat->mkztype){
    case 2:
	return lightColAndMkz(mat, j);
    case 1:
	return pureMkz(mat, j);
    case 0:
    default:
	return Cavallar(mat, j);
    }
}

void
MkzPopQueue(int32_t *dj, int32_t *mkz,  filter_matrix_t *mat)
{
  int32_t *Q = mat->MKZQ;
  int32_t *A = mat->MKZA;

    *dj = MkzGet(Q, 1, 0);
    *mkz = MkzGet(Q, 1, 1);
    if (MkzCount (mat, *dj) != *mkz)
      {
        fprintf (stderr, "Error, cost do not match for j=%d, expected %d, got %d\n", *dj, MkzCount (mat, *dj), *mkz);
        exit (1);
      }
#if 0 // to see what happens
    {
	int i;

	for(i = 2; i <= Q[0]; i++)
	    if(MkzGet(Q, i, 1) != *mkz)
		break;
	printf("N(mkz=%d)=%d\n", *mkz, i);
    }
#endif
    A[*dj] = MKZ_INF;
    MkzAssign(Q, A, 1, Q[0]);
    Q[0]--;
    MkzDownQueue(Q, A, 1);
}

void
MkzInit(filter_matrix_t *mat)
{
    int32_t j, mkz;
    int sz = 0;

#if MKZ_TIMINGS
    tmkzup = tmkzdown = tmkzupdown = tmkzcount = 0.0;
#endif
    fprintf(stderr, "Entering initMarkowitz");
    fprintf(stderr, " (wmstmax=%d, rnd=%d", mat->wmstmax, mat->mkzrnd);
    fprintf(stderr, ", type=%d)\n", mat->mkztype);
    // compute number of elligible columns in the heap
    for(j = mat->jmin; j < mat->jmax; j++)
	if(mat->wt[GETJ(mat, j)] > 0)
	    sz++;
    fprintf(stderr, "Allocating heap for %d columns\n", sz);
    mat->MKZQ = (int32_t *)malloc((sz+1) * 2 * sizeof(int32_t));
    mat->MKZQ[0] = 0;
    mat->MKZQ[1] = sz; // why not?
    // every j needs a pointer
    mat->MKZA = (int32_t *)malloc((mat->jmax - mat->jmin + 1) * sizeof(int32_t));
    for(j = mat->jmin; j < mat->jmax; j++)
	if(mat->wt[GETJ(mat, j)] > 0){
	    mkz = MkzCount(mat, j);
#if MKZ_DEBUG >= 1
	    printf("j=%d wt=%d", j, mat->wt[GETJ(mat, j)]);
	    printf(" => mkz=%d\n", mkz);
#endif
	    MkzInsert(mat->MKZQ, mat->MKZA, GETJ(mat, j), mkz);
	}
        else
	    mat->MKZA[GETJ(mat, j)] = MKZ_INF;
    // TODO: overflow in mkz???
    MkzCheck(mat);
#if MKZ_DEBUG >= 1
    fprintf(stderr, "Initial queue is\n");
    MkzPrintQueue(mat->MKZQ);
#endif
}

void
MkzClose(filter_matrix_t *mat)
{
    fprintf(stderr, "Max Markowitz count: %d\n",
	    MkzGet(mat->MKZQ, mat->MKZQ[0], 1));
#if MKZ_TIMINGS
    fprintf(stderr, "MKZT: up=%d down=%d updown=%d count=%d\n",
	    (int)tmkzup, (int)tmkzdown, (int)tmkzupdown, (int)tmkzcount);
#endif
    free(mat->MKZQ);
    free(mat->MKZA);
}

int
MkzIncrCol(filter_matrix_t *mat, int32_t j)
{
    int ind;

#if MKZ_DEBUG >= 1
    fprintf(stderr, "Incr: wt(%d) was %d\n", j, mat->wt[GETJ(mat, j)]);
#endif
    ind = mat->wt[GETJ(mat, j)] = incrS(mat->wt[GETJ(mat, j)]);
    return ind;
}

// Row[i] has been adjoined to column j, so that we can incrementally
// change the Markowitz count.
void
MkzUpdate(filter_matrix_t *mat, int32_t i MAYBE_UNUSED, int32_t j)
{
    int32_t adr = mat->MKZA[GETJ(mat, j)];
    int mkz;

    if(adr == -1){
#if MKZ_DEBUG >= 1
	fprintf(stderr, "Prevented use of adr[%d]=-1 in MkzUpdate\n", j);
#endif
	return;
    }
#if MKZ_DEBUG >= 1
    if((mat->wt[GETJ(mat, j)] == 0) || (mat->wt[GETJ(mat, j)] == 1))
	fprintf(stderr, "W: wt[%d] = %d\n", j, mat->wt[GETJ(mat, j)]);
#endif
#if 1
    // costly?
    mkz = MkzCount(mat, j);
#else
    // old_count = min (r_ii-1)*(w-2) = mu * (w-2)
    // new_count = min(mu, r_i-1)*(w-1)
    mkz = MkzGet(mat->MKZQ, adr, 1)/(mat->wt[j]-2); // mu
    if(mat->rows[i][0] < mkz+1)
	mkz = mat->rows[i][0]-1;
    mkz *= (mat->wt[j]-1);
#endif
#if MKZ_DEBUG >= 1
    fprintf(stderr, "Updating j=%d (old=%d, new=%d)\n", j,
	    MkzGet(mat->MKZQ, adr, 1), mkz);
#endif
    // nothing to do if new count == old count
    if(mkz < MkzGet(mat->MKZQ, adr, 1)){
	// add new count
	MkzSet(mat->MKZQ, adr, 1, mkz);
	// a variant of delete is needed...!
	MkzMoveUpOrDown(mat->MKZQ, mat->MKZA, adr);
    }
}

/*
   Updates:
   - mat->wt[j] (weight of column j)

   We arrive here when mat->wt[j] > 0.

*/
void
MkzDecreaseColWeight(filter_matrix_t *mat, int32_t j)
{
    int32_t dj = GETJ(mat, j);

#if MKZ_DEBUG >= 1
    fprintf(stderr, "Decreasing col %d; was %d\n", j, mat->wt[dj]);
#endif
    mat->wt[dj] = decrS(mat->wt[dj]);
}

void
MkzRemoveJ(filter_matrix_t *mat, int32_t j)
{
    int32_t dj = GETJ(mat, j);

    if(mat->MKZA[dj] == MKZ_INF){
#if MKZ_DEBUG >= 1
	fprintf(stderr, "Prevented use of adr[%d]=-1 in MkzRemoveJ\n", j);
#endif
	return;
    }
#if MKZ_DEBUG >= 1
    fprintf(stderr, "Removing col %d of weight %d\n", j, mat->wt[dj]);
#endif
    mat->wt[dj] = 0;
    // remove j from the QA structure
    MkzDelete(mat->MKZQ, mat->MKZA, mat->MKZA[dj]);
    mat->MKZA[dj] = MKZ_INF;
#if MKZ_DEBUG >= 1
    MkzIsHeap(mat->MKZQ);
#endif
}

// let's say we remove some columns with the highest Mkz count.
// Not very pertinent right now.
int
MkzDeleteHeavyColumns(report_t *rep MAYBE_UNUSED, filter_matrix_t *mat MAYBE_UNUSED)
{
#if 1
    return 0;
#else
    int32_t j;
    int nj, njmax, w;

    if(MkzQueueCardinality(mat->MKZQ) < 5000)
	return 0;
    njmax = 10; // humf
    for(nj = 0; nj < njmax; nj++){
	j = MkzGet(mat->MKZQ, Q[0], 0);
	w = mat->wt[GETJ(mat, j)]; // make a copy of the weight
	MkzRemoveJ(mat, j);
        // mat->wt[j] was put to 0...
        mat->wt[GETJ(mat, j)] = -w; // restore and update
	Q[0]--;
    }
    return njmax;
#endif
}

#endif /* USE_MARKOWITZ */
