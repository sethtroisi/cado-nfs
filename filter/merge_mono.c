// TODO: use unified compact lists...!
// TODO: reintroduce mat->weight

/* If one wants to change the R data structure, please check the diff of
   revision 1568, which clearly identifies the places where R is used. */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "utils.h"

#include "merge_opts.h"
#include "sparse.h"
#include "dclist.h"
#include "filter_matrix.h"
#include "mst.h"
#include "report.h"

#ifndef USE_MARKOWITZ
# include "swar.h"
#else
# include "markowitz.h"
#endif
#include "merge_mono.h"

#define DEBUG 0

#ifdef USE_MARKOWITZ
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#if DEBUG >= 1
/* total number of row additions in merge */
static unsigned long row_additions = 0;
#endif

// not mallocing to speed up(?).
static int
findBestIndex(filter_matrix_t *mat, int m, int32_t *ind)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], i, j, imin, wmin, w;

    ASSERT(m <= MERGE_LEVEL_MAX);
    if(m == 2)
	return 0;
    fillRowAddMatrix(A, mat, m, ind);
    // iterate over all vertices
    imin = -1;
    wmin = 0;
    for(i = 0; i < m; i++){
	// compute the new total weight if i is used as pivot
	w = 0;
	for(j = 0; j < m; j++)
	    // note A[i][i] = 0
	    w += A[i][j];
#if DEBUG >= 1
	fprintf(stderr, "W[%d]=%d\n", i, w);
#endif
	if((imin == -1) || (w < wmin)){
	    imin = i;
	    wmin = w;
	}
    }
    return imin;
}

#if DEBUG >= 2
// We should have mat->wt[j] == m == nchk
static void
checkCoherence(filter_matrix_t *mat, int m, int j)
{
    int nchk = 0, k;

#ifndef USE_COMPACT_R
    for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
	if(mat->R[GETJ(mat, j)][k] != -1)
	    nchk++;
    ASSERT(nchk == (mat->wt[GETJ(mat, j)] >= 0 ? mat->wt[GETJ(mat, j)] : -mat->wt[GETJ(mat, j)]));
    if(m != -1){
	if(nchk != m){
	    fprintf(stderr, "HYPERCHECK:");
	    fprintf(stderr, "mat->R[%d][0]=%ld, m=%d\n", j,
                    (long int) mat->R[GETJ(mat, j)][0], m);
	    fprintf(stderr, "Gasp: nchk=%d\n", nchk);
	}
	ASSERT(nchk == m);
    }
#else
    fprintf(stderr, "R: NYI in checkCoherence\n");
    exit(1);
#endif
}
#endif

//////////////////////////////////////////////////////////////////////
// making things independent of the real data structure used
//////////////////////////////////////////////////////////////////////

static void
removeColDefinitely(report_t *rep, filter_matrix_t *mat, int32_t j)
{
    int32_t k;

#ifndef USE_COMPACT_R
    for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++)
	if(mat->R[GETJ(mat, j)][k] != -1){
# if TRACE_COL >= 0
	    if(j == TRACE_COL)
		fprintf(stderr, "deleteAllCols: row is %d\n",mat->R[GETJ(mat, j)][k]);
# endif
	    remove_j_from_row(mat, mat->R[GETJ(mat, j)][k], j);
	    removeRowDefinitely(rep, mat, mat->R[GETJ(mat, j)][k]);
	    mat->rem_ncols--;
	}
    mat->wt[GETJ(mat, j)] = 0;
#else
    fprintf(stderr, "R: NYI in removeColDefinitely\n");
    exit(1);
#endif
}

static void
removeColumnAndUpdate(filter_matrix_t *mat, int j)
{
#ifndef USE_MARKOWITZ
    remove_j_from_SWAR(mat, j);
#else
    MkzRemoveJ(mat, j);
#endif
    mat->wt[GETJ(mat, j)] = 0;
    freeRj(mat, j);
}

// The cell [i, j] may be incorporated to the data structure, at least
// if j is not too heavy, etc.
static void
addCellAndUpdate(filter_matrix_t *mat, int i, int32_t j)
{
#if TRACE_ROW >= 0
    if(i == TRACE_ROW)
	fprintf(stderr, "TRACE_ROW: addCellSWAR i=%d j=%d\n", i, j);
#endif
#ifndef USE_MARKOWITZ
    if(addColSWAR(mat, j) >= 0)
	// update R[j] by adding i
	add_i_to_Rj(mat, i, j);
#else
    {
    int w = MkzIncrCol(mat, j);
    if(w >= 0){
	if(w > mat->cwmax){
	    // first time...
	    removeColumnAndUpdate(mat, j);
	    mat->wt[j] = -w;
	}
	else{
	    // update R[j] by adding i
	    add_i_to_Rj(mat, i, j);
	    MkzUpdate(mat, i, j);
	}
    }
    }
#endif
}

// remove the cell (i,j), and updates matrix correspondingly.
static void
removeCellAndUpdate(filter_matrix_t *mat, int i, int32_t j)
{
#if TRACE_ROW >= 0
    if(i == TRACE_ROW){
	fprintf(stderr, "TRACE_ROW: removeCellAndUpdate i=%d j=%d\n", i, j);
    }
#endif
    if(mat->wt[GETJ(mat, j)] < 0){
	// if mat->wt[j] is already < 0, we don't care about
	// decreasing, updating, etc. except when > 2
	return;
    }
#ifndef USE_MARKOWITZ
    decreaseColWeightSWAR(mat, j);
#else
    MkzDecreaseColWeight(mat, j);
#endif
    // update R[j] by removing i
    remove_i_from_Rj(mat, i, j);
}

// These columns are simply removed from the current squeleton matrix, but
// not from the small matrix.
#ifndef USE_MARKOWITZ
static int
deleteAllColsFromStack(report_t *rep, filter_matrix_t *mat, int iS)
{
#ifdef USE_MARKOWITZ
    fprintf(stderr, "What's the point to have deleteAllColsFromStack?\n");
    return 0;
#else
    dclist dcl;
    int32_t j;
    int k, njrem = 0;

    while(1){
        dcl = dclistFirst (mat->S[iS]);
	if (dcl == NULL)
	    break;
	j = dcl->j;
	njrem++;
#if DEBUG >= 1
	fprintf(stderr, "Removing column %d from S[%d]\n", j, iS);
#endif
#if TRACE_COL >= 0
	if(j == TRACE_COL)
	    fprintf(stderr, "Removing column %d from S[%d]\n", j, iS);
#endif
	// we destroy column j
	// TODO: do faster?
	if(iS == 0)
	    mat->rem_ncols--;
	if(iS == 1)
	    removeColDefinitely(rep, mat, j);
	k = mat->wt[GETJ(mat, j)]; // make a copy of the weight
	removeColumnAndUpdate(mat, j);
	// mat->wt[j] was put to 0...
	mat->wt[GETJ(mat, j)] = -k; // restore and update
    }
#if DEBUG >= 1
    if(njrem > 0)
	fprintf(stderr, "deleteAllColsFromStack[%d]: %d\n", iS, njrem);
#endif
    return njrem;
#endif
}
#endif

//////////////////////////////////////////////////////////////////////
// now, these guys are generic...!

/* weight of buried columns for row i */
#if BURIED_MODEL == 0
#define BURIEDROWWEIGHT(mat,i) mat->wburied[i]
#else
#define BURIEDROWWEIGHT(mat,i) (int) (mat->bdensity * (double) mat->nburied + 0.5)
#endif

// for all j in row[i], removes j and update data
static void
removeRowAndUpdate(filter_matrix_t *mat, int i)
{
    int k;

#if TRACE_ROW >= 0
    if(i == TRACE_ROW)
	fprintf(stderr, "TRACE_ROW: removeRowAndUpdate i=%d\n", i);
#endif
    mat->weight -= lengthRow(mat, i) + BURIEDROWWEIGHT(mat, i);
    for(k = 1; k <= lengthRow(mat, i); k++){
#if TRACE_COL >= 0
	if(cell(mat, i, k) == TRACE_COL){
	    fprintf(stderr, "removeRowAndUpdate removes %d from R_%d\n",
		    TRACE_COL, i);
	}
#endif
	removeCellAndUpdate(mat, i, cell(mat, i, k));
    }
}

// All entries M[i, j] are potentially added to the structure.
static void
addOneRowAndUpdate(filter_matrix_t *mat, int i)
{
    int k;

    mat->weight += lengthRow(mat, i) + BURIEDROWWEIGHT(mat,i);
    for(k = 1; k <= lengthRow(mat, i); k++)
	addCellAndUpdate(mat, i, cell(mat, i, k));
}

// realize mat[i1] += mat[i2] and update the data structure.
// len could be the real length of row[i1]+row[i2] or -1.
void
addRowsAndUpdate(filter_matrix_t *mat, int i1, int i2, int len)
{
    // cleaner one, that shares addRowsData() to prepare the next move...!
    // i1 is to disappear, replaced by a new one
    removeRowAndUpdate(mat, i1);
    // we know the length of row[i1]+row[i2]
    addRows(mat->rows, i1, i2, len);
    addOneRowAndUpdate(mat, i1);
}

static int
removeSingletons(report_t *rep, filter_matrix_t *mat)
{
#ifndef USE_MARKOWITZ
    return deleteAllColsFromStack(rep, mat, 1);
#else
    int32_t j;
    int njrem = 0;

    for(j = mat->jmin; j < mat->jmax; j++)
	if(mat->wt[GETJ(mat, j)] == 1){
	    removeColDefinitely(rep, mat, j);
	    njrem++;
	}
    return njrem;
#endif
}

int
deleteHeavyColumns(report_t *rep, filter_matrix_t *mat)
{
#ifndef USE_MARKOWITZ
    return deleteAllColsFromStack(rep, mat, mat->cwmax+1);
#else
    return MkzDeleteHeavyColumns(rep, mat);
#endif
}

void
removeRowDefinitely(report_t *rep, filter_matrix_t *mat, int32_t i)
{
    removeRowAndUpdate(mat, i);
    destroyRow(mat, i);
    report1(rep, i);
    mat->rem_nrows--;
}

// try all combinations to find the smaller one; resists to m==1
static void
tryAllCombinations(report_t *rep, filter_matrix_t *mat, int m, int32_t *ind)
{
    int i, k;

    if(m == 1)
	fprintf(stderr, "Warning: m==1 in tryAllCombinations\n");
    i = findBestIndex(mat, m, ind);
#if DEBUG >= 1
    fprintf(stderr, "Minimal is i=%d (%d %d)\n", i, ind[0], ind[1]);
#endif
    for(k = 0; k < m; k++)
	if(k != i){
	    addRowsAndUpdate(mat, ind[k], ind[i], -1);
#if DEBUG >= 1
	    fprintf(stderr, "new row[%d]=", ind[k]);
	    print_row(mat, ind[k]);
	    fprintf(stderr, "\n");
#endif
	}
    if(i > 0){
	// put ind[i] in front
	int itmp = ind[i];
	for(k = i; k > 0; k--)
	    ind[k] = ind[k-1];
	ind[0] = itmp;
    }
#if DEBUG >= 1
    fprintf(stderr, "=> new_ind: %d %d\n", ind[0], ind[1]);
#endif
    reportn(rep, ind, m);
    removeRowAndUpdate(mat, ind[0]);
    destroyRow(mat, ind[0]);
}

// add u to its sons; we save the history in the history array, so that
// we can report all at once and prepare for MPI.
// A[i][j] contains the estimated weight/length of R[ind[i]]+R[ind[j]].
static int
addFatherToSonsRec(int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		   filter_matrix_t *mat, int m, int *ind,
		   int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX],
		   int *father,
		   int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		   int u, int level0)
{
    int k, itab = 1, i1, i2, level = level0, len;

    if(sons[u][0] == 0)
	// nothing to do for a leaf...!
	return -1;
    i2 = ind[u];
    if(u == father[u])
	history[level0][itab++] = i2; // trick for the root!!!
    else
	// the usual trick not to destroy row
	history[level0][itab++] = -(i2+1);
    for(k = 1; k <= sons[u][0]; k++){
	i1 = addFatherToSonsRec(history, mat, m, ind, A,
				father, sons, sons[u][k], level+1);
	if(i1 != -1)
	    level = i1;
	i1 = ind[sons[u][k]];
	// add u to its son
#if 1
	// recover true length of non-buried part for R[i1]+R[i2]
	len = A[sons[u][k]][u] - buriedRowsWeight (mat, i1, i2);
	addRowsAndUpdate(mat, i1, i2, len);
#else
	addRowsAndUpdate(mat, i1, i2, -1);
#endif
	history[level0][itab++] = i1;
    }
    history[level0][0] = itab-1;
    return level;
}

int
addFatherToSons(int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		filter_matrix_t *mat, int m, int *ind,
		int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX],
		int *father,
                int *height MAYBE_UNUSED, int hmax MAYBE_UNUSED,
		int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1])
{
    return addFatherToSonsRec(history, mat, m, ind, A, father, sons, 0, 0);
}

void
MSTWithA(report_t *rep, filter_matrix_t *mat, int m, int32_t *ind, double *tMST,
	 int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
    int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];
    int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];
    int father[MERGE_LEVEL_MAX], height[MERGE_LEVEL_MAX], hmax, i, w;

    *tMST = seconds();
    hmax = minimalSpanningTree(&w, father, height, sons, m, A);
    *tMST = seconds()-*tMST;
#if DEBUG >= 1
    printMST(father, sons, m);
#endif
    hmax = addFatherToSons(history, mat, m, ind, A, father, height, hmax,sons);
    for(i = hmax; i >= 0; i--)
#if 0
	reporthis(rep, history, i);
#else
        reportn(rep, history[i]+1, history[i][0]);
#endif
    removeRowAndUpdate(mat, ind[0]);
    destroyRow(mat, ind[0]);
}

static void
useMinimalSpanningTree(report_t *rep, filter_matrix_t *mat, int m,
		       int32_t *ind, double *tfill, double *tMST)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX];

    *tfill = seconds();
    fillRowAddMatrix(A, mat, m, ind);
    *tfill = seconds()-*tfill;
    MSTWithA(rep, mat, m, ind, tMST, A);
}

static void
findOptimalCombination(report_t *rep, filter_matrix_t *mat, int m,
		       int32_t *ind, double *tfill, double *tMST, int useMST)
{
    if((m <= 2) || !useMST){
	*tfill = *tMST = 0;
	tryAllCombinations(rep, mat, m, ind);
    }
    else
	useMinimalSpanningTree(rep, mat, m, ind, tfill, tMST);
}

#if DEBUG >= 1
static void
checkWeight(filter_matrix_t *mat, int32_t j)
{
    int i, w = 0;

    fprintf(stderr, "Rows containing %ld:", (long int) j);
    for(i = 0; i < mat->nrows; i++)
	if(!isRowNull(mat, i))
	    if(hasCol(mat->rows, i, j)){
		fprintf(stderr, " %d", i);
		w++;
	    }
    fprintf(stderr, "\n");
    ASSERT(w == (mat->wt[GETJ(mat, j)] >= 0 ? mat->wt[GETJ(mat, j)] : -mat->wt[GETJ(mat, j)]));
}

static void
checkMatrixWeight(filter_matrix_t *mat)
{
    unsigned long w = 0;
    int i;

    for(i = 0; i < mat->nrows; i++)
	if(!isRowNull(mat, i))
	    w += lengthRow(mat, i);
    ASSERT(w == mat->weight);
}
#endif

// j has weight m, which should coherent with mat->wt[j] == m
static void
mergeForColumn(report_t *rep, double *tt, double *tfill, double *tMST,
	       filter_matrix_t *mat, int m, int32_t j, int useMST)
{
    int32_t ind[MERGE_LEVEL_MAX];
    int ni, k;

    /* each m-merge leads to m-1 additions of rows */
#if DEBUG >= 1
    row_additions += m - 1;
#endif
#if BURIED_MODEL == 1
    /* Let c be the average density. When we add a row j to a row i:
       (a) either i has 0 and j has 1 with probability c*(1-c),
           the weight of row i increases by 1
       (b) both i and j have 1 with probablity c^2
           the weight of row i decreases by 1
       (c) in the remaining cases, the weight of row i remains unchanged
       Thus the number of non-zero elements in the column increases by
       c*(1-c) - c^2 = c - 2c^2, and the average density increases by
       (c - 2c^2) / nrows */
    for (k = 0; k < m - 1; k++)
      mat->bdensity += mat->bdensity * (1.0 - 2.0 * mat->bdensity)
        / (double) mat->nrows;
#endif

    if(m > MERGE_LEVEL_MAX){
	fprintf(stderr, "PB: m=%d > MERGE_LEVEL_MAX=%d\n", m, MERGE_LEVEL_MAX);
	exit(-1);
    }
#ifdef USE_MARKOWITZ
# if 0
    // let's be cautious...
    if(mat->wt[GETJ(mat, j)] != m){
	fprintf(stderr, "GASP: wt[%d]=%d != %d\n", j, mat->wt[j], m);
    }
# endif
#endif

#if DEBUG >= 1
    fprintf(stderr, "Treating column %d of weight %d",j,mat->wt[GETJ(mat, j)]);
    fprintf(stderr, "\n");
#endif
#if DEBUG >= 2
    fprintf(stderr, "Status before next j=%d to start\n", j);
    printSWAR(mat, mat->ncols);
#endif
    // the corresponding rows are in R[j], skipping 1st cell and -1's
#if DEBUG >= 2
    checkCoherence(mat, m, j);
#endif
    ni = 0;
#ifndef USE_COMPACT_R
    for(k = 1; k <= mat->R[GETJ(mat, j)][0]; k++){
	if(mat->R[GETJ(mat, j)][k] != -1){
	    ind[ni++] = mat->R[GETJ(mat, j)][k];
	    if(ni == m)
		break;
	}
    }
#else
    fprintf(stderr, "R: NYI in mergeForColumn\n");
    exit(1);
#endif
#if DEBUG >= 1
    fprintf(stderr, " %d", j);
    fprintf(stderr, " => the %d rows are:\n", m);
    for(k = 0; k < m; k++){
	fprintf(stderr, "row[%d]=", ind[k]);
	print_row(mat, ind[k]);
	fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
#endif
    *tt = seconds();
    findOptimalCombination(rep, mat, m, ind, tfill, tMST, useMST);
    *tt = seconds()-(*tt);
    mat->rem_nrows--;
    mat->rem_ncols--;
    removeColumnAndUpdate(mat, j);
}

#ifndef USE_MARKOWITZ
// maxdo is 0 if we want to perform a non-bounded number of operations;
// an integer >= 1 otherwise.
// Default for the monoproc version is 0.
static int
merge_m_fast(report_t *rep, filter_matrix_t *mat, int m, int maxdo, int verbose)
{
    double totopt=0, tot=seconds(), tt, totfill=0, totMST=0, tfill, tMST;
    double totdel = 0;
    int njproc = 0, njrem = 0;
    int32_t j;
    dclist dcl;
    int report = (verbose == 0) ? 10000 : 1000;

#if DEBUG >= 1
    fprintf(stderr, "Weight %d:", m);
#endif

    /* checkWeight(mat, 132461); */

    while(1){
#if DEBUG >= 1
	checkMatrixWeight(mat);
	checkData(mat);
#endif
	// get next j
	dcl = dclistFirst (mat->S[m]);
	if (dcl == NULL)
	    break;
	j = dcl->j;
	njproc++;
	if(!(njproc % report))
	  fprintf (stderr, "# %d columns of weight %d processed at %2.2lf\r",
		   njproc, m, seconds ());
	mergeForColumn(rep, &tt, &tfill, &tMST, mat, m, j, 1);
#if TRACE_COL >= 0
	fprintf(stderr, "wt[%d]=%d after findOptimalCombination\n",
		TRACE_COL, mat->wt[GETJ(mat, TRACE_COL)]);
#endif
	totopt += tt;
	totfill += tfill;
	totMST += tMST;
	tt = seconds();
	njrem += deleteHeavyColumns(rep, mat);
	totdel += (seconds()-tt);
#if TRACE_COL >= 0
        fprintf(stderr, "wt[%d]=%d after deleteHeavyColumns\n",
                TRACE_COL, mat->wt[GETJ(mat, TRACE_COL)]);
#endif
	maxdo--;
	if(maxdo == 0)
	    // what a trick: maxdo = 0 => never == 0 until infinity
	    return njrem;
    }
    fprintf (stderr, "\n"); /* to keep last output on screen */
    tot = seconds()-tot;
    fprintf(stderr, "m=%d nj=%d", m, njproc);
    fprintf(stderr, " findopt=%2.2lf (fill=%2.2lf mst=%2.2lf) tot=%2.2lf",
	    totopt, totfill, totMST, tot);
    fprintf(stderr, " del=%2.2lf\n", totdel);
	
#if DEBUG >= 1
    fprintf(stderr, "Status at the end of the process\n");
    Sparse2Tex(mat);
#endif
    return njrem;
}

static int
merge_m(report_t *rep, filter_matrix_t *mat, int m, int maxdo, int verbose)
{
  return merge_m_fast(rep, mat, m, maxdo, verbose);
}

static int
minColWeight(filter_matrix_t *mat)
{
    // in this case, locating minw is easy: inspect the stacks!
    int m;

    for (m = 1; m <= mat->cwmax; m++)
      if (dclistFirst(mat->S[m]) != NULL)
	return m;
    return -1;
}
#endif

#if !defined(USE_MPI)
static int
isRowToBeDeleted(filter_matrix_t *mat, int i)
{
    int heavy = 0;

    int k;

    heavy = 1;
    for(k = 1; k <= lengthRow(mat, i); k++)
	if((cell(mat,i,k) != -1) && (mat->wt[GETJ(mat, cell(mat,i,k))] >= 0)){
	    heavy = 0;
	    break;
	}
    return heavy;
}
#endif

static int
inspectRowWeight(report_t *rep, filter_matrix_t *mat)
{
    int i, nirem = 0, niremmax = 128;
#if !defined(USE_MPI)
    int useless = 0;
#endif

    for(i = 0; i < mat->nrows; i++){
	if((mat->rem_nrows - mat->rem_ncols) <= mat->keep)
	    return nirem;
	if(!isRowNull(mat, i)){
	    if(lengthRow(mat, i) > mat->rwmax){
#if DEBUG >= 1
		fprintf(stderr, "Removing too heavy row[%d]: %d\n",
			i, lengthRow(mat, i));
#endif
#if TRACE_ROW >= 0
		if(i == TRACE_ROW)
		    fprintf(stderr,
			    "TRACE_ROW: removing too heavy row[%d]: %d\n",
			    i, lengthRow(mat, i));
#endif
		removeRowDefinitely(rep, mat, i);
		nirem++;
		if(!(nirem % 10000))
		    fprintf(stderr, "#removed_rows=%d at %2.2lf\n",
			    nirem, seconds());
		if(nirem > niremmax)
		    break;
	    }
	}
    }
#if !defined(USE_MPI)
    /* only activated rarely...! */
    for(i = 0; i < mat->nrows; i++)
	if(!isRowNull(mat, i) && isRowToBeDeleted(mat, i)){
# if DEBUG >= 1
	    fprintf(stderr, "Row %d is no longer useful in merge\n", i);
# endif
	    removeRowAndUpdate(mat, i);
	    destroyRow(mat, i);
	    useless++;
	}
    fprintf(stderr, "#useless rows=%d\n", useless);
#endif
    return nirem;
}

// It could be fun to find rows s.t. at least one j is of small weight.
// Or components of rows with some sharing?
// We could define sf(row) = sum_{j in row} w(j) and remove the
// rows of smallest value of sf?
// For the time being, remove heaviest rows.
static int
deleteScore(filter_matrix_t *mat, int32_t i)
{
#if 1
    // plain weight to remove heaviest rows
    return lengthRow(mat, i);
#endif
#if 0
    // -plain weight to remove lightest rows
    return -lengthRow(mat, i);
#endif
#if 0
    // not using rows with too many heavy column part
    int k, s = 0;

    for(k = 1; k <= lengthRow(mat, i); k++)
	s += abs(mat->wt[GETJ(mat, mat->rows[i][k])]);
    return s;
#endif
#if 0
    // not using rows with too many light columns
    int k, s = 0;

    for(k = 1; k <= lengthRow(mat, i); k++)
	s += abs(mat->wt[GETJ(mat, mat->rows[i][k])]);
    return -s;
#endif
#if 0
    // return the weight of the lightest column
    int k, s = abs(mat->wt[GETJ(mat, mat->rows[i][1])]);

    for(k = 2; k <= lengthRow(mat, i); k++)
	if(abs(mat->wt[GETJ(mat, mat->rows[i][k])]) > s)
	    s = abs(mat->wt[GETJ(mat, mat->rows[i][k])]);
    return s;
#endif
}

// locate columns of weight 3 and delete one of the rows, but just one!
static int
findSuperfluousRowsFor2(int *tmp, int ntmp, filter_matrix_t *mat)
{
#ifdef USE_MARKOWITZ
# if DEBUG >= 1
    fprintf(stderr, "What's the use of findSuperfluousRowsFor2 for MKZ?\n");
# endif
    return 0;
#else
    dclist dcl;
    int32_t *Rj;
    int itmp = 0, j, k, kmin;

    for (dcl = dclistFirst (mat->S[3]); dcl != NULL; dcl = dcl->next){
	j = dcl->j;
#ifndef USE_COMPACT_R
	Rj = mat->R[GETJ(mat, j)];
	// be semi-stupid
	for(kmin = 1; Rj[kmin] == -1; kmin++); // should stop at some point...!
	for(k = kmin+1; k <= mat->R[GETJ(mat, j)][0]; k++)
	    if(Rj[k] != -1)
		if(lengthRow(mat, Rj[k]) > lengthRow(mat, Rj[kmin]))
		    kmin = k;
	tmp[itmp++] = 0; // horrible trick not to change things elsewhere...!
	tmp[itmp++] = Rj[kmin];
	if(itmp >= ntmp)
	    // should not happen, but...
	    break;
#else
	fprintf(stderr, "R: NYI in findSuperfluousRowsFor2\n");
	exit(1);
#endif
    }
    return itmp;
#endif
}

// this guy is linear => total is quadratic, be careful!
static int
findSuperfluousRows(int *tmp, int ntmp, filter_matrix_t *mat, int m)
{
    int i, itmp;

    if(m == 2){
	itmp = findSuperfluousRowsFor2(tmp, ntmp, mat);
	if(itmp != 0)
	    return itmp;
    }
    for(i = 0, itmp = 0; i < mat->nrows; i++)
        if(!isRowNull(mat, i)){
	    tmp[itmp++] = deleteScore(mat, i);
	    tmp[itmp++] = i;
	}
    // rows with largest score will be at the end
    qsort(tmp, itmp>>1, 2 * sizeof(int), cmp);
    return itmp;
}

// Delete superfluous rows s.t. nrows-ncols >= keep.
// Let's say we are at "step" m.
static int
deleteSuperfluousRows(report_t *rep, filter_matrix_t *mat, int keep, int niremmax, int m)
{
    int nirem = 0, *tmp, ntmp, i;

    if((mat->rem_nrows - mat->rem_ncols) <= keep)
	return 0;
    ntmp = mat->rem_nrows << 1;
    tmp = (int *)malloc(ntmp * sizeof(int));
    ntmp = findSuperfluousRows(tmp, ntmp, mat, m);
#if 0
    fprintf(stderr, "Number of removable rows[%d]: %d\n", m, ntmp>>1);
#endif
    // remove rows with largest score
    for(i = ntmp-1; i >= 0; i -= 2){
	if((nirem >= niremmax) || (mat->rem_nrows - mat->rem_ncols) <= keep)
	    break;
	if(mat->rows[tmp[i]] != NULL){
	    removeRowDefinitely(rep, mat, tmp[i]);
	    nirem++;
	}
    }
    free(tmp);
    return nirem;
}

#ifndef USE_MARKOWITZ
void
merge(report_t *rep, filter_matrix_t *mat, int maxlevel, int verbose, int forbw)
{
    unsigned long bwcostmin = 0, oldcost = 0, cost;
    int old_nrows, old_ncols, m, mm, njrem = 0, ncost = 0, ncostmax;

    ncostmax = 20; // was 5
    m = 2;
    while(1){
	cost = ((unsigned long)mat->rem_ncols) * ((unsigned long)mat->weight);
	fprintf(stderr, "w(M)=%lu, w(M)*ncols=%lu", mat->weight, cost);
	fprintf(stderr, " w(M)/ncols=%2.2lf\n",
		((double)mat->weight)/((double)mat->rem_ncols));
	if((bwcostmin == 0) || (cost < bwcostmin))
	    bwcostmin = cost;
	if(forbw)
	    // what a trick!!!!
	    fprintf(rep->outfile, "BWCOST: %lu\n", cost);
	if(forbw && (oldcost != 0) && (cost > oldcost)){
	    ncost++;
	    fprintf(stderr, "WARNING: New cost > old cost (%2.6e) [%d/%d]\n",
		    ((double)cost)-((double)oldcost), ncost, ncostmax);
	    if(ncost >= ncostmax){
		int nirem;

		fprintf(stderr, "WARNING: New cost > old cost %d times",
                        ncost);
		fprintf(stderr, " in a row:");
		nirem = deleteSuperfluousRows(rep, mat, mat->keep, 128, m);
		if(nirem == 0){
		    fprintf(stderr, " stopping\n");
		    break;
		}
		else{
		    fprintf(stderr, " try again after removing %d rows!\n",
			    nirem);
		    continue;
		}
	    }
	}
	else
	    ncost = 0;
	oldcost = cost;
	fprintf(stderr, "Performing merge %d at %2.2lf\n", m, seconds());
	old_nrows = mat->rem_nrows;
	old_ncols = mat->rem_ncols;
	if(m == 1)
	    njrem += removeSingletons(rep, mat);
	else
            njrem += merge_m (rep, mat, m, 0, verbose);
#if DEBUG >= 1
	checkData(mat);
#endif
#if 0
    {
	int i, ni = 0, j0 = 191342;

	for(i = 0; i < mat->nrows; i++)
	    if((mat->data[i].val != NULL) && hasCol(mat->rows, i, j0))
		ni++;
	printf("CHECK: %d %d\n", ni, mat->wt[GETJ(mat, j0)]);
    }
#endif
	fprintf(stderr, "=> nrows=%d ncols=%d (%d) njrem=%d\n",
		mat->rem_nrows, mat->rem_ncols,
		mat->rem_nrows - mat->rem_ncols, njrem);
	inspectRowWeight(rep, mat);
	deleteEmptyColumns(mat);
	mm = minColWeight(mat);
	fprintf(stderr, "Min col weight = %d\n", mm);
#if M_STRATEGY == 2
	// jump to the next minimal merge level immediately
	m = mm;
	if((m > maxlevel) || (m <= 0))
	    break;
#else
#  if M_STRATEGY == 1
	if(mm < m)
	    // something new happened, anyway
	    m = mm;
	else
#  endif
	    if((old_nrows == mat->rem_nrows) && (old_ncols == mat->rem_ncols)){
		// nothing happened this time and mm > m
		m = mm;
		if((m > maxlevel) || (m <= 0))
		    break;
	    }
#endif
    }
    if(forbw){
	fprintf(rep->outfile, "BWCOSTMIN: %lu\n", bwcostmin);
	fprintf(stderr, "Minimal bwcost found: %lu\n", bwcostmin);
    }
}
#endif

static uint64_t
my_cost(unsigned long N, unsigned long w, int forbw)
{
    if(forbw == 2){
	double K1 = .19e-9, K2 = 3.4e-05, K3 = 1.4e-10; // kinda average
	double dN = (double)N, dw = (double) w;

	return (uint64_t)((K1+K3)*dN*dw+K2*dN*log(dN)*log(dN));
    }
    else if(forbw == 3)
	return (uint64_t)(w/N);
    else if(forbw <= 1)
	return ((uint64_t)N)*((uint64_t)w);
    return (uint64_t)0;
}

static void
mergeForColumn2(report_t *rep, filter_matrix_t *mat, int *njrem,
		double *totopt, double *totfill, double *totMST,
		double *totdel, int useMST, int32_t j)
{
    double tt, tfill, tMST;

    mergeForColumn(rep, &tt, &tfill, &tMST, mat, mat->wt[j], j, useMST);
    *totopt += tt;
    *totfill += tfill;
    *totMST += tMST;
    tt = seconds();
    *njrem += deleteHeavyColumns(rep, mat);
    *totdel += (seconds()-tt);
}

// njrem is the number of columns removed.
void
doOneMerge(report_t *rep, filter_matrix_t *mat, int *njrem, double *totopt, double *totfill, double *totMST, double *totdel, int m, int maxdo, int useMST, int verbose)
{
#ifdef USE_MARKOWITZ
    fprintf(stderr, "What's the use of doOneMerge for MKZ?\n");
#else
    dclist dcl;
    int32_t j;

    if(m == 1){
	// do all of these
	*njrem += removeSingletons(rep, mat);
	fprintf(stderr, "T=%d after m=1\n", (int)(seconds()));
    }
    else if (m == 2)
      *njrem += merge_m (rep, mat, m, maxdo, verbose);
    else{
	//	    fprintf(stderr, "Performing one merge for m=%d\n", m);
        dcl = dclistFirst (mat->S[m]);
	j = dcl->j;
	mergeForColumn2(rep, mat, njrem, totopt, totfill, totMST, totdel,
			useMST, j);
    }
#endif
}

int
number_of_superfluous_rows(filter_matrix_t *mat)
{
    int kappa = (mat->rem_nrows-mat->rem_ncols) / mat->keep, ni2rem;

    if(kappa <= (1<<4))
	ni2rem = mat->keep/2;
    else if(kappa <= (1<<5))
	ni2rem = mat->keep * (kappa/4);
    else if(kappa <= (1<<10))
	ni2rem = mat->keep * (kappa/8);
    else if(kappa <= (1<<15))
	ni2rem = mat->keep * (kappa/16);
    else if(kappa <= (1<<20))
	ni2rem = mat->keep * (kappa/32);
    else
	ni2rem = mat->keep * (kappa/64);
    return ni2rem;
}

int
deleteEmptyColumns(filter_matrix_t *mat)
{
#ifndef USE_MARKOWITZ
    return deleteEmptyColumnsSWAR(mat);
#else
    // nothing to do; the pb has to be solved elsewhere...!
    return 0;
#endif
}

void
mergeOneByOne(report_t *rep, filter_matrix_t *mat, int maxlevel, int verbose, int forbw, double ratio, int coverNmax)
{
    double totopt = 0.0, totfill = 0.0, totMST = 0.0, totdel = 0.0;
    uint64_t bwcostmin = 0, oldbwcost = 0, bwcost = 0;
    int old_nrows, old_ncols, m = 2, njrem = 0, ncost = 0, ncostmax, njproc;
    int target = 10000, ni2rem;
#ifndef USE_MARKOWITZ
    int mmax = 0;
#else
    int32_t dj, j, mkz;

    // clean things
    njrem = removeSingletons(rep, mat);
#endif

    fprintf(stderr, "Using mergeOneByOne\n");
    ncostmax = 20; // was 5
    njproc = 0;
    while(1){
	if(mat->itermax && (njproc >= mat->itermax)){
	    fprintf(stderr, "itermax=%d reached, stopping!\n", mat->itermax);
	    break;
	}
	oldbwcost = bwcost;
	old_nrows = mat->rem_nrows;
	old_ncols = mat->rem_ncols;
#ifndef USE_MARKOWITZ
	m = minColWeight(mat);
	if(m > mmax)
	    mmax = m;
	if(mmax > maxlevel){
	    fprintf(stderr, "maxlevel reached, stopping!\n");
	    break;
	}
	if(m == -1){
	    fprintf(stderr, "All stacks empty, stopping!\n");
	    break;
	}
	doOneMerge(rep, mat, &njrem, &totopt, &totfill, &totMST, &totdel,
		   m, 0, 1, verbose);
#else
	if(MkzQueueCardinality(mat->MKZQ) <= 1){
	    // rare event, I must say!
	    fprintf(stderr, "Q is almost empty: rare!!!\n");
	    break;
	}
	if((mat->mkzrnd == 0) || (mat->mkzrnd > mat->MKZQ[0]))
	    MkzPopQueue(&dj, &mkz, mat->MKZQ, mat->MKZA);
	else
	    MkzRemove(&dj, &mkz, mat->MKZQ, mat->MKZA, mat->mkzrnd);
	j = dj + mat->jmin;
# if DEBUG >= 1
	fprintf(stderr, "I popped j=%d wt=%d mkz=%d (#Q=%d)",
		j, mat->wt[dj], mkz, mat->MKZQ[0]);
	fprintf(stderr, " nrows=%d ncols=%d\n",mat->rem_nrows,mat->rem_ncols);
# endif
	if(mat->wt[dj] == 1){
# if DEBUG >= 1
	    fprintf(stderr,"Popped j=%d with w=%d\n", j, mat->wt[dj]);
#endif
	    removeColDefinitely(rep, mat, j);
	}
	else if(mat->wt[dj] > 0){
	    m = mat->wt[dj]; // for deleteSuperfluousRows below
	    mergeForColumn2(rep, mat, &njrem,
			    &totopt, &totfill, &totMST, &totdel, 1, j);
	}
# if 0
	else{
	    int w = mat->wt[dj];
	    removeColumnAndUpdate(mat, j);
	    mat->wt[dj] = w;
	}
# endif
#endif
	// number of columns removed
	njproc += old_ncols - mat->rem_ncols;
	deleteEmptyColumns(mat);
#ifndef USE_MARKOWITZ
	if((old_nrows == mat->rem_nrows) && (old_ncols == mat->rem_ncols)){
	    // is this supposed to be activated at some point? Really?
	    if((m > maxlevel) || (m <= 0)){
		fprintf(stderr, "Stopping, since m=%d // maxlevel=%d\n",
			m, maxlevel);
		break;
	    }
	}
#endif
	bwcost = my_cost ((unsigned long) mat->rem_nrows,
                          (unsigned long) mat->weight, forbw);
	if(njproc >= target){ // somewhat arbitrary...!
#ifdef USE_MARKOWITZ
	    njrem = removeSingletons(rep, mat);
#endif
	    ni2rem = number_of_superfluous_rows(mat);
	    deleteSuperfluousRows(rep, mat, mat->keep, ni2rem, m);
	    inspectRowWeight(rep, mat);
#ifndef USE_MARKOWITZ
	    fprintf(stderr, "mmax=%d ", mmax);
#endif	
	    fprintf(stderr, "N=%d (%d) w=%lu",
		    mat->rem_nrows, mat->rem_nrows - mat->rem_ncols,
		    mat->weight);
	    if(forbw == 2)
		fprintf(stderr, " bw=%"PRIu64"", bwcost);
	    else if(forbw == 3)
		fprintf(stderr, " w*N=%"PRIu64"",
			((uint64_t)mat->rem_nrows)
			*((uint64_t)mat->weight));
	    else if(forbw <= 1)
		fprintf(stderr, " w*N=%"PRIu64"", bwcost);
	    fprintf(stderr, " w/N=%2.2lf\n",
		    ((double)mat->weight)/((double)mat->rem_nrows));
	    // njrem=%d at %2.2lf\n",
	    if((forbw != 0) && (forbw != 3))
		// what a trick!!!!
		fprintf(rep->outfile, "BWCOST: %"PRIu64"\n", bwcost);
	    target = njproc + 10000;
	}
	if((bwcostmin == 0) || (bwcost < bwcostmin)){
	    bwcostmin = bwcost;
	    if((forbw != 0) && (forbw != 3))
		// what a trick!!!!
		fprintf(rep->outfile, "BWCOST: %"PRIu64"\n", bwcost);
	}
	// to be cleaned one day...
	if((forbw == 0) || (forbw == 2)){
	    double r = ((double)bwcost)/((double)bwcostmin);
	    if(r > ratio){
		if(mat->rem_nrows-mat->rem_ncols > mat->keep){
		    // drop all remaining columns at once
		    ni2rem = mat->rem_nrows-mat->rem_ncols+mat->keep;
		    fprintf(stderr, "Dropping %d rows at once\n", ni2rem);
		    deleteSuperfluousRows(rep, mat, mat->keep, ni2rem, -1);
		}
		else{
		    if(forbw == 0)
			fprintf(stderr, "cN too high, stopping [%2.2lf]\n", r);
		    else
			fprintf(stderr, "bw too high, stopping [%2.2lf]\n", r);
		    break;
		}
	    }
	}
	else if(forbw == 3)
          {
            if (bwcost >= (uint64_t) coverNmax)
              {
		fprintf (stderr, "w/N too high (%"PRIu64"), stopping\n",
                         bwcost);
		break;
              }
          }
	if((forbw != 0) && (oldbwcost != 0) && (bwcost > oldbwcost)){
	    ncost++;
#if 0
	    fprintf(stderr, "New cost > old cost (%2.6e) [%d/%d]\n",
		    ((double)bwcost)-((double)oldbwcost), ncost, ncostmax);
#endif
	    if(ncost >= ncostmax){
		int nirem;

		fprintf(stderr, "New cost > old cost %d times", ncost);
		fprintf(stderr, " in a row:");
		nirem = deleteSuperfluousRows(rep, mat, mat->keep, 128, m);
		if(nirem == 0){
		    fprintf(stderr, " stopping\n");
		    break;
		}
		else{
		    fprintf(stderr, " try again after removing %d rows!\n",
			    nirem);
		    njproc += nirem; // humf: odd name for njproc...!
		    continue;
		}
	    }
	}
	else
	    ncost = 0;
    }
    if(mat->itermax == 0)
	deleteSuperfluousRows(rep, mat, mat->keep,
			      mat->rem_nrows-mat->rem_ncols+mat->keep, -1);
    if((forbw != 0) && (forbw != 3)){
	fprintf(rep->outfile, "BWCOSTMIN: %"PRIu64"\n", bwcostmin);
	fprintf(stderr, "Minimal bwcost found: %"PRIu64"\n", bwcostmin);
    }
#if DEBUG >= 1
    fprintf(stderr, "Total number of row additions: %lu\n", row_additions);
#endif
}

//////////////////////////////////////////////////////////////////////
//
// Resume section: very much inspired by replay.c, of course...!
//
//////////////////////////////////////////////////////////////////////

// A line is "i i1 ... ik".
static int
indicesFromString(int32_t *ind, char *str)
{
    int ni = 0;
    char *tok = strtok(str, " ");

    while(tok != NULL){
	ind[ni++] = atoi(tok);
	tok = strtok(NULL, " ");
    }
    return ni;
}

// A line is "i i1 ... ik".
// If i >= 0 then
//     row[i] is to be added to rows i1...ik and destroyed at the end of
//     the process.
//     Works also is i is alone (hence: destroyed row).
// If i < 0 then
//     row[-i-1] is to be added to rows i1...ik and NOT destroyed.
//
static void
doAllAdds(report_t *rep, filter_matrix_t *mat, char *str)
{
    int32_t ind[MERGE_LEVEL_MAX], i0;
    int ni, sg, k;

    ni = indicesFromString(ind, str);
    if(ind[0] < 0){
	i0 = -ind[0]-1;
	sg = -1;
    }
    else{
	i0 = ind[0];
	sg = 1;
    }
    for(k = 1; k < ni; k++)
	addRowsAndUpdate(mat, ind[k], i0, -1);
    reportn(rep, ind, ni);
    if(sg > 0){
	// when ni == 1, then the corresponding row was too heavy and dropped
	// unless m = 1 was used, in which case S[1] will contain j and
	// can dispensed of next time (?) In the case of S[0], it'll be done
	// later on also (?)
	removeRowAndUpdate(mat, i0);
	destroyRow(mat, i0);
	mat->rem_nrows--;
	// the number of active j's is recomputed, anyway, later on
    }
}

// resumename is a file of the type mergehis.
// TODO: Compiles, but not really tested with Markowitz...!
void
resume(report_t *rep, filter_matrix_t *mat, char *resumename)
{
    FILE *resumefile = fopen(resumename, "r");
    char str[RELATION_MAX_BYTES];
    unsigned long addread = 0;
    int nactivej;
    int32_t j;
    char * rp;

    fprintf(stderr, "Resuming computations from %s\n", resumename);
    // skip first line containing nrows ncols
    rp = fgets(str, RELATION_MAX_BYTES, resumefile);
    ASSERT_ALWAYS(rp);

    fprintf(stderr, "Reading row additions\n");
    while(fgets(str, RELATION_MAX_BYTES, resumefile)){
	addread++;
	if((addread % 100000) == 0)
	    fprintf(stderr, "%lu lines read at %2.2lf\n", addread, seconds());
	if(str[strlen(str)-1] != '\n'){
	    fprintf(stderr, "Gasp: not a complete a line!");
	    fprintf(stderr, " I stop reading and go to the next phase\n");
	    break;
	}
	if(strncmp(str, "BWCOST", 6) != 0)
	    doAllAdds(rep, mat, str);
    }
    fclose(resumefile);
    for(j = mat->jmin; j < mat->jmax; j++)
#ifndef USE_MARKOWITZ
	if((mat->wt[GETJ(mat, j)] == 0) && (mat->A[GETJ(mat, j)] != NULL))
#else
	if((mat->wt[GETJ(mat,j)] == 0) && MkzIsAlive(mat->MKZA, GETJ(mat,j)))
#endif
	    // be sure j was removed...
	    removeColumnAndUpdate(mat, j);
    nactivej = 0;
    for(j = mat->jmin; j < mat->jmax; j++)
	if(mat->wt[GETJ(mat, j)] != 0)
	    nactivej++;
    mat->rem_ncols = nactivej;
    fprintf(stderr, "At the end of resume, we have");
    fprintf(stderr, " nrows=%d ncols=%d (%d)\n",
	    mat->rem_nrows, mat->rem_ncols, mat->rem_nrows-mat->rem_ncols);
}
