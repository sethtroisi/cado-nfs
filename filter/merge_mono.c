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

#include "portability.h"
#include "utils.h"

#include "merge_opts.h"
#include "filter_matrix.h"
#include "sparse.h"
#include "mst.h"
#include "report.h"

#include "markowitz.h"
#include "merge_mono.h"

#define DEBUG 0

#pragma GCC diagnostic ignored "-Wunused-parameter"

#if DEBUG >= 1
/* total number of row additions in merge */
static unsigned long row_additions = 0;
#endif

// not mallocing to speed up(?).
static int
findBestIndex(filter_matrix_t *mat, int m, int32_t *ind, int32_t ideal)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], i, j, imin, wmin, w;

    ASSERT(m <= MERGE_LEVEL_MAX);
    if(m == 2)
	return 0;
    fillRowAddMatrix(A, mat, m, ind, ideal);
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
	printf ("W[%d]=%d\n", i, w);
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

    for(k = 1; k <= mat->R[j][0]; k++)
	if(mat->R[j][k] != -1)
	    nchk++;
    ASSERT(nchk == (mat->wt[j] >= 0 ? mat->wt[j] : -mat->wt[j]));
    if(m != -1){
	if(nchk != m){
	    printf ("HYPERCHECK: mat->R[%d][0]=%ld, m=%d\n", j,
                     (long int) mat->R[j][0], m);
	    printf ("Gasp: nchk=%d\n", nchk);
	}
	ASSERT(nchk == m);
    }
}
#endif

//////////////////////////////////////////////////////////////////////
// making things independent of the real data structure used
//////////////////////////////////////////////////////////////////////

static void
removeColDefinitely(report_t *rep, filter_matrix_t *mat, int32_t j)
{
  unsigned int k;

    for(k = 1; k <= mat->R[j][0]; k++)
	if(mat->R[j][k] != UMAX(index_t)){
# if TRACE_COL >= 0
	    if(j == TRACE_COL)
		printf ("deleteAllCols: row is %d\n",mat->R[j][k]);
# endif
	    remove_j_from_row(mat, mat->R[j][k], j);
	    removeRowDefinitely(rep, mat, mat->R[j][k]);
	    mat->rem_ncols--;
	}
    mat->wt[j] = 0;
}

/* remove column j and update matrix */
static void
removeColumnAndUpdate(filter_matrix_t *mat, int j)
{
    MkzRemoveJ (mat, j);
    freeRj(mat, j);
}

// The cell [i, j] may be incorporated to the data structure, at least
// if j is not too heavy, etc.
static void
addCellAndUpdate(filter_matrix_t *mat, int i, int32_t j)
{
    int w = MkzIncrCol(mat, j);

    if (w >= 0)
      {
	if(w > mat->cwmax)
          {
	    // the weight of column w exceeds cwmax, thus we remove it
	    removeColumnAndUpdate(mat, j);
	    mat->wt[j] = -w;
          }
	else
          {
	    // update R[j] by adding i
	    add_i_to_Rj(mat, i, j);
	    MkzUpdate(mat, i, j);
          }
      }
}

// remove the cell (i,j), and updates matrix correspondingly.
// if final, also update the Markowitz counts
static void
removeCellAndUpdate(filter_matrix_t *mat, int i, int32_t j, int final)
{
#if TRACE_ROW >= 0
    if(i == TRACE_ROW){
	printf ("TRACE_ROW: removeCellAndUpdate i=%d j=%d\n", i, j);
    }
#endif
    if(mat->wt[j] < 0){
	// if mat->wt[j] is already < 0, we don't care about
	// decreasing, updating, etc. except when > 2
	return;
    }
    MkzDecreaseColWeight(mat, j);
    // update R[j] by removing i
    remove_i_from_Rj(mat, i, j);
    if (final)
      MkzUpdateDown (mat, i, j);
}

//////////////////////////////////////////////////////////////////////
// now, these guys are generic...!

// for all j in row[i], removes j and update data
// if final, also update the Markowitz counts
static void
removeRowAndUpdate(filter_matrix_t *mat, int i, int final)
{
  unsigned int k;

#if TRACE_ROW >= 0
    if(i == TRACE_ROW)
	printf ("TRACE_ROW: removeRowAndUpdate i=%d\n", i);
#endif
    mat->weight -= matLengthRow(mat, i);
    for(k = 1; k <= matLengthRow(mat, i); k++){
#if TRACE_COL >= 0
	if(matCell(mat, i, k) == TRACE_COL){
	    printf ("removeRowAndUpdate removes %d from R_%d\n", TRACE_COL, i);
	}
#endif
	removeCellAndUpdate(mat, i, matCell(mat, i, k), final);
    }
}

// All entries M[i, j] are potentially added to the structure.
static void
addOneRowAndUpdate(filter_matrix_t *mat, int i)
{
  unsigned int k;

    mat->weight += matLengthRow(mat, i);
    for(k = 1; k <= matLengthRow(mat, i); k++)
	addCellAndUpdate(mat, i, matCell(mat, i, k));
}

// realize mat[i1] += mat[i2] and update the data structure.
// len could be the real length of row[i1]+row[i2] or -1.
void
addRowsAndUpdate(filter_matrix_t *mat, int i1, int i2, int32_t j)
{
    // cleaner one, that shares addRowsData() to prepare the next move...!
    // i1 is to disappear, replaced by a new one
    removeRowAndUpdate(mat, i1, 0);
    // we know the length of row[i1]+row[i2]
    addRows(mat->rows, i1, i2, j);
    addOneRowAndUpdate(mat, i1);
}

static int
removeSingletons(report_t *rep, filter_matrix_t *mat)
{
  index_t j;
    int njrem = 0;

    for(j = 0; j < mat->ncols; j++)
	if(mat->wt[j] == 1){
	    removeColDefinitely(rep, mat, j);
	    njrem++;
	}
    return njrem;
}

int
deleteHeavyColumns(report_t *rep, filter_matrix_t *mat)
{
    return MkzDeleteHeavyColumns(rep, mat);
}

void
removeRowDefinitely(report_t *rep, filter_matrix_t *mat, int32_t i)
{
    removeRowAndUpdate(mat, i, 1);
    destroyRow(mat, i);
    report1(rep, i, -1);
    mat->rem_nrows--;
}

// try all combinations to find the smaller one; resists to m==1
static void
tryAllCombinations(report_t *rep, filter_matrix_t *mat, int m, int32_t *ind,
                   int32_t j)
{
    int i, k;

    if(m == 1)
      printf ("Warning: m=1 in tryAllCombinations\n");
    i = findBestIndex(mat, m, ind, j);
#if DEBUG >= 1
    printf ("Minimal is i=%d (%d %d)\n", i, ind[0], ind[1]);
#endif
    for(k = 0; k < m; k++)
	if(k != i){
	    addRowsAndUpdate(mat, ind[k], ind[i], j);
#if DEBUG >= 1
	    printf ("new row[%d]=", ind[k]);
	    print_row (mat, ind[k]);
	    printf ("\n");
#endif
	}
    if(i > 0){
	// put ind[i] in front
        // FIXME: can we simply swap ind[0] and ind[i]?
	int itmp = ind[i];
	for(k = i; k > 0; k--)
	    ind[k] = ind[k-1];
	ind[0] = itmp;
    }
#if DEBUG >= 1
    printf ("=> new_ind: %d %d\n", ind[0], ind[1]);
#endif
    reportn(rep, ind, m, j);
    removeRowAndUpdate(mat, ind[0], 1);
    destroyRow(mat, ind[0]);
}

// add u to its sons; we save the history in the history array, so that
// we can report all at once
// A[i][j] contains the estimated weight/length of R[ind[i]]+R[ind[j]].
static int
addFatherToSonsRec(int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		   filter_matrix_t *mat, int m, int *ind, int32_t j,
		   int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX],
		   int *father,
		   int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		   int u, int level0)
{
    int k, itab = 1, i1, i2, level = level0;

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
	i1 = addFatherToSonsRec(history, mat, m, ind, j, A,
				father, sons, sons[u][k], level+1);
	if(i1 != -1)
	    level = i1;
	i1 = ind[sons[u][k]];
	// add u to its son
	addRowsAndUpdate(mat, i1, i2, j);
	history[level0][itab++] = i1;
    }
    history[level0][0] = itab-1;
    return level;
}

int
addFatherToSons(int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		filter_matrix_t *mat, int m, int *ind, int32_t j,
		int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX],
		int *father,
                int *height MAYBE_UNUSED, int hmax MAYBE_UNUSED,
		int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1])
{
    return addFatherToSonsRec(history, mat, m, ind, j, A, father, sons, 0, 0);
}

void
MSTWithA(report_t *rep, filter_matrix_t *mat, int m, int32_t *ind, int32_t j, 
         double *tMST, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
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
    hmax = addFatherToSons(history, mat, m, ind, j,A, father, height, hmax,sons);
    for(i = hmax; i >= 0; i--)
#if 0
	reporthis(rep, history, i);
#else
        reportn(rep, history[i]+1, history[i][0], j);
#endif
    removeRowAndUpdate(mat, ind[0], 1);
    destroyRow(mat, ind[0]);
}

static void
useMinimalSpanningTree(report_t *rep, filter_matrix_t *mat, int m, int32_t *ind,
                       int32_t j, double *tfill, double *tMST)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX];

    *tfill = seconds();
    fillRowAddMatrix(A, mat, m, ind, j);
    *tfill = seconds()-*tfill;
    MSTWithA(rep, mat, m, ind, j, tMST, A);
}

static void
findOptimalCombination(report_t *rep, filter_matrix_t *mat, int m, int32_t *ind,
                       int32_t j, double *tfill, double *tMST, int useMST)
{
  if ((m <= 2) || (useMST == 0))
    {
      *tfill = *tMST = 0;
      tryAllCombinations (rep, mat, m, ind, j);
    }
  else
    useMinimalSpanningTree (rep, mat, m, ind, j, tfill, tMST);
}

#if DEBUG >= 1
static void
checkWeights (filter_matrix_t *mat)
{
  int *W, j, i, k, minw = INT_MAX;

  W = (int*) malloc (mat->ncols * sizeof(int));
  for (j = 0; j < mat->ncols; j++)
    W[j] = 0;
  for(i = 0; i < mat->nrows; i++)
    if(!isRowNull(mat, i))
      {
        for(k = 1; k <= lengthRow(mat, i); k++)
          {
            j = cell(mat, i, k);
            ASSERT_ALWAYS(0 <= j && j < mat->ncols);
            W[j] ++;
          }
      }
  for (j = 0; j < mat->ncols; j++)
    {
      static int count = 0;
      if (W[j] != mat->wt[j])
        {
          if (mat->wt[j] >= 0)
            {
              fprintf (stderr, "Wrong weight for column %d: expected %d, got %d\n",
                       j, mat->wt[j], W[j]);
              exit (1);
            }
          else if (W[j] <= mat->mergelevelmax)
            {
              if (count++ <= 10)
                fprintf (stderr, "Warning, column %d: wt=%d W=%d\n",
                         j, mat->wt[j], W[j]);
            }
        }
      if (W[j] == 1 && mat->wt[j] == 1)
        {
          fprintf (stderr, "Error, weight=1 for column %d\n", j);
          exit (1);
        }
      if (W[j] != 0 && mat->wt[j] >= 0 && W[j] < minw)
        minw = W[j];
      if (W[j] == 2 && mat->wt[j] == 2)
        {
          static int count = 0;
          if (count ++ <= 10)
            fprintf (stderr, "Warning, column %d: wt=%d W=%d\n",
                     j, mat->wt[j], W[j]);
        }
    }
  printf ("Minimum non-zero column weight: %d\n", minw);
  free (W);
}
#endif

#if 0
static void
checkWeight(filter_matrix_t *mat, int32_t j)
{
    int i, w = 0;

    printf ("Rows containing %ld:", (long int) j);
    for(i = 0; i < mat->nrows; i++)
	if(!isRowNull(mat, i))
	    if(hasCol(mat->rows, i, j)){
		printf (" %d", i);
		w++;
	    }
    printf ("\n");
    ASSERT(w == (mat->wt[j] >= 0 ? mat->wt[j] : -mat->wt[j]));
}
#endif

// j has weight m, which should coherent with mat->wt[j] == m
static void
mergeForColumn (report_t *rep, double *tt, double *tfill, double *tMST,
                filter_matrix_t *mat, int m, int32_t j, int useMST)
{
    int32_t ind[MERGE_LEVEL_MAX];
    unsigned int k;
    int ni;

# if 0
    // let's be cautious...
    if(mat->wt[j] != m){
	printf ("GASP: wt[%d]=%d != %d\n", j, mat->wt[j], m);
    }
# endif

    /* each m-merge leads to m-1 additions of rows */
#if DEBUG >= 1
    row_additions += m - 1;
    if (m > MERGE_LEVEL_MAX)
      {
	fprintf (stderr, "Error: m=%d > MERGE_LEVEL_MAX=%d\n",
                 m, MERGE_LEVEL_MAX);
	exit(1);
      }
    printf ("Treating column %d of weight %d\n",j,mat->wt[j]);
#if DEBUG >= 2
    printf ("Status before next j=%d to start\n", j);
    // the corresponding rows are in R[j], skipping 1st cell and -1's
    checkCoherence(mat, m, j);
#endif
#endif

    for(ni = 0, k = 1; k <= mat->R[j][0]; k++){
	if(mat->R[j][k] != UMAX(index_t)){
	    ind[ni++] = mat->R[j][k];
      if (ni == m)
              break; /* early abort, since we know there are m rows */
	}
    }
    /* now ind[0], ..., ind[m-1] are the m rows containing j */

#if DEBUG >= 1
    printf (" %d", j);
    printf (" => the %d rows are:\n", m);
    for(k = 0; k < m; k++){
	printf ("row[%d]=", ind[k]);
	print_row (mat, ind[k]);
	printf ("\n");
    }
    printf ("\n");
#endif

    *tt = seconds();
    findOptimalCombination (rep, mat, m, ind, j, tfill, tMST, useMST);
    *tt = seconds()-(*tt);
    mat->rem_nrows--;
    mat->rem_ncols--;
    removeColumnAndUpdate (mat, j);
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
    return matLengthRow(mat, i);
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

// this guy is linear => total is quadratic, be careful!
static int
findSuperfluousRows(int *tmp, int ntmp, filter_matrix_t *mat)
{
    index_t i, itmp;

    for(i = 0, itmp = 0; i < mat->nrows; i++)
        if(!isRowNull(mat, i)){
	    tmp[itmp++] = deleteScore(mat, i);
	    tmp[itmp++] = i;
	}
    // rows with largest score will be at the end
    qsort(tmp, itmp>>1, 2 * sizeof(int), cmp);
    return itmp;
}

// Delete at most niremmax superfluous rows such that nrows-ncols >= keep.
// Let's say we are at merge-level m.
// Return the number of removed rows.
static int
deleteSuperfluousRows (report_t *rep, filter_matrix_t *mat,
                       int niremmax, int m)
{
  int keep = mat->keep;
  int nirem = 0, *tmp, ntmp, i;

  if (m <= 2)
    return 0; /* it is not worth removing rows during level-2 merges */

  if ((int) (mat->rem_nrows - mat->rem_ncols) <= keep)
    return 0;

  if (niremmax > (int) (mat->rem_nrows - mat->rem_ncols) - keep)
    niremmax = (mat->rem_nrows - mat->rem_ncols) - keep;

  ntmp = mat->rem_nrows << 1;
  tmp = (int *) malloc (ntmp * sizeof(int));
  ntmp = findSuperfluousRows (tmp, ntmp, mat);
  // remove rows with largest score
  for (i = ntmp - 1; i >= 0 && nirem < niremmax; i -= 2)
    if (mat->rows[tmp[i]] != NULL)
      {
        removeRowDefinitely(rep, mat, tmp[i]);
        nirem++;
      }
  free (tmp);
  return nirem;
}

static double
my_cost (double N, double w, int forbw)
{
  if (forbw == 2)
    {
      double K1 = .19e-9, K2 = 3.4e-05, K3 = 1.4e-10; // kinda average
      return (K1+K3)*N*w+K2*N*log(N)*log(N);
    }
  else if (forbw == 3)
    return w / N;
  else if (forbw <= 1)
    return N * w;
  return 0.0;
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

int
number_of_superfluous_rows(filter_matrix_t *mat)
{
    int kappa, ni2rem;
    if (mat->keep > 0)
      kappa = (mat->rem_nrows-mat->rem_ncols) / mat->keep;
    else
      kappa = (mat->rem_nrows-mat->rem_ncols);

    if(kappa <= (1<<4))
      ni2rem = mat->keep / 2;
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

void
print_report (report_t *rep, filter_matrix_t *mat, int forbw, double bwcost)
{
  printf ("N=%d (%d) w=%lu", mat->rem_nrows, mat->rem_nrows - mat->rem_ncols,
                             mat->weight);
  if (forbw == 2)
    printf (" bw=%e", bwcost);
  else if (forbw == 3)
    printf (" w*N=%" PRIu64 "",
                        ((uint64_t)mat->rem_nrows) * ((uint64_t)mat->weight));
  else if (forbw <= 1)
    printf (" w*N=%e", bwcost);

  printf (" w/N=%2.2lf\n", ((double)mat->weight)/((double)mat->rem_nrows));

  if((forbw != 0) && (forbw != 3)) // what a trick!!!!
    fprintf (rep->outfile, "BWCOST: %1.0f\n", bwcost);

  fflush (stdout);
}


void
mergeOneByOne (report_t *rep, filter_matrix_t *mat, int maxlevel,
               int forbw, double ratio, double coverNmax)
{
    double totopt = 0.0, totfill = 0.0, totMST = 0.0, totdel = 0.0;
    double bwcostmin = 0.0, oldbwcost = 0.0, bwcost = 0.0;
    int old_ncols, m = 2, njrem = 0, ncost = 0, ncostmax, njproc;
    int ni2rem;
    int *nb_merges;
    int32_t dj, j, mkz;
    int useMST = 1; /* non-zero if we use minimal spanning tree */
    double REPORT = 20.0; /* threshold of w/N from which reports are done */
                          /* on stdout                                    */
    double FREQ_REPORT = 5.0; /* Once the threshold is exceeded, this is added */

    // clean things
    njrem = removeSingletons(rep, mat);

    nb_merges = (int*) malloc ((maxlevel + 1) * sizeof (int));
    for (m = 1; m <= maxlevel; m++)
      nb_merges[m] = 0;
    printf ("Using mergeOneByOne\n");
    ncostmax = 20; // was 5
    njproc = 0;

    print_report (rep, mat, forbw, bwcost);

    while(1){
	if(mat->itermax && (njproc >= mat->itermax)){
	    printf ("itermax=%d reached, stopping!\n", mat->itermax);
	    break;
	}
	oldbwcost = bwcost;
	old_ncols = mat->rem_ncols;
        if (MkzPopQueue(&dj, &mkz, mat) == 0)
          {
            printf ("Warning: heap is empty, increase maxlevel\n");
            break;
          }
	j = dj + 0;
        m = mat->wt[dj];
	if (m == 1) /* singleton ideal */
          removeColDefinitely(rep, mat, j);
	else if (m > 0) /* m=0 can happen for already merged ideals */
          mergeForColumn2(rep, mat, &njrem,
                          &totopt, &totfill, &totMST, &totdel, useMST, j);
        if (nb_merges[m]++ == 0 && m > 1)
          printf ("First %d-merge, cost %d (#Q=%d)\n", m, mkz,
                  MkzQueueCardinality(mat->MKZQ));
	// number of columns removed
	njproc += old_ncols - mat->rem_ncols;
	bwcost = my_cost ((double) mat->rem_nrows, (double) mat->weight,
                          forbw);
	//if (mat->rem_nrows % REPORT == 0)
  if ((((double)mat->weight)/((double)mat->rem_nrows)) > REPORT)
  {
      REPORT += FREQ_REPORT;
	    njrem = removeSingletons (rep, mat);
	    ni2rem = number_of_superfluous_rows (mat);
	    deleteSuperfluousRows (rep, mat, ni2rem, m);
      print_report (rep, mat, forbw, bwcost);
	}
	if((bwcostmin == 0.0) || (bwcost < bwcostmin)){
	    bwcostmin = bwcost;
	    if((forbw != 0) && (forbw != 3))
		// what a trick!!!!
		fprintf(rep->outfile, "BWCOST: %1.0f\n", bwcost);
	}
	// to be cleaned one day...
	if((forbw == 0) || (forbw == 2)){
          double r = bwcost / bwcostmin;
	    if(r > ratio){
		if((int) (mat->rem_nrows-mat->rem_ncols) > mat->keep){
		    // drop all remaining columns at once
		    ni2rem = mat->rem_nrows-mat->rem_ncols+mat->keep;
		    printf ("Dropping %d rows at once\n", ni2rem);
		    deleteSuperfluousRows(rep, mat, ni2rem, INT_MAX);
		}
		else{
		    if(forbw == 0)
			printf ("cN too high, stopping [%2.2lf]\n", r);
		    else
			printf ("bw too high, stopping [%2.2lf]\n", r);
		    break;
		}
	    }
	}
	else if (forbw == 3 && bwcost >= coverNmax)
          {
            index_t nrows = mat->rem_nrows;

            /* if the excess is still larger than what is wanted,
               remove the heaviest rows and loop again */
            printf ("remains %d rows\n", mat->rem_nrows);
            deleteSuperfluousRows (rep, mat,
                       (mat->rem_nrows - mat->rem_ncols) - mat->keep, INT_MAX);
            printf ("after deleteSuperfluousRows, remains %d rows\n",
                    mat->rem_nrows);
            removeSingletons (rep, mat);
            printf ("after removeSingletons, remains %d rows\n",
                    mat->rem_nrows);
            if (mat->rem_nrows == nrows)
              {
                printf ("w/N too high (%1.2f), stopping\n", bwcost);
                break;
              }
          }
	if((forbw == 1) && (oldbwcost != 0.0) && (bwcost > oldbwcost)){
	    ncost++;
#if 0
	    printf ("New cost > old cost (%.16e > %.16e) [%d/%d]\n",
		    bwcost, oldbwcost, ncost, ncostmax);
#endif
	    if(ncost >= ncostmax){
		int nirem;

		printf ("New cost > old cost %d times in a row:", ncost);
		nirem = deleteSuperfluousRows(rep, mat, 128, m);
		if(nirem == 0){
		    printf (" stopping\n");
		    break;
		}
		else{
		    printf (" try again after removing %d rows!\n", nirem);
		    njproc += nirem; // humf: odd name for njproc...!
		    continue;
		}
	    }
	}
	else
	    ncost = 0;
    }
    if (mat->itermax == 0)
      {
        printf ("Removing final excess, nrows=%d\n", mat->rem_nrows);
	deleteSuperfluousRows(rep, mat,
			      (mat->rem_nrows - mat->rem_ncols) - mat->keep,
                              INT_MAX);
        printf ("Removing singletons, nrows=%d\n", mat->rem_nrows);
        removeSingletons (rep, mat);
      }

#if DEBUG >= 1
    checkWeights (mat);
#endif

    if((forbw != 0) && (forbw != 3)){
	fprintf(rep->outfile, "BWCOSTMIN: %1.0f\n", bwcostmin);
	printf ("Minimal bwcost found: %1.0f\n", bwcostmin);
    }
#if DEBUG >= 1
    printf ("Total number of row additions: %lu\n", row_additions);
#endif
    for (m = 1; m <= maxlevel; m++)
      if (nb_merges[m] > 0)
        printf ("Number of %d-merges: %d\n", m, nb_merges[m]);
    free (nb_merges);
}

//////////////////////////////////////////////////////////////////////
//
// Resume section: very much inspired by replay.c, of course...!
//
//////////////////////////////////////////////////////////////////////

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
  int32_t j;
  int32_t ind[MERGE_LEVEL_MAX], i0;
  int ni, sg, k;

  ni = parse_hisfile_line (ind, str, &j);  
  
  if (ind[0] < 0)
    {
      sg = -1;
      i0 = -ind[0]-1;
    }
  else
    {
      sg = 1;
      i0 = ind[0];
    }

  for (k = 1; k < ni; k++)
      addRowsAndUpdate(mat, ind[k], i0, j);

  reportn(rep, ind, ni, j);
  
  if (sg > 0)
    {
      removeRowAndUpdate(mat, i0, 1);
      destroyRow(mat, i0);
      mat->rem_nrows--;
    }
}

// resumename is a file of the type mergehis.
// TODO: Compiles, but not really tested with Markowitz...!
void
resume(report_t *rep, filter_matrix_t *mat, const char *resumename)
{
    FILE *resumefile = fopen(resumename, "r");
    char str[RELATION_MAX_BYTES];
    unsigned long addread = 0;
    int nactivej;
    index_t j;
    char * rp;
    int REPORT = 10000;

    printf ("Resuming computations from %s\n", resumename);
    // skip first line containing nrows ncols
    rp = fgets(str, RELATION_MAX_BYTES, resumefile);
    ASSERT_ALWAYS(rp);

    printf ("Reading row additions\n");
    while(fgets(str, RELATION_MAX_BYTES, resumefile)){
	addread++;
	if((addread % (10 * REPORT)) == 0)
	    printf ("%lu lines read at %2.2lf\n", addread, seconds());
	if(str[strlen(str)-1] != '\n'){
	    printf ("Gasp: not a complete a line!");
	    printf (" I stop reading and go to the next phase\n");
	    break;
	}
	if(strncmp(str, "BWCOST", 6) != 0)
	    doAllAdds(rep, mat, str);
    }
    fclose(resumefile);
    for(j = 0; j < mat->ncols; j++)
	if((mat->wt[j] == 0) && MkzIsAlive(mat->MKZA, j))
	    // be sure j was removed...
	    removeColumnAndUpdate(mat, j);
    nactivej = 0;
    for(j = 0; j < mat->ncols; j++)
	if(mat->wt[j] != 0)
	    nactivej++;
    mat->rem_ncols = nactivej;
    printf ("At the end of resume, we have");
    printf (" nrows=%d ncols=%d (%d)\n",
	    mat->rem_nrows, mat->rem_ncols, mat->rem_nrows-mat->rem_ncols);
}
