/* 
 * Program: merge
 * Author : F. Morain
 * Purpose: merging relations
 * 
 * Algorithm: Cavallar++
 *
 */

// TODO: use compact lists...!

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>
#include <string.h>

#define DEBUG 0

#define USE_USED 0

#if USE_USED >= 1
char *used;
#endif

typedef struct {
  int len;    /* number of non-zero entries */
  int *val;   /* array of size len containing column indices of non-zero
                 entries */
} rel_t;

typedef struct {
  int nrows;
  int ncols;
  int rem_weight;
  int rem_nrows;     /* number of remaining rows */
  int rem_ncols;     /* number of remaining columns */
  rel_t *data;
  unsigned long* wt; /* weight of prime j, <= 1 for a deleted prime */
  unsigned long* ad; /* ad[j] contains the sum of indices of all rows
                        containing prime j; when wt[j]=2, and we know
                        that row i contains j, then the other row is
                        simply k = ad[j] - i. Overflow is not a problem,
                        since the value is always exact mod 2^32 or 2^64. */
} sparse_mat_t;

static int 
cmp(const void *p, const void *q) {
    int x = *((int *)p);
    int y = *((int *)q);
    return (x <= y ? -1 : 1);
}

void
inspectWeight(int *usecol, FILE *purgedfile, int nrows, int ncols, int wmax)
{
    int i, j, ret, x, nwlarge, nc;

    memset(usecol, 0, ncols * sizeof(int));
    for(i = 0; i < nrows; i++){
	ret = fscanf(purgedfile, "%d", &j); // unused index to rels file
	assert (ret == 1);
	ret = fscanf (purgedfile, "%d", &nc);
	assert (ret == 1);
	for(j = 0; j < nc; j++){
	    ret = fscanf(purgedfile, "%d", &x);
	    assert (ret == 1);
	    usecol[x]++;
	}
    }
    nwlarge = 0;
    for(j = 0; j < ncols; j++)
	if(usecol[j] > wmax){
	    usecol[j] = 0;
	    nwlarge++;
	}
    fprintf(stderr, "#(w > %d) = %d\n", wmax, nwlarge);
}

/* Reads a matrix file, and puts in mat->wt[j], 0 <= j < ncols, the
   weight of column j (adapted from matsort.c), and in mat->ad[j],
   0 <= j < ncols, the sum of row indices where j appears. 

   We skip columns that are too heavy, that is columns for which usecol[j]==0.

*/
void
readmat(FILE *file, sparse_mat_t *mat, int *usecol, int wmax)
{
  int ret;
  int i, j, l2, *w;
  int nc, x, buf[100], ibuf;

  w = (int *)malloc((wmax+1) * sizeof(int));
  memset(w, 0, (wmax+1) * sizeof(int));

  ret = fscanf (file, "%d %d", &(mat->nrows), &(mat->ncols));
  assert (ret == 2);

  fprintf (stderr, "Reading matrix of %d rows and %d columns: excess is %d\n",
           mat->nrows, mat->ncols, mat->nrows - mat->ncols);
  mat->rem_nrows = mat->nrows;

  mat->wt = (unsigned long*) malloc (mat->ncols * sizeof (unsigned long));
  mat->ad = (unsigned long*) malloc (mat->ncols * sizeof (unsigned long));
  for (j = 0; j < mat->ncols; j++)
      mat->wt[j] = mat->ad[j] = 0;

  mat->rem_weight = 0;
  mat->data = (rel_t*) malloc (mat->nrows * sizeof (rel_t));

  for (i = 0; i < mat->nrows; i++){
      ret = fscanf(file, "%d", &j); // unused index to rels file
      assert (ret == 1);
      ret = fscanf (file, "%d", &nc);
      assert (ret == 1);
      if(nc == 0){
	  mat->data[i].len = 0;
	  mat->data[i].val = NULL;
	  mat->rem_nrows--;
      }
      else{
	  for(j = 0, ibuf = 0; j < nc; j++){
	      ret = fscanf(file, "%d", &x);
#if DEBUG >= 1
	      fprintf(stderr, "i = %d, j = %d, x = %d\n", i, j, x);
#endif
	      assert (ret == 1);
	      assert (0 <= x && x < mat->ncols);
	      if(usecol[x]){
		  buf[ibuf++] = x;
		  mat->wt[x] ++;
		  mat->rem_weight++;
		  mat->ad[x] += i; /* computed mod 2^32 or 2^64 */
	      }
	  }
	  mat->data[i].len = ibuf;
	  mat->data[i].val = (int*) malloc (ibuf * sizeof (int));
	  memcpy(mat->data[i].val, buf, ibuf * sizeof(int));
	  // sort indices in val to ease row merges
	  qsort(mat->data[i].val, ibuf, sizeof(int), cmp);
      }
    }

  /* count the number of primes appearing at least twice */
  for(j = l2 = 0; j < mat->ncols; j++){
      if(mat->wt[j] >= 2){
	  l2++;
	  if(mat->wt[j] <= wmax)
	      w[mat->wt[j]]++;
      }
      else if(usecol[j])
	  fprintf(stderr, "Warning: prime %d appears only %lu time(s)\n", j,
		  mat->wt[j]);
  }
  fprintf (stderr, "Found %d primes appearing at least twice\n", l2);
  for(j = 2; j <= wmax; j++)
      if(w[j] > 0)
	  fprintf (stderr, "Found %d primes of weight %d\n", w[j], j);
  mat->rem_ncols = l2;
  free(w);
}

void
print_row(sparse_mat_t *mat, int i)
{
    int j;
    
    for(j = 0; j < mat->data[i].len; j++)
	fprintf(stderr, " %d", mat->data[i].val[j]);
}

void
matrix2tex(sparse_mat_t *mat)
{
    int i, j, k, *tab;

    tab = (int *)malloc(mat->ncols * sizeof(int));
    fprintf(stderr, "\\begin{array}{l|");
    for(j = 0; j < mat->ncols; j++)
	fprintf(stderr, "c");
    fprintf(stderr, "}\n");
    for(i = 0; i < mat->nrows; i++){
	memset(tab, 0, mat->ncols * sizeof(int));
	for(k = 0; k < mat->data[i].len; k++)
	    tab[mat->data[i].val[k]] = 1;
	fprintf(stderr, "R_{%d}", i);
	for(j = 0; j < mat->ncols; j++)
	    fprintf(stderr, "& %d", tab[j]);
	fprintf(stderr, "\\\\\n");
    }
    fprintf(stderr, "\\hline\n");
    fprintf(stderr, "\\text{weight}");
    for(j = 0; j < mat->ncols; j++)
	fprintf(stderr, "& %lu", mat->wt[j]);
    fprintf(stderr, "\\\\\n");
    fprintf(stderr, "\\end{array}\n");
    free(tab);
}

void
destroy_row(sparse_mat_t *mat, int i)
{
    free(mat->data[i].val);
    mat->data[i].val = NULL;
    mat->data[i].len = 0;
}

void
removeWeight(sparse_mat_t *mat, int i)
{
    int k;

    for(k = 0; k < mat->data[i].len; k++)
	mat->wt[mat->data[i].val[k]]--;
}

// i1 += i2, mat->wt is updated at the same time.
void
addrows(sparse_mat_t *mat, int i1, int i2)
{
    int k1, k2, k, len, *tmp, *tmp2;

#if DEBUG >= 2
    fprintf(stderr, "row[%d] =", i1); print_row(mat, i1);
    fprintf(stderr, "\n");
    fprintf(stderr, "row[%d] =", i2); print_row(mat, i2);
    fprintf(stderr, "\n");
#endif
    // merge row[i1] and row[i2] in i1...
    len = mat->data[i1].len + mat->data[i2].len;
    tmp = (int *)malloc(len * sizeof(tmp));
    k = k1 = k2 = 0;

    // TODO: how do we update ad????

    // i1 is to disappear, replaced by a new one
    removeWeight(mat, i1);
    // loop while everybody is here
    while((k1 < mat->data[i1].len) && (k2 < mat->data[i2].len)){
	if(mat->data[i1].val[k1] < mat->data[i2].val[k2]){
	    tmp[k++] = mat->data[i1].val[k1++];
	    mat->wt[mat->data[i1].val[k1-1]]++;
	}
	else if(mat->data[i1].val[k1] > mat->data[i2].val[k2]){
            tmp[k++] = mat->data[i2].val[k2++];
	    mat->wt[mat->data[i2].val[k2-1]]++;
	}
	else{
#if DEBUG >= 1
	    fprintf(stderr, "ADD[%d=(%ld,%lu), %d=(%ld,%lu)]: new w[%d]=%lu\n", 
		    i1, mat->data[i1].a, mat->data[i1].b,
		    i2, mat->data[i2].a, mat->data[i2].b,
		    mat->data[i1].val[k1], 
		    mat->wt[mat->data[i1].val[k1]]);
#endif
	    mat->rem_weight -= 2; // enough?
	    k1++;
	    k2++;
	}
    }
    // finish with k1
    for( ; k1 < mat->data[i1].len; k1++){
	tmp[k++] = mat->data[i1].val[k1];
	mat->wt[mat->data[i1].val[k1]]++;
    }
    // finish with k2
    for( ; k2 < mat->data[i2].len; k2++){
	tmp[k++] = mat->data[i2].val[k2];
	mat->wt[mat->data[i2].val[k2]]++;
    }
    // destroy and copy back
    free(mat->data[i1].val);
    tmp2 = (int *)malloc(k * sizeof(int));
    memcpy(tmp2, tmp, k * sizeof(int));
    mat->data[i1].len = k;
    mat->data[i1].val = tmp2;
    free(tmp);
#if DEBUG >= 2
    fprintf(stderr, "row[%d]+row[%d] =", i1, i2);
    print_row(mat, i1); fprintf(stderr, "\n");
#endif
}

void
merge2rows(sparse_mat_t *mat, int j)
{
    int i1, i2, ok = 0, k1, k2;

    for(i1 = 0; i1 < mat->nrows; i1++){
	if(mat->data[i1].val == NULL)
	    continue;
	for(k1 = 0; k1 < mat->data[i1].len; k1++)
	    if(mat->data[i1].val[k1] == j){
		ok = 1;
		break;
	    }
	if(ok)
	    break;
    }
    if(!ok){
	printf("Pb in merge2rows: no row containing %d\n", j);
	return;
    }
#if 0 // one day...
    i2 = mat->ad[j]-i1;
#else
    ok = 0;
    for(i2 = i1+1; i2 < mat->nrows; i2++){
	if(mat->data[i2].val == NULL)
	    continue;
	for(k2 = 0; k2 < mat->data[i2].len; k2++)
            if(mat->data[i2].val[k2] == j){
                ok = 1;
                break;
            }
	if(ok)
	    break;
    }
    if(!ok){
	printf("Pb in merge2rows: no 2nd row containing %d\n", j);
	return;
    }
#endif
#if USE_USED >= 1
    if(used[i1])
	fprintf(stderr, "R%d used %d times\n", i1, (int)used[i1]);
    if(used[i2])
	fprintf(stderr, "R%d used %d times\n", i2, (int)used[i2]);
    used[i1]++;
    used[i2] = -1;
#endif
#if DEBUG >= 1
    fprintf(stderr, "Merging rows: %d += %d [j=%d]\n", i1, i2, j);
#endif
    printf("%d %d\n", i2, i1); // note the order!!!
    addrows(mat, i1, i2);
    removeWeight(mat, i2);
    destroy_row(mat, i2);
    mat->rem_nrows--;
    mat->rem_ncols--;
#if DEBUG >= 1
    matrix2tex(mat);
#endif
}

void
merge2(sparse_mat_t *mat, int nb_merge_max)
{
    int j, nb_merge = 0;

#if DEBUG >= 1
    matrix2tex(mat);
#endif
    for(j = 0; j < mat->ncols; j++){
	if(nb_merge == nb_merge_max){
	    fprintf(stderr, "Warning: nb_merge has reached its maximum:");
	    fprintf(stderr, " %d\n", nb_merge_max);
	    break; // useful for debugging and more...!
	}
	if(mat->wt[j] == 2){
#if DEGUB >= 1
	    fprintf(stderr, "j = %d has weight 2\n", j);
#endif
	    nb_merge++;
	    merge2rows(mat, j);
	    if(!(nb_merge % 1000))
		fprintf(stderr, "nrows=%d ncols=%d weight=%d\n",
			mat->rem_nrows, mat->rem_ncols, mat->rem_weight);
	}
    }
}

// what is the weight of the sum of Ra and Rb?
int
weightSum(sparse_mat_t *mat, int i1, int i2)
{
    int k1 = 0, k2 = 0, w = 0;

    while((k1 < mat->data[i1].len) && (k2 < mat->data[i2].len)){
	if(mat->data[i1].val[k1] < mat->data[i2].val[k2]){
	    k1++;
	    w++;
	}
	else if(mat->data[i1].val[k1] > mat->data[i2].val[k2]){
	    k2++;
	    w++;
	}
	else{
	    k1++; k2++;
	}
    }
    w += (k1 > mat->data[i1].len ? 0 : mat->data[i1].len-k1+1);
    w += (k2 > mat->data[i2].len ? 0 : mat->data[i2].len-k2+1);
    return w;
}

int
findAllRowsWithGivenj(int *ind, sparse_mat_t *mat, int j, int nb)
{
    int i, k, r = 0;

    for(i = 0; i < mat->nrows; i++){
	if(mat->data[i].val == NULL)
	    continue;
	for(k = 0; k < mat->data[i].len-1; k++) // trick!
	    if(mat->data[i].val[k] >= j)
		break;
	if(mat->data[i].val[k] == j){
	    ind[r++] = i;
	    if(r == nb)
		return 1;
	}
    }
    return 0;
}

int
findBestIndex3(sparse_mat_t *mat, int *ind)
{
    int w01, w02, w12, wtot, w0, w1, w2, k, i;

    w01 = weightSum(mat, ind[0], ind[1]);
    w02 = weightSum(mat, ind[0], ind[2]);
    w12 = weightSum(mat, ind[1], ind[2]);
#if DEBUG >= 1
    fprintf(stderr,"W(0+1)=%d\n", w01);
    fprintf(stderr,"W(0+2)=%d\n", w02);
    fprintf(stderr,"W(1+2)=%d\n", w12);
#endif
    // first scenario r1 += r0, r2 += r0, r0 disappears
    wtot = 0;
    for(k = 0; k < 3; k++)
	wtot += mat->data[ind[k]].len;
    w0 = w01+w02-wtot;
    w1 = w01+w12-wtot;
    w2 = w02+w12-wtot;
#if DEBUG >= 1
    fprintf(stderr, "Using r0 = r[%d]: %d\n", ind[0], w0);
    fprintf(stderr, "Using r1 = r[%d]: %d\n", ind[1], w1);
    fprintf(stderr, "Using r2 = r[%d]: %d\n", ind[2], w2);
#endif
    if(w0 < w1)
	i = (w0 < w2 ? 0 : 2);
    else
	// w0 >= w1
	i = (w1 < w2 ? 1 : 2);
    return i;
}

int
findBestIndex(sparse_mat_t *mat, int m, int *ind)
{
    int **A, i, j;

    if(m == 3)
	return findBestIndex3(mat, ind);
    fprintf(stderr, "Being implemented...!\n");
    A = (int **)malloc(m * sizeof(int *));
    for(i = 0; i < m; i++){
	A[i] = (int *)malloc(m * sizeof(int));
	A[i][i] = 0;
    }
    // A[i][j] <- Weight(R[ind[i]]+R[ind[j]]);
    for(i = 0; i < m; i++)
	for(j = i+1; j < m; j++){
	    A[i][j] = weightSum(mat, ind[i], ind[j]);
	    A[j][i] = A[i][j];
	}
    for(i = 0; i < m; i++)
	free(A[i]);
    free(A);
    return -1;
}

void
merge_m(sparse_mat_t *mat, int m)
{
    int *ind, j, k, i;

#if DEBUG >= 1
    fprintf(stderr, "Weight 3:");
#endif
    ind = (int *)malloc(m * sizeof(int));
    for(j = 0; j < mat->ncols; j++){
	if(mat->wt[j] != m)
	    continue;
	// we need to find the three rows and then
	if(!findAllRowsWithGivenj(ind, mat, j, m))
	    fprintf(stderr, "Could not find the %d required rows\n", m);
	// try all combinations to find the smaller one
#if DEBUG >= 1
	fprintf(stderr, " %d", j);
	fprintf(stderr, "=> the %d rows are:\n", j, m);
	for(k = 0; k < m; k++){
	    fprintf(stderr, "row[%d]=", ind[k]);
	    print_row(mat, ind[k]);
	    fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
#endif
	i = findBestIndex(mat, m, ind);
#if DEBUG >= 1
	fprintf(stderr, "Minimal is i=%d\n", i);
#endif
	// terrific hack: everybody on the same line
	// the first is to be destroyed in replay!!!
	printf("%d", ind[i]);
	for(k = 0; k < m; k++)
	    if(k != i){
		addrows(mat, ind[k], ind[i]);
#if DEBUG >= 1
		fprintf(stderr, "new row[%d]=", ind[k]);
		print_row(mat, ind[k]);
		fprintf(stderr, "\n");
#endif
		printf(" %d", ind[k]);
	    }
	printf("\n");
	removeWeight(mat, ind[i]);
	destroy_row(mat, ind[i]);
	mat->rem_nrows--;
	mat->rem_ncols--;
	mat->wt[j] = 0;
    }
    free(ind);
}

void
mergeGe3(sparse_mat_t *mat, int m)
{
    int j, nbm;

    for(j = 0, nbm = 0; j < mat->ncols; j++)
	if((mat->wt[j] > 0) && (mat->wt[j] < m))
	    fprintf(stderr, "#W# wt[%d] = %ld\n", j, mat->wt[j]);
        else if(mat->wt[j] == m){
	    //	    fprintf(stderr, "# wt[%d] = %d\n", j, m);
	    nbm++;
	}
    fprintf(stderr, "There are %d column(s) of weight %d\n", nbm, m);
    if(m == 3)
    	merge_m(mat, m);
#if DEBUG >= 1
    matrix2tex(mat);
#endif
}

void
merge(sparse_mat_t *mat, int nb_merge_max, int maxlevel)
{
    int old_nrows, old_ncols;

    printf("%d %d\n", mat->nrows, mat->ncols);
    do{
	old_nrows = mat->rem_nrows;
	old_ncols = mat->rem_ncols;
	merge2(mat, nb_merge_max);
	fprintf(stderr, "=> nrows=%d ncols=%d weight=%d\n",
		mat->rem_nrows, mat->rem_ncols, mat->rem_weight);
	break; // TODO: destroy this?
    } while((old_nrows != mat->rem_nrows) && (old_ncols != mat->rem_ncols));
    if(maxlevel >= 3)
	mergeGe3(mat, maxlevel);
}

void
checkmat(sparse_mat_t *mat)
{
    int i;

    for(i = 0; i < mat->nrows; i++){
	fprintf(stderr, "Row %d =", i);
	print_row(mat, i);
	fprintf(stderr, "\n");
    }
    exit(-1);
}

// not activated anymore...!
void
dumpSparse(FILE *ofile, sparse_mat_t *mat, int *code)
{
    int i, j, k, buf[1000], ibuf, new_nrows = 0;

    fprintf(ofile, "%d %d\n", mat->rem_nrows, mat->rem_ncols);
    for(i = 0; i < mat->nrows; i++){
	if(mat->data[i].val == NULL)
	    continue;
	new_nrows++;
	ibuf = 0;
	for(k = 0; k < mat->data[i].len; k++){
	    j = mat->data[i].val[k];
	    if(code[j])
		buf[ibuf++] = code[j]-1;
	}
	fprintf(ofile, "%d", ibuf);
	for(k = 0; k < ibuf; k++)
	    fprintf(ofile, " %d", buf[k]);
	fprintf(ofile, "\n");
    }
    assert(new_nrows == mat->rem_nrows);
}

int
main (int argc, char *argv[])
{
  FILE *purgedfile;
  sparse_mat_t mat;
  char *purgedname = NULL, *hisname = NULL;
  int nb_merge_max = 0, domerge = 0, *usecol, nrows, ncols;
  int wmax = 20, maxlevel = 2;

  fprintf (stderr, "%s rev. %s\n", argv[0], REV);

  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 2 && strcmp (argv[1], "-merge") == 0)
	{
	  domerge = 1;
	  nb_merge_max = atoi(argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-mat") == 0)
	{
	  purgedname = argv[2];
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-rebuild") == 0)
	{
	  hisname = argv[2];
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-wmax") == 0)
	{
	  wmax = atoi(argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-maxlevel") == 0)
	{
	  maxlevel = atoi(argv[2]);
	  argc -= 2;
	  argv += 2;
	}
      else 
	break;
    }

  purgedfile = fopen(purgedname, "r");
  fscanf(purgedfile, "%d %d", &nrows, &ncols);
  usecol = (int *)malloc(ncols *sizeof(int));
  inspectWeight(usecol, purgedfile, nrows, ncols, wmax);

  rewind(purgedfile);
  readmat(purgedfile, &mat, usecol, wmax);
  fclose(purgedfile);
#if DEBUG >= 3
  checkmat(&mat);
#endif
  if(domerge){
#if USE_USED >= 1
      used = (char *)malloc(mat.nrows * sizeof(char));
      memset(used, 0, mat.nrows * sizeof(char));
#endif
      merge(&mat, nb_merge_max, maxlevel);
  }
  return 0;
}
