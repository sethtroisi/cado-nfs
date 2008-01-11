/* 
 * Program: gauss
 * Author : P. Gaudry
 * Purpose: solve linear systems over GF(2)
 * 
 * Algorithm: Gaussian elimination as described in ...
 *
 */


/*
Compile with one of these:

gcc -O4 -DNDEBUG -funroll-loops -DMULTI_ROW=3 -DNB_TEST=10 gauss.c -lgmp

gcc gauss.c -lgmp


Note: for small sizes (100-200), MULTI_ROW=2 seems to be fine,
      for larger sizes (also for small sizes on alpha), try MULTI_ROW=3

      With NO_MAIN option, a .o file is generated, exporting only the
      kernel function.  In that case, the code does not use GMP at all.
      
The exported function is the following:

 int kernel(mp_limb_t* mat, mp_limb_t** ker, int nrows, int ncols,
            int limbs_per_row, int limbs_per_col);

which returns the dimension of the kernel and fill-in the pointer ker.
If just the dimension of kernel is wanted, set ker=NULL.

*/

/*===========================================================================*/
/*    includes, defines, and declarations...                                 */
/* Warning: there are some (static) global variables                         */
/*===========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "gmp.h"   /* only used for setting a random matrix */
#include "utils.h"  /* for seconds() */
#include "gauss.h"

#ifndef VERBOSE
#define VERBOSE 0
#endif
#ifndef NB_TEST
#define NB_TEST 100
#endif
#ifndef MULTI_ROW
#define MULTI_ROW 2
#endif


/* Look on which computer we are */
#ifdef __alpha__
#define MACHINE_WORD_SIZE (64)
#define BYTE_PER_WORD (8)
#else
#ifdef __i386__
#define MACHINE_WORD_SIZE (32)
#define BYTE_PER_WORD (4)
#else
#ifdef __x86_64__
#define MACHINE_WORD_SIZE (64)
#define BYTE_PER_WORD (8)
#else
#error unknown cpu
#endif
#endif
#endif

#ifndef __GMP_H__
typedef unsigned long int mp_limb_t;
#endif

#define INLINE __inline__


/* some global variables storing the data */
static int NROWS = 0;
static int NCOLS = 0;
static int LIMBS_PER_ROW = 0;
static mp_limb_t* matrix = NULL;
static mp_limb_t** ptr_rows = NULL;

/* functions declaration */

static void check_soundness();

#if 0 /* ununsed currently */
static void addPartialRows(int row, int pivot, mp_limb_t mask,
			   int j_cur, mp_limb_t **ptr);
#endif
INLINE static void addRows(int row, int pivot, mp_limb_t mask,
			   /*int j_cur,*/ mp_limb_t **ptr);
INLINE static void add2Rows(int row, int row2, int pivot, mp_limb_t mask,
			    /*int j_cur,*/ mp_limb_t **ptr_current);
INLINE static void add3Rows(int row, int row2, int row3, int pivot,
			    mp_limb_t mask, /* int j_current, */
			    mp_limb_t **ptr_current);
INLINE static int getPivot(mp_limb_t **ptr, mp_limb_t mask);
#if VERBOSE
static void printVector(mp_limb_t *V, int n);
static void printMatrix(mp_limb_t *M, int nrows, int ncols, int limbs_per_row);
#endif


/*===========================================================================*/
/*  functions definition                                                     */
/*===========================================================================*/

#ifndef NO_MAIN
int main(int argc, char **argv) {
  mp_limb_t* mat;
  mp_limb_t** ker;
  int nrows, ncols, limbs_per_row, limbs_per_col;
  int i, j, justrank = 0;
  

  if (argc > 1) {
    nrows = atoi(argv[1]);
    if ((nrows < 1) || (nrows > 1000000)) {
      printf("usage: %s n [justrank]\n", argv[0]);
      printf("   where n is the size of the matrix\n");
      exit(1);
    }
    ncols = nrows;
    if (argc > 2) {
      justrank = 1;
      if (strcmp(argv[2], "justrank")!=0) {
	printf("usage: %s n [justrank]\n", argv[0]);
	printf("   where n is the size of the matrix\n");
	exit(1);
      }
    }
  } else {
    nrows = 163;
    ncols = 163;
  }
  limbs_per_row = (ncols / MACHINE_WORD_SIZE) + 1;
  limbs_per_col = (nrows / MACHINE_WORD_SIZE) + 1;

  mat = (mp_limb_t *)malloc(limbs_per_row*nrows*sizeof(mp_limb_t));
  ker = (mp_limb_t **)malloc(nrows*sizeof(mp_limb_t *));
  for (i = 0; i < nrows; ++i)
    ker[i] = (mp_limb_t *)malloc(limbs_per_col*sizeof(mp_limb_t));

  for (i = 0; i < NB_TEST; ++i) {
    int dim;
    mpn_random(mat, limbs_per_row*nrows);

#if VERBOSE
    printf("Matrix is :\n");
    printMatrix(mat, nrows, ncols, limbs_per_row);
#endif

    if (!justrank)
      dim = kernel(mat, ker, nrows, ncols, limbs_per_row, limbs_per_col);
    else
      dim = kernel(mat, NULL, nrows, ncols, limbs_per_row, limbs_per_col);

#if VERBOSE
    printf("Echelonized matrix is :\n");
    printMatrix(mat, nrows, ncols, limbs_per_row);
#endif
    
    printf("dim of ker = %d\n", dim);

#if VERBOSE
    if ((!justrank) && (dim > 0)) 
      printf("Basis of kernel given by:\n");
      for (j = 0; j < dim; ++j) {
	printVector(ker[j], nrows); printf("\n");
      }
#endif

  }

  free(mat);
  free(ker);
}
#endif

static void check_soundness(){
#ifndef NDEBUG
  mp_limb_t x = 2342326;
  int i;
  assert (LIMBS_PER_ROW*MACHINE_WORD_SIZE >= NCOLS);
  /* touch everywhere in the matrix to look for SEGV's */
  for (i = 0; i < NROWS; ++i) 
    x += matrix[LIMBS_PER_ROW*i];
  if (x == matrix[LIMBS_PER_ROW*NROWS-1])
    assert(0);  /* very unlikely to happen, but forces no optimizations */
#endif
}

/* 
 * Kernel computation.
 *
 * Put the vectors of the kernel in the (preallocated) table of cols ker.
 * Return the dimension of the kernel.
 * If just the dimension is wanted, then pass ker = NULL, and in that case,
 * limbs_per_cols is meaningless.
 * 
 * Variables description:
 *   ptr_current : table containing for each row the address of
 *                 the limb which we are currently dealing with.
 *                 These are incremented after MACHINE_WORD_SIZE 
 *                 pivot elimination.
 *   col_current : index of the current column
 *   j_current   : the index of the limb containing the current column
 *   mask1       : mask with only one 1 set at the current column
 *   mask2       : mask with 1's set at column higher than the current
 *   set_pivot   : table. at index j, contains the index of the row that
 *                 used for pivoting column j. Initialized with -1.
 *   set_used    : table. at index i, contains 1 if the i-th row has
 *                 been used for pivoting, at some point, 0 otherwise.
 *   rank        : the rank of the matrix
 *   i, j        : used in for(;;) loops
 */


int kernel(mp_limb_t* mat, mp_limb_t** ker, int nrows, int ncols,
	   int limbs_per_row, int limbs_per_col) {
  int i, j, rank;
  mp_limb_t **ptr_current;
  int *set_pivot;
  int *set_used;
  int j_current, col_current;
  mp_limb_t mask1, mask2;
  double st0 = seconds(), st;

  fprintf (stderr, "Using MULTI_ROW=%d\n", MULTI_ROW);

  /* store the data in the global variables */
  matrix = mat;
  NROWS = nrows;
  NCOLS = ncols;
  LIMBS_PER_ROW = limbs_per_row;
  ptr_rows = (mp_limb_t **)malloc((NROWS+1)*sizeof(mp_limb_t *));
  for (i = 0; i < NROWS+1; ++i)
    ptr_rows[i] = matrix + LIMBS_PER_ROW*i;
  check_soundness();
  
  /* initialize ptr_current, set_used and set_pivot */
  ptr_current = (mp_limb_t **)malloc((NROWS)*sizeof(mp_limb_t *));
  for (i = 0; i < NROWS; ++i)
    ptr_current[i] = ptr_rows[i];
  set_pivot = (int *)malloc(NCOLS*sizeof(int));
  set_used  = (int *)malloc(NROWS*sizeof(int));
  for (i = 0; i < NCOLS; ++i)
    set_pivot[i] = -1;
  for (i = 0; i < NROWS; ++i)
    set_used[i] = 0;

  /* Main Loop: for each column, find the pivot and eliminate */
  j_current = 0;
  mask1 = 1UL;
  mask2 = (-(1UL)) ^ (1UL);
  col_current = 0;
  while (col_current < NCOLS) {
    int pivot;

    /* Get the pivot */
    pivot = getPivot(ptr_current, mask1);
    if (pivot != -1) {
#ifndef NDEBUG
      for (i = 0; i < pivot; ++i)
	assert (!((*ptr_current[i]) & mask1 ));
      assert (((*ptr_current[pivot]) & mask1 ));
      assert (set_used[pivot] == 0);
#endif
      set_pivot[col_current] = pivot;
      set_used[pivot] = 1;

      /* Eliminate with the pivot */
#if MULTI_ROW == 1
      for (i = pivot + 1; i < NROWS; ++i) {
	/* is there a 1 in position (i, col_current) ? */
	if ( (*ptr_current[i]) & mask1 )   
	  addRows(i, pivot, mask1, /*j_current,*/ ptr_current);
      }      
#elif MULTI_ROW == 2
      /* try with two rows at a time */
      { 
	int i1;
	i = pivot + 1;
	i1 = -1;
	while (i < NROWS) {
	  if ( (*ptr_current[i]) & mask1 ) {
	    if (i1 >= 0) {
	      add2Rows(i1, i, pivot, mask1, /*j_current,*/ ptr_current);
	      i1 = -1;
	    } else {
	      i1 = i;
	    }
	  } /* end if */
	  ++i;
	} /* end while */
	if (i1 >= 0) 
	  addRows(i1, pivot, mask1, /*j_current,*/ ptr_current);
      } /* end block */
#elif MULTI_ROW == 3
      /* try with three rows at a time */
      { 
	int i1, i2;
	i = pivot + 1;
	i1 = -1;  i2 = -1;
	while (i < NROWS) {
	  if ( (*ptr_current[i]) & mask1 ) {
	    if (i1 >= 0) {
	      if (i2 >= 0) {
		add3Rows(i1, i2, i, pivot, mask1, /*j_current,*/ ptr_current);
		i1 = -1; i2 = -1;
	      }
	      else {
		i2 = i;
	      }
	    } else {
	      i1 = i;
	    }
	  } /* end if */
	  ++i;
	} /* end while */
	if (i1 >= 0)  {
	  if (i2 >= 0)
	    add2Rows(i1, i2, pivot, mask1, /*j_current,*/ ptr_current);
	  else
	    addRows(i1, pivot, mask1, /*j_current,*/ ptr_current);
	}
      } /* end block */
#else
#error MULTI_ROW must be 1, 2, or 3.
#endif 

      /* Purge the pivot line */
      addRows(pivot, pivot, mask1, /*j_current,*/ ptr_current);
    }

    /* increment */
    ++col_current;
    mask1 <<= 1;
    mask2 <<= 1;
    if (!mask1) { /* we have done MACHINE_WORD_SIZE operations */
      mask1 = 1UL;
      mask2 = (-(1UL)) ^ (1UL);
      j_current++;
      assert (j_current <= LIMBS_PER_ROW);
      for (i = 0; i < NROWS; ++i)
	ptr_current[i]++;
    }      

    /* some verbosity... */
    if ((col_current % 128) == 0)
      {
        st = (seconds () - st0); /* time in seconds */
        fprintf (stderr, "done %d pivots in %1.0fs (est. %1.0fs)\n",
			col_current,
                 st, (double) NCOLS * st / (double) col_current);
      }

  } /* end while */

#if VERBOSE
  printf("Pivots are: [");
  for (i = 0; i < NCOLS-1; ++i)
    printf("%d,", set_pivot[i]);
  printf("%d]\n", set_pivot[NCOLS-1]);

  printf("Used rows are: [");
  for (i = 0; i < NROWS-1; ++i)
    printf("%d,", set_used[i]);
  printf("%d]\n", set_used[NROWS-1]);
#endif     
 
  assert (limbs_per_col*MACHINE_WORD_SIZE >= NROWS);
  /* Explore the set of unused rows, get the rank */
  rank = NROWS;
  for (i = 0; i < NROWS; ++i)
    if (!set_used[i]) {
      rank--;
      if (ker != NULL) {
	/* each unused row gives one element of the kernel */
#if SAVE_KERNEL_MEMORY
	ker[0] = (mp_limb_t *) malloc(limbs_per_col * sizeof(mp_limb_t));
#endif

	/* ker[0] contains the space for the new element of the kernel */
	mp_limb_t *ptr = ptr_rows[i];
	
	for (j = 0; j < limbs_per_col; ++j)
	  ker[0][j] = 0;
	ker[0][i / MACHINE_WORD_SIZE] = 1UL << (i % MACHINE_WORD_SIZE);

	mask1 = 1UL;
	j = 0;
	while (j < NCOLS) {
	  if ((*ptr) & mask1) { /* have to recover the pivoting */
	    int pivot;
	    pivot = set_pivot[j];
	    (ker[0])[pivot / MACHINE_WORD_SIZE] |= 
	      (1UL << (pivot % MACHINE_WORD_SIZE));
	  } /* end if */
	  /* increment */
	  ++j;
	  mask1 <<= 1;
	  if (!mask1) { /* we have done MACHINE_WORD_SIZE operations */
	    mask1 = 1UL;
	    ptr++;
	  }
	} /* end while */
	ker++;
      } /* end if ker != NULL */
    } /* end if / for, loop over all unused */

  free(ptr_rows);
  free(ptr_current);
  free(set_pivot);
  free(set_used);

  return NROWS - rank;
} /* end function kernel */



/* add the pivot row to the given row, expect for the given column 
 *   mask gives the bit we have to keep
 *   j_current is the number of the limb which contains the column
 */

INLINE static void addRows(int row, int pivot, mp_limb_t mask,
                           /*int j_current,*/  mp_limb_t **ptr_current) {
  int i;
  mp_limb_t *ptr1, *ptr2;

  ptr1 = ptr_rows[row];
  ptr2 = ptr_rows[pivot];

  for (i = 0; i < LIMBS_PER_ROW; ++i) 
    ptr1[i] ^= ptr2[i];
  *ptr_current[row] |= mask; 
} /* end function addRows */


/* add the pivot row to the two given rows, expect for the given column 
 *   mask gives the bit we have to keep
 *   j_current is the number of the limb which contains the column
 */

INLINE static void add2Rows(int row, int row2, int pivot, mp_limb_t mask,
			    /*int j_current,*/ mp_limb_t **ptr_current)
{
  // int i;
  mp_limb_t *ptr1, *ptr2, *ptr3, *ptr_lim;

  ptr1 = ptr_rows[row];
  ptr2 = ptr_rows[row2];
  ptr3 = ptr_rows[pivot];
  ptr_lim = ptr_rows[pivot + 1];
  while (ptr3 != ptr_lim) {
    *ptr1++ ^= *ptr3;
    *ptr2++ ^= *ptr3++;
  }

  *ptr_current[row] |= mask; 
  *ptr_current[row2] |= mask; 
} /* end function add2Rows */


INLINE static void add3Rows(int row, int row2, int row3, int pivot,
			    mp_limb_t mask, /* int j_current, */
			    mp_limb_t **ptr_current) {
  // int i;
  mp_limb_t *ptr1, *ptr2, *ptr3, *ptr_piv, *ptr_lim;

  ptr1 = ptr_rows[row];
  ptr2 = ptr_rows[row2];
  ptr3 = ptr_rows[row3];
  ptr_piv = ptr_rows[pivot];
  ptr_lim = ptr_rows[pivot + 1];
  while (ptr_piv != ptr_lim) {
    mp_limb_t aux = *ptr_piv++;
    *ptr1++ ^= aux;
    *ptr2++ ^= aux;
    *ptr3++ ^= aux;
  }

  *ptr_current[row] |= mask; 
  *ptr_current[row2] |= mask; 
  *ptr_current[row3] |= mask; 
} /* end function add3Rows */

#if 0 /* unused currently */
/* add the pivot row to the given row, starting with the given column+1
   Note: the mask and ptr_current allow to do this quickly. */
/*
  TODO: do several elimination simultaneously to reduce the cost of the loop
*/

static void addPartialRows(int row, int pivot, mp_limb_t mask, int j_current,
			   mp_limb_t **ptr_current) {
  // int i;
  mp_limb_t *ptr1, *ptr2, *ptrlim;

  *ptr_current[row] ^= (*ptr_current[pivot] & mask);
  
  ptr1 = ptr_current[row] + 1;
  ptr2 = ptr_current[pivot] + 1;
  ptrlim = ptr_rows[row+1];
  while (ptr1 < ptrlim)
    *ptr1++ ^= *ptr2++;
} /* end function addPartialRows */
#endif

/* Return the index of the first row with a one in the desired column. */
/* If the column is empty, return -1 */

INLINE static int getPivot(mp_limb_t **ptr, mp_limb_t mask) {
  int i;
  for (i = 0; i < NROWS; ++i)
    if (*ptr[i]&mask)
      return i;
 
  return -1;
}


#if	VERBOSE
static void printVector(mp_limb_t *V, int n) {
  int i;
  printf("[");
  for (i = 0; i < n; ++i) {
    if (V[i / MACHINE_WORD_SIZE] & (1UL << (i % MACHINE_WORD_SIZE)))
      printf("1");
    else
      printf("0");
    if (i < n - 1)
      printf(",");
  }
  printf("]");
}

static void printMatrix(mp_limb_t *M, int nrows, int ncols,
			int limbs_per_row) {
  int i;
  printf("[\n");
  for (i = 0; i < nrows; ++i) {
    printVector(M + i*limbs_per_row, ncols);
    if (i < nrows - 1)
      printf(",\n");
  }
  printf("\n]\n");
}
#endif
