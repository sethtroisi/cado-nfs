//#define USE_MARKOWITZ // says it...!

/* rows correspond to relations, and columns to primes (or prime ideals) */
typedef struct {
  int nrows;
  int ncols;
  int rem_nrows;     /* number of remaining rows */
  int rem_ncols;     /* number of remaining columns */
  INT jmin, jmax;    /* we are interested in columns [jmin..jmax[ */
  INT **rows;
  int *wt;           /* weight w of column j, if w <= cwmax,
                        else <= 1 for a deleted column
                        (trick: we store -w if w > cwmax) */
  unsigned long *ad;
  unsigned long weight;
  int cwmax;         /* bound on weight of j to enter the SWAR structure */
  int rwmax;         /* if a weight(row) > rwmax, kill that row */
  int delta;         /* bound for nrows-ncols */
  int mergelevelmax; /* says it */
  INT **R;           /* R[j][k] contains the rows of the non-empty elements
                        of column j, 0 <= j < ncols, 1 <= k <= R[j][0], for
                        weight(j) <= cwmax.
                        R[j][k] = -1 if the corresponding row has been deleted.
                        R[j]=NULL for weight(j) > cwmax. */
#ifndef USE_MARKOWITZ
  dclist *S;         /* S[w] is a doubly-chained list with columns of weight w,
                        for w <= cwmax. For technical reasons (to avoid NULL
                        arguments that we can't modify), S[w] always contains
                        a first cell with value -1, thus the real list is
                        S[w]->next. */
  dclist *A;         /* A[j] points to the unique cell in S[w] containing j,
                        where w = weight(j), for w <= cwmax. A[j]=NULL for
                        w > cwmax. */
#else
  INT *MKZQ;         /* priority queue for Markowitz stuff */    
  INT *MKZA;         /* MKZA[j] gives u s.t. MKZQ[u] = j */ 
#endif
} sparse_mat_t;

