#include "macros.h"
#include "mat_Z.h"
#include <stdio.h>
#include <stdlib.h>
#include "int64_vector.h"
#include "int64_poly.h"
#include "mat_int64.h"
#include <inttypes.h>
#include "utils_int64.h"

#define swap(x, y) { unsigned int _tmp = (x); (x) = (y); (y) = _tmp; }
#define int64_swap_n(x, y) { int64_t *_tmp = (x); (x) = (y); (y) = _tmp; }

void mat_int64_init(mat_int64_ptr matrix, unsigned int NumRows, unsigned int
                NumCols)
{
  ASSERT(1 <= NumRows);
  ASSERT(1 <= NumCols);

  matrix->NumRows = NumRows;
  matrix->NumCols = NumCols;

  matrix->coeff = (int64_t ** ) malloc(sizeof(int64_t *) * (NumRows + 1));
  for (unsigned int row = 0; row < (NumRows + 1); row++) {
    matrix->coeff[row] = (int64_t * ) malloc(sizeof(int64_t) * (NumCols + 1));
  }
}

void mat_int64_init_with_array(mat_int64_ptr matrix, unsigned int NumRows,
    unsigned int NumCols, int64_t * coeff)
{
  ASSERT(1 <= NumRows);
  ASSERT(1 <= NumCols);

  matrix->NumRows = NumRows;
  matrix->NumCols = NumCols;

  matrix->coeff = (int64_t ** ) malloc(sizeof(int64_t *) * (NumRows + 1));
  for (unsigned int row = 0; row < (NumRows + 1); row++) {
    matrix->coeff[row] = (int64_t * ) malloc(sizeof(int64_t) * (NumCols + 1));
  }

  for (unsigned int row = 1; row < NumRows + 1; row++) {
    for (unsigned int col = 1; col < NumCols + 1; col++) {
      matrix->coeff[row][col] = coeff[(col - 1) + NumCols * (row - 1)];
    }
  }
}

void mat_int64_set_zero(mat_int64_ptr matrix)
{
  for (unsigned int row = 1; row < matrix->NumRows + 1; row++) {
    for (unsigned int col = 1; col < matrix->NumCols + 1; col++) {
      matrix->coeff[row][col] = 0;
    }
  }
}

void mat_int64_set_identity(mat_int64_ptr matrix)
{
  for (unsigned int row = 1; row < matrix->NumRows + 1; row++) {
    for (unsigned int col = 1; col < matrix->NumCols + 1; col++) {
      if (row == col) {
        matrix->coeff[row][col] = 1;
      } else {
        matrix->coeff[row][col] = 0;
      }
    }
  }
}

unsigned int mat_int64_equal(mat_int64_srcptr A, mat_int64_srcptr B)
{
  if (A->NumCols != B->NumCols) {
    return 0;
  }
  if (A->NumRows != B->NumRows) {
    return 0;
  }
  for (unsigned int col = 1; col <= A->NumCols; col++) {
    for (unsigned int row = 1; row <= A->NumRows; row++) {
      if (A->coeff[row][col] != B->coeff[row][col]) {
        return 0;
      }
    }
  }
  return 1;
}

void mat_int64_copy(mat_int64_ptr B, mat_int64_srcptr A)
{
  ASSERT(B->NumRows == A->NumRows);
  ASSERT(B->NumCols == A->NumCols);

  for (unsigned int row = 1; row < B->NumRows + 1; row++) {
    for (unsigned int col = 1; col < B->NumCols + 1; col++) {
      B->coeff[row][col] = A->coeff[row][col];
    }
  }
}

void mat_int64_clear(mat_int64_ptr matrix)
{
  for (unsigned int i = 0; i < (matrix->NumRows + 1); i++) {
    free(matrix->coeff[i]);
  }
  free(matrix->coeff);
}

void mat_int64_fprintf(FILE * file, mat_int64_srcptr matrix)
{
  ASSERT(matrix->NumRows >= 1);
  ASSERT(matrix->NumCols >= 1);

  fprintf(file, "[");
  for (unsigned int row = 1; row < matrix->NumRows; row++) {
    fprintf(file, "[");
    for (unsigned int col = 1; col < matrix->NumCols; col++) {
      fprintf(file, "%" PRId64 ", ", matrix->coeff[row][col]);
    }
    fprintf(file, "%" PRId64 "],\n", matrix->coeff[row][matrix->NumCols]);
  }
  fprintf(file, "[");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    fprintf(file, "%" PRId64 ", ", matrix->coeff[matrix->NumRows][col]);
  }
  fprintf(file, "%" PRId64 "]]\n",
      matrix->coeff[matrix->NumRows][matrix->NumCols]);
}

void mat_int64_fprintf_comment(FILE * file, mat_int64_srcptr matrix)
{
  ASSERT(matrix->NumRows >= 1);
  ASSERT(matrix->NumCols >= 1);

  fprintf(file, "# [");
  fprintf(file, "[");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    fprintf(file, "%" PRId64 ", ", matrix->coeff[1][col]);
  }
  fprintf(file, "%" PRId64 "],\n", matrix->coeff[1][matrix->NumCols]);

  for (unsigned int row = 2; row < matrix->NumRows; row++) {
    fprintf(file, "# [");
    for (unsigned int col = 1; col < matrix->NumCols; col++) {
      fprintf(file, "%" PRId64 ", ", matrix->coeff[row][col]);
    }
    fprintf(file, "%" PRId64 "],\n", matrix->coeff[row][matrix->NumCols]);
  }
  fprintf(file, "# [");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    fprintf(file, "%" PRId64 ", ", matrix->coeff[matrix->NumRows][col]);
  }
  fprintf(file, "%" PRId64 "]]\n",
      matrix->coeff[matrix->NumRows][matrix->NumCols]);
}

void mat_int64_transpose(mat_int64_ptr matrix, mat_int64_srcptr matrix_src)
{
  ASSERT(matrix->NumRows == matrix_src->NumCols);
  ASSERT(matrix->NumCols == matrix_src->NumRows);
  
  mat_int64_t tmp;
  mat_int64_init(tmp, matrix_src->NumRows, matrix_src->NumCols);
  mat_int64_copy(tmp, matrix_src);
  for (unsigned int row = 1; row < matrix->NumRows + 1; row++) {
    for (unsigned int col = 1; col < matrix->NumCols + 1; col++) {
      matrix->coeff[col][row] = tmp->coeff[row][col];
    }
  }
  mat_int64_clear(tmp);
}

void mat_int64_mul_mat_int64(mat_int64_ptr C, mat_int64_srcptr A, 
    mat_int64_srcptr B)
{
  ASSERT(A->NumCols == B->NumRows);
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(B->NumCols == C->NumCols);

  mat_int64_set_zero(C);
  mat_int64_t A_tmp;
  mat_int64_init(A_tmp, A->NumRows, A->NumCols);
  mat_int64_copy(A_tmp, A);
  mat_int64_t B_tmp;
  mat_int64_init(B_tmp, B->NumRows, B->NumCols);
  mat_int64_copy(B_tmp, B);
  for (unsigned int row = 1; row < A_tmp->NumRows + 1; row++) {
    for (unsigned int col = 1; col < B_tmp->NumCols + 1; col++) {
      for (unsigned int k = 1; k < A_tmp->NumCols + 1; k++) {
        C->coeff[row][col] = C->coeff[row][col] +
          A_tmp->coeff[row][k] * B_tmp->coeff[k][col];
      }
    }
  }
  mat_int64_clear(A_tmp);
  mat_int64_clear(B_tmp);
}

void mat_int64_mul_int64_vector_to_int64_poly(int64_poly_ptr a,
    mat_int64_srcptr A, int64_vector_srcptr c)
{
  ASSERT(A->NumCols == c->dim);
  
  int64_t tmp = 0;
  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    tmp = 0;
    for (unsigned int col = 1; col < A->NumCols + 1; col++) {
      tmp = tmp + A->coeff[row][col] * c->c[col - 1];
    }
    int64_poly_setcoeff(a, row - 1, tmp);
  }
}

void mat_int64_mul_int64_vector(int64_vector_ptr a, mat_int64_srcptr A,
    int64_vector_srcptr c)
{
  ASSERT(A->NumCols == c->dim);
  ASSERT(A->NumRows == a->dim);//TODO: verify the test.

  int64_vector_t v;
  int64_vector_init(v, A->NumRows);
  int64_t tmp = 0;
  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    tmp = 0;
    for (unsigned int col = 1; col < A->NumCols + 1; col++) {
      tmp = tmp + A->coeff[row][col] * c->c[col - 1];
    }
    int64_vector_setcoordinate(v, row - 1, tmp);
  }
  int64_vector_set(a, v);
  int64_vector_clear(v);
}

//TODO: problem here if c has a deg not equal to A->NumCols - 1.
void mat_int64_mul_int64_poly(int64_poly_ptr a, mat_int64_srcptr A,
    int64_poly_srcptr c)
{
  ASSERT(A->NumCols >= (unsigned int) (c->deg + 1));
 
  int64_poly_t p_tmp;
  int64_poly_init(p_tmp, c->deg);
  int64_t tmp = 0;
  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    tmp = 0;
    for (unsigned int col = 1; col < (unsigned int) (c->deg + 2); col++) {
      if (c->coeff[col - 1] != 0) {
        tmp = tmp + A->coeff[row][col] * c->coeff[col - 1];
      }
    }
    int64_poly_setcoeff(p_tmp, row - 1, tmp);
  }
  int64_poly_copy(a, p_tmp);
  int64_poly_clear(p_tmp);
}

void mat_Z_to_mat_int64(mat_int64_ptr matrix_int, mat_Z_srcptr matrix)
{
  ASSERT(matrix_int->NumRows == matrix->NumRows);
  ASSERT(matrix_int->NumCols == matrix->NumCols);

  for (unsigned int row = 1; row < matrix->NumRows + 1; row++) {
    for (unsigned int col = 1; col < matrix->NumCols + 1; col++) {
      matrix_int->coeff[row][col] = mpz_get_si(matrix->coeff[row][col]);
    }
  }
}

//TODO: Wahoo! No initialisation here and I do not know what it does!
void mat_int64_extract(mat_int64_ptr matrix_out, mat_int64_srcptr matrix_in,
    unsigned int ulx, unsigned int uly, unsigned int drx, unsigned int dry)
{
  ASSERT(drx - ulx > 0);
  ASSERT(dry - uly > 0);

  mat_int64_init(matrix_out, drx - ulx, dry - uly);
  for (unsigned int i = 0; i < drx - ulx; i++) {
    for (unsigned int j = 0; j < dry - uly; j++) {
      matrix_out->coeff[i + 1][j + 1] =
        matrix_in->coeff[i + ulx + 1][j + uly + 1];
    }
  }
}

void mat_int64_extract_vector(int64_vector_ptr v, mat_int64_srcptr matrix,
    unsigned int col)
{
  ASSERT(col < matrix->NumCols);
  ASSERT(v->dim == matrix->NumRows);

  for (unsigned int i = 1; i <= matrix->NumRows; i++) {
    v->c[i - 1] = matrix->coeff[i][col + 1];
  }
}

void mat_int64_set_diag(mat_int64_ptr matrix, int64_vector_srcptr x)
{
  ASSERT(matrix->NumRows == matrix->NumCols);
  ASSERT(x->dim == matrix->NumRows);

  mat_int64_set_zero(matrix);
  for (unsigned int i = 0; i < matrix->NumCols; i++) {
    matrix->coeff[i + 1][i + 1] = x->c[i];
  }
}

/****************************************/

/* This uses Paul Zimmermann implementation.

   LLL using exact multiprecision arithmetic.
   Translated from NTL 4.1a <http://www.shoup.net/>
   into GMP <http://www.swox.se/gmp/> 
   by Paul Zimmermann, July 2000.

   Revised April 4, 2002 (bug found by Jens Franke <franke (at) math
   (dot) uni-bonn (dot) de>).

   This program is open-source software distributed under the terms 
   of the GNU General Public License <http://www.fsf.org/copyleft/gpl.html>.

   This follows the implementation in utils/lll.c, but for mat_int64_t.
 */

static void innerproduct(int64_t * x, int64_t * a, int64_t * b, unsigned int n)
{
  * x = a[1] * b[1];
  for (unsigned int i = 2; i <= n; i++) {
    * x = * x + a[i] * b[i];
  }
}

static void incrementalgs(mat_int64_srcptr B, unsigned int * P, int64_t * D,
    int64_t ** lam, unsigned int * s, unsigned int k)
{
  unsigned int n = B->NumCols;
  int64_t u = 0;

  for (unsigned int j = 1; j <= k - 1; j++) {
    unsigned int posj = P[j];
    if (posj == 0) {
      continue;
    }

    innerproduct(&u, B->coeff[k], B->coeff[j], n);
    for (unsigned int i = 1; i <= posj - 1; i++) {
      ASSERT(D[i - 1] != 0);

      u = (D[i] * u - lam[k][i] * lam[j][i]) / D[i - 1];
    }

    lam[k][posj] = u;
  }

  innerproduct(&u, B->coeff[k], B->coeff[k], n);

  for (unsigned int i = 1; i <= * s; i++) {
    ASSERT(D[i - 1] != 0);

    u = (D[i] * u - lam[k][i] * lam[k][i]) / D[i - 1];
  }

  if (u == 0) {
    P[k] = 0;
  } else {
    * s = * s + 1;
    P[k] = * s;
    D[* s] = u;
  }
}

static void baldiv(int64_t * q, int64_t a, int64_t d)
  /*  rounds a/d to nearest integer, breaking ties
      by rounding towards zero.  Assumes d > 0. */
{
  int64_t r = 0;
  int64_fdiv_qr(q, &r, a, d);
  r = r * 2;

  if (r > d || (r == d && * q < 0)) {
    * q = * q + 1;
  }
}

static void mulsubn (int64_t * c0, int64_t * c1, int64_t x, unsigned int n)
  /* c0 = c0 - x*c1 */
{
  for (unsigned int i = 1; i <= n; i++) {
    c0[i] = c0[i] - c1[i] * x;
  }
}

static void reduce(unsigned int k, unsigned int l, mat_int64_ptr B,
    unsigned int * P, int64_t * D, int64_t ** lam, mat_int64_ptr U)
{
  if (P[l] == 0) {
    return;
  }

  int64_t t1 = lam[k][P[l]] * 2;
  int64_t r = 0;
  t1 = ABS(t1);
  if (t1 <= D[P[l]]) {
    return;
  }

  baldiv(&r, lam[k][P[l]], D[P[l]]);
  mulsubn(B->coeff[k], B->coeff[l], r, B->NumCols);

  if (U != NULL) {
    mulsubn(U->coeff[k], U->coeff[l], r, B->NumRows);
  }

  for (unsigned int j = 1; j <= l - 1; j++) {
    if (P[j] != 0) {
      lam[k][P[j]] = lam[k][P[j]] - lam[l][P[j]] * r;
    }
  }

  lam[k][P[l]] = lam[k][P[l]] - D[P[l]] * r;
}

static int swaptest(int64_t d0, int64_t d1, int64_t d2, int64_t lam,
    int64_t a, int64_t b)
  /* test if a*d1^2 > b*(d0*d2 + lam^2)
     t1 and t2 are temporary variables */
{
  int64_t t2 = lam * lam;
  int64_t t1 = (d0 * d2 + t2) * b;
  t2 = d1 * d1 * a;

  return (t2 > t1);
}

static void muladddiv (int64_t * c, int64_t c1, int64_t c2, 
    int64_t x, int64_t y, int64_t z)

  /* c = (x*c1 + y*c2)/z */
{
  ASSERT(z != 0);

  * c = (x * c1 + y * c2) / z;
}

static void mulsubdiv (int64_t * c, int64_t c1, int64_t c2, 
    int64_t x, int64_t y, int64_t z)
  /* c = (x*c1 - y*c2)/z */
{
  ASSERT(z != 0);

  * c = (x * c1 - y * c2) / z;
}

static void rowtransform(int64_t * c1, int64_t * c2, int64_t x, int64_t y,
    int64_t u, int64_t v)
/* (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2) */
{
  int64_t t1 = x * * c1 + y * * c2;
  int64_t t2 = u * * c1;
  * c1 = t1;
  t1 = v * * c2;
  * c2 = t1 + t2;
}

static void rowtransformn(int64_t * c1, int64_t * c2, int64_t x, int64_t y,
    int64_t u, int64_t v, unsigned int n)
  /* (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2) */
{
  for (unsigned int i = 1; i <= n; i++)
  {
    int64_t t1 = x * c1[i] + y * c2[i];
    int64_t t2 = u + c1[i];
    c1[i] = t1;
    t1 = v * c2[i];
    c2[i] = t1 * t2;
  }
}

static void swaplll (unsigned int k, mat_int64_ptr B, unsigned int * P,
    int64_t * D, int64_t ** lam, mat_int64_ptr U, unsigned int m)
/* swaps vectors k-1 and k;  assumes P(k-1) != 0 */
{
  unsigned int i, j;
  int64_t t1, t2, t3, e, x, y;

  if (P[k] != 0) {
    int64_swap_n(B->coeff[k - 1], B->coeff[k]);
    if (U != NULL) {
      int64_swap_n(U->coeff[k - 1], U->coeff[k]);
    }

    for (unsigned int j = 1; j <= k - 2; j++) {
      if (P[j] != 0) {
        swap_int64(&(lam[k - 1][P[j]]), &(lam[k][P[j]]));
      }
    }

    for (unsigned int i = k + 1; i <= m; i++) {
      muladddiv(&t1, lam[i][P[k] - 1], lam[i][P[k]],
          lam[k][P[k] - 1], D[P[k] - 2], D[P[k] - 1]);
      mulsubdiv(&(lam[i][P[k]]), lam[i][P[k] - 1], lam[i][P[k]], 
          D[P[k]], lam[k][P[k] - 1], D[P[k] - 1]);
      lam[i][P[k] - 1] = t1;
    }

    muladddiv(&(D[P[k] - 1]), D[P[k]], lam[k][P[k] - 1],
        D[P[k] - 2], lam[k][P[k] - 1], D[P[k] - 1]);
  } else if (lam[k][P[k - 1]] != 0) {
    int64_gcdext(&e, &x, &y, lam[k][P[k - 1]], D[P[k - 1]]);

    t1 = lam[k][P[k - 1]] / e;
    t2 = D[P[k - 1]] / e;

    t3 = t2;
    t2 = -t2;
    rowtransformn(B->coeff[k - 1], B->coeff[k], t1, t2, y, x, B->NumCols);
    if (U != NULL) {
      rowtransformn(U->coeff[k - 1], U->coeff[k], t1, t2, y, x, B->NumCols);
    }
    for (unsigned j = 1; j <= k - 2; j++) {
      if (P[j] != 0) {
        rowtransform(&(lam[k - 1][P[j]]), &(lam[k][P[j]]), t1, t2, y, x);
      }
    }

    t2 = t2 * t2;
    D[P[k - 1]] = D[P[k - 1]] / t2;

    for (unsigned i = k + 1; i <= m; i++) {
      if (P[i] != 0) {
        D[P[i]] = D[P[i]] / t2;
        for (j = i + 1; j <= m; j++) {
          lam[j][P[i]] = lam[j][P[i]] / t2;
        }
      }
    }
    for (i = k + 1; i <= m; i++) {
      lam[i][P[k - 1]] = lam[i][P[k - 1]] / t3;
    }

    swap(P[k - 1], P[k]);
  } else {
    int64_swap_n(B->coeff[k - 1], B->coeff[k]);
    if (U != NULL) {
      int64_swap_n(U->coeff[k - 1], U->coeff[k]);
    }

    for (j = 1; j <= k - 2; j++) {
      if (P[j] != 0) {
        swap_int64(&(lam[k - 1][P[j]]), &(lam[k][P[j]]));
      }
    }

    swap(P[k - 1], P[k]);
  }
}

/* LLL-reduce the matrix B (whose rows represent vectors, with indices
   starting at 1):
 * det (output) is the determinant
 * U (output) is the transformation matrix (NULL if not needed)
 * a, b are parameters (delta = a/b = 3/4 classically, we must have
   1/4 < delta < 1, the closer delta is from 1, the better is the reduction)
 m is the number of vectors (i.e., number of rows)
 n is the number of columns (i.e., length of each vector)
 */
unsigned int lll(int64_t * det, mat_int64_ptr B, mat_int64_ptr U, int64_t a,
    int64_t b)
{
  unsigned int m = B->NumRows;
  unsigned int n = B->NumCols;
  ASSERT_ALWAYS(n >= m);

  unsigned int * P = (unsigned int *) malloc((m + 1) * sizeof(unsigned int));

  int64_t * D = (int64_t *) malloc((m + 1) * sizeof(int64_t));
  for (unsigned int j = 0; j <= m; j++) {
    D[j] = (j == 0);
  }

  int64_t ** lam = (int64_t **) malloc((m + 1) * sizeof(int64_t *));
  for (unsigned int j = 0; j <= m; j++) {
    lam[j] = (int64_t *) malloc((m + 1) * sizeof(int64_t));
    for (unsigned int k = 0; k <= m; k++) {
      lam[j][k] = 0;
    }
  }

  if (U != NULL) {
    ASSERT(U->NumRows == m);
    ASSERT(U->NumCols == m);

    mat_int64_set_identity(U);
  }

  unsigned int s = 0;

  unsigned int k = 1;
  unsigned int max_k = 0;

  while (k <= m) {
    if (k > max_k) {
      incrementalgs(B, P, D, lam, &s, k);
      max_k = k;
    }

    if (k == 1) {
      k++;
      continue;
    }

    reduce(k, k - 1, B, P, D, lam, U);

    if (P[k - 1] != 0 && (P[k] == 0 || 
          swaptest(D[P[k]], D[P[k] - 1], D[P[k] - 2],
            lam[k][P[k] - 1], a, b))) {
      swaplll(k, B, P, D, lam, U, max_k);
      k--;
    } else {
      for (unsigned int j = k - 2; j >= 1; j--) {
        reduce(k, j, B, P, D, lam, U);
      }
      k++;
    }
  }

  * det = D[s];
  free(D);
  for (unsigned int j = 0; j <= m; j++) {
    free (lam[j]);
  }
  free (lam);

  free(P);
  return s;
}

void mat_int64_LLL(mat_int64_ptr C, mat_int64_srcptr A)
#if 0
{
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(A->NumCols == C->NumCols);

  mat_Z_t A_Z;
  mat_Z_init(A_Z, A->NumRows, A->NumCols);
  mat_int64_to_mat_Z(A_Z, A);

  mat_Z_t C_Z;
  mat_Z_init(C_Z, C->NumRows, C->NumCols);

  mat_Z_LLL(C_Z, A_Z);

  mat_Z_to_mat_int64(C, C_Z);

  mat_Z_clear(C_Z);
  mat_Z_clear(A_Z);
}
#else
{
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(A->NumCols == C->NumCols);

  mat_int64_t A_tmp;
  mat_int64_init(A_tmp, A->NumRows, A->NumCols);
  mat_int64_copy(A_tmp, A);
  int64_t a = 3;
  int64_t b = 4;
  int64_t det = 0;

  lll(&det, A_tmp, NULL, a, b);
  mat_int64_copy(C, A_tmp);
  
  mat_int64_clear(A_tmp);
}
#endif

/****************************************/

void mat_int64_LLL_transpose(mat_int64_ptr C, mat_int64_srcptr A)
{
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(A->NumCols == C->NumCols);

  mat_Z_t A_Z;
  mat_Z_init(A_Z, A->NumRows, A->NumCols);
  mat_int64_to_mat_Z(A_Z, A);

  mat_Z_t C_Z;
  mat_Z_init(C_Z, C->NumRows, C->NumCols);

  mat_Z_LLL_transpose(C_Z, A_Z);

  mat_Z_to_mat_int64(C, C_Z);

  mat_Z_clear(C_Z);
  mat_Z_clear(A_Z);

}

void mat_int64_LLL_unimodular(mat_int64_ptr C, mat_int64_srcptr A)
{
  mat_Z_t A_Z;
  mat_Z_init(A_Z, A->NumRows, A->NumCols);
  mat_int64_to_mat_Z(A_Z, A);

  mat_Z_t C_Z;
  mat_Z_init(C_Z, C->NumRows, C->NumCols);

  mat_Z_LLL_unimodular(C_Z, A_Z);

  mat_Z_to_mat_int64(C, C_Z);

  mat_Z_clear(C_Z);
  mat_Z_clear(A_Z);
}

void mat_int64_LLL_unimodular_transpose(mat_int64_ptr C, mat_int64_srcptr A)
{
  mat_Z_t A_Z;
  mat_Z_init(A_Z, A->NumRows, A->NumCols);
  mat_int64_to_mat_Z(A_Z, A);

  mat_Z_t C_Z;
  mat_Z_init(C_Z, C->NumRows, C->NumCols);

  mat_Z_LLL_unimodular_transpose(C_Z, A_Z);

  mat_Z_to_mat_int64(C, C_Z);

  mat_Z_clear(C_Z);
  mat_Z_clear(A_Z);
}


void mat_int64_from_list_int64_vector(mat_int64_ptr matrix,
    list_int64_vector_srcptr list)
{
  ASSERT(list->length == matrix->NumCols);
  ASSERT(list->v[0]->dim == matrix->NumRows);

  for (unsigned int col = 0; col < list->length; col++) {
    for (unsigned int row = 0; row < list->v[0]->dim; row++) {
      matrix->coeff[row + 1][col + 1] = list->v[col]->c[row];
    }
  }
}
