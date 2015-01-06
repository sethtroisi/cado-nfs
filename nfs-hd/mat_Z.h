#ifndef MAT_Z_H
#define MAT_Z_H

#include <gmp.h>
#include <stdint.h>
#include "cado.h"
#include "utils.h"

typedef struct {
  mpz_t **coeff;
  unsigned int NumRows;
  unsigned int NumCols;
} s_mat_Z_t;

typedef s_mat_Z_t mat_Z_t[1];
typedef s_mat_Z_t * mat_Z_ptr;
typedef const s_mat_Z_t * mat_Z_srcptr;

/*
  WARNING: to correspond to the LLL algorithm in cado-nfs, all the matrix have
  the index of the rows and the columns begin at 1.
*/

/*
  Initialise a matrix.

  matrix: the matrix.
  NumRows: number of rows in the matrix.
  NumCols: number of cols in the matrix.
*/
void mat_Z_init(mat_Z_ptr matrix, unsigned int NumRows, unsigned int NumCols);

/*
  Initialise a matrix and set the coefficients of the array.

  matrix: the matrix.
  NumRows: number of rows in the matrix.
  NumCols: number of cols in the matrix.
  coeff: an array with the coefficient of the matrix, given in the order to fill
   the row.
*/
void mat_Z_init_with_array(mat_Z_ptr matrix, unsigned int NumRows, unsigned int
                           NumCols, mpz_t * coeff);

/*
  Copy A in B. Assume that B and A have the same size.

  B: the new matrix.
  A: the root matrix.
 */
void mat_Z_copy(mat_Z_ptr B, mat_Z_srcptr A);

/*
  Set a coefficient to i, this coefficient is in the row (row) and the column
   (col).

  matrix: the matrix.
  i: the new coefficient.
  row: number for row.
  col: number for column.
*/
void mat_Z_set_coeff(mat_Z_ptr matrix, mpz_t i, unsigned int row, unsigned int
                     col);

/*
  Set a coefficient to i, this coefficient is in the row (row) and the column
   (col).

  matrix: the matrix.
  i: the new coefficient.
  row: number for row.
  col: number for column.
*/
void mat_Z_set_coeff_int64(mat_Z_ptr matrix, int64_t i, unsigned int row,
                            unsigned int col);

/*
  Set a coefficient to i, this coefficient is in the row (row) and the column
   (col).

  matrix: the matrix.
  i: the new coefficient.
  row: number for row.
  col: number for column.
*/
void mat_Z_set_coeff_uint64(mat_Z_ptr matrix, uint64_t i, unsigned int row,
                             unsigned int col);

/*
  Delete a matrix.

  matrix: mat_Z_ptr, the matrix.
*/
void mat_Z_clear(mat_Z_ptr matrix);

/*
  Print a matrix.

  matrix: mat_Z_ptr, the matrix.
*/
void mat_Z_printf(mat_Z_srcptr matrix);

/*
  Write a matrix in a file.

  file: the file.
  matrix: the matrix.
*/
void mat_Z_fprintf(FILE * file, mat_Z_srcptr matrix);

/*
  Transpose a matrix_src, and set the result in matrix. Assume that matrix have
   the good size and is already initialized.

  matrix: the matrix equals to the transpose of matrix_src.
  matrix_src: the root matrix.
*/
void mat_Z_transpose(mat_Z_ptr matrix, mat_Z_srcptr matrix_src);

/*
  Compute naively the product A * B = C. Assume that C have the good size.

  C: the result.
  A: a matrix.
  B: a matrix.
 */
void mat_Z_mul_mat_Z(mat_Z_ptr C, mat_Z_srcptr A, mat_Z_srcptr B);


/*
  Compute the product A * c and set the coefficient in a polynomial a.

  a: a polynomial.
  A: a matrix.
  c: a vector.
*/
void mat_Z_mul_mpz_vector_to_mpz_poly(mpz_poly_ptr a, mat_Z_srcptr A,
                                      mpz_vector_srcptr c);

/*
  Compute the product A * c (the coefficients of c is considered as coefficients
   of a vector) and set the coefficient in a polynomial a.

  a: a polynomial.
  A: a matrix.
  c: a polynomial.
*/
void mat_Z_mul_mpz_poly(mpz_poly_ptr a, mat_Z_srcptr A, mpz_poly_srcptr c);

/*
  LLL on a matrix whose columns represent vectors.

  matrix: perform LLL on the matrix and set the result in this matrix.
*/
void mat_Z_LLL_transpose(mat_Z_ptr matrix);

/* Come from cado-nfs. */

/* LLL-reduce the matrix A (whose rows represent vectors, with indices
   starting at 1) and put the LLL reduction of A in B:
   * det (output) is the determinant
   * U (output) is the transformation matrix (NULL if not needed)
   * a, b are parameters (a/b=3/4 classically, we should have 1/4 < a/b < 1)
   m is the number of vectors (i.e., number of rows)
   n is the number of columns (i.e., length of each vector)
*/
long full_mat_Z_LLL(mat_Z_ptr B, mat_Z_ptr U, mpz_t det, mat_Z_srcptr A, mpz_t a, mpz_t b);

void mat_Z_LLL(mat_Z_ptr B, mat_Z_srcptr A);

#endif  /* MAT_Z_H */
