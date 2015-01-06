#ifndef MAT_INT64_H
#define MAT_INT64_H

#include <stdint.h>
#include "mat_Z.h"
#include "int64_vector.h"
#include "int64_poly.h"

typedef struct {
  int64_t **coeff;
  unsigned int NumRows;
  unsigned int NumCols;
} s_mat_int64_t;

typedef s_mat_int64_t mat_int64_t[1];
typedef s_mat_int64_t * mat_int64_ptr;
typedef const s_mat_int64_t * mat_int64_srcptr;

/*
  Initialize a matrix.

  matrix: the matrix.
  NumRows: number of rows in the matrix.
  NumCols: number of cols in the matrix.
*/
void mat_int64_init(mat_int64_ptr matrix, unsigned int NumRows, unsigned int NumCols);

/*
  Initialize a matrix and set the coefficients.

  matrix: the matrix.
  NumRows: number of rows in the matrix.
  NumCols: number of cols in the matrix.
  coeff: an array with the coefficient of the matrix, given in the order to fill
   the row.
*/
void mat_int64_init_with_array(mat_int64_ptr matrix, unsigned int NumRows, unsigned int
                           NumCols, int64_t * coeff);

/*
  Copy A in B. Assume that B and A have the same size.

  B: the new matrix.
  A: the root matrix.
 */
void mat_int64_copy(mat_int64_ptr B, mat_int64_srcptr A);

/*
  Set a coefficient to i, this coefficient is in the row (row) and the column
   (col).

  matrix: the matrix.
  i: the new coefficient.
  row: number for row.
  col: number for column.
*/
void mat_int64_set_coeff(mat_int64_ptr matrix, int64_t i, unsigned int row, unsigned int
                     col);

/*
  Set a coefficient to i, this coefficient is in the row (row) and the column
   (col).

  matrix: the matrix.
  i: the new coefficient.
  row: number for row.
  col: number for column.
*/
void mat_int64_set_coeff_int64t(mat_int64_ptr matrix, int64_t i, unsigned int row,
                            unsigned int col);

/*
  Delete a matrix.

  matrix: mat_int64_ptr, the matrix.
*/
void mat_int64_clear(mat_int64_ptr matrix);

/*
  Write a matrix in a file.

  file: the file.
  matrix: the matrix.
*/
void mat_int64_fprintf(FILE * file, mat_int64_srcptr matrix);

/*
  Transpose a matrix_src, and set the result in matrix. Assume that matrix have
   the good size and is already initialized.

  matrix: the matrix equals to the transpose of matrix_src.
  matrix_src: the root matrix.
*/
void mat_int64_transpose(mat_int64_ptr matrix, mat_int64_srcptr matrix_src);

/*
  Compute naively the product A * B = C. Assume that C have the good size.

  C: the result.
  A: a matrix.
  B: a matrix.
 */
void mat_int64_mul_mat_int64(mat_int64_ptr C, mat_int64_srcptr A, mat_int64_srcptr B);


/*
  Compute the product A * c and set the coefficient in a polynomial a.

  a: a polynomial.
  A: a matrix.
  c: a vector.
*/
void mat_int64_mul_int64_vector_to_int64_poly(int64_poly_ptr a, mat_int64_srcptr A,
					      int64_vector_srcptr c);

void mat_int64_mul_int64_vector(int64_vector_ptr a, mat_int64_srcptr A,
				int64_vector_srcptr c);

/*
  Compute the product A * c (the coefficients of c is considered as coefficients
   on a vector) and set the coefficient in a polynomial a.

  a: a polynomial.
  A: a matrix.
  c: a polynomial.
*/
void mat_int64_mul_mpz_poly(mpz_poly_ptr a, mat_int64_srcptr A, mpz_poly_srcptr c);

void mat_Z_to_mat_int64(mat_int64_ptr matrix_int, mat_Z_srcptr matrix);

#endif // MAT_INT64_H
