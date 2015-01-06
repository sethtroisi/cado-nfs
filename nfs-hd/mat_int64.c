#include "macros.h"
#include "mat_Z.h"
#include <stdio.h>
#include <stdlib.h>
#include "int64_vector.h"
#include "int64_poly.h"
#include "mat_int64.h"
#include <inttypes.h>

void mat_int64_init(mat_int64_ptr matrix, unsigned int NumRows, unsigned int
                NumCols)
{
  ASSERT(1 <= NumRows);
  ASSERT(1 <= NumCols);

  matrix->NumRows = NumRows;
  matrix->NumCols = NumCols;

  matrix->coeff = malloc(sizeof(int64_t *) * (NumRows + 1));
  for (unsigned int row = 0; row < (NumRows + 1); row++) {
    matrix->coeff[row] = malloc(sizeof(int64_t) * (NumCols + 1));
  }

  for (unsigned int row = 0; row < NumRows + 1; row++) {
    for (unsigned int col = 0; col < NumCols + 1; col++) {
      matrix->coeff[row][col] = 0;
    }
  }
}

void mat_int64_init_with_array(mat_int64_ptr matrix, unsigned int NumRows, unsigned int
                           NumCols, int64_t * coeff)
{
  ASSERT(1 <= NumRows);
  ASSERT(1 <= NumCols);

  matrix->NumRows = NumRows;
  matrix->NumCols = NumCols;

  matrix->coeff = malloc(sizeof(int64_t *) * (NumRows + 1));
  for (unsigned int row = 0; row < (NumRows + 1); row++) {
    matrix->coeff[row] = malloc(sizeof(int64_t) * (NumCols + 1));
  }

  for (unsigned int row = 0; row < NumRows + 1; row++) {
    for (unsigned int col = 0; col < NumCols + 1; col++) {
      matrix->coeff[row][col] = 0;
    }
  }

  for (unsigned int row = 1; row < NumRows + 1; row++) {
    for (unsigned int col = 1; col < NumCols + 1; col++) {
      //Test how the matrix is set.
      matrix->coeff[row][col] = coeff[(col - 1) + NumCols * (row - 1)];
    }
  }
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

void mat_int64_set_coeff(mat_int64_ptr matrix, int64_t i, unsigned int row, unsigned int
                         col)
{
  ASSERT(row <= matrix->NumRows);
  ASSERT(col <= matrix->NumCols);
  ASSERT(1 <= row);
  ASSERT(1 <= col);
  matrix->coeff[row][col] = i;
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
  fprintf(file, "[");
  for (unsigned int row = 1; row < matrix->NumRows; row++) {
    fprintf(file, "[");
    for (unsigned int col = 1; col < matrix->NumCols; col++) {
      fprintf(file, "%" PRId64 ", ", matrix->coeff[row][col]);
    }
    gmp_fprintf(file, "%" PRId64 "],\n", matrix->coeff[row][matrix->NumCols]);
  }
  fprintf(file, "[");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    gmp_fprintf(file, "%" PRId64 ", ", matrix->coeff[matrix->NumRows][col]);
  }
  gmp_fprintf(file, "%" PRId64 "]]\n", matrix->coeff[matrix->NumRows][matrix->NumCols]);
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

void mat_int64_mul_mat_int64(mat_int64_ptr C, mat_int64_srcptr A, mat_int64_srcptr B)
{
  ASSERT(A->NumCols == B->NumRows);
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(B->NumCols == C->NumCols);

  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    for (unsigned int col = 1; col < B->NumCols + 1; col++) {
      for (unsigned int k = 1; k < A->NumCols + 1; k++) {
        C->coeff[row][col] = C->coeff[row][col] + A->coeff[row][k] * B->coeff[k][col];
      }
    }
  }
}

void mat_int64_mul_int64_vector_to_int64_poly(int64_poly_ptr a, mat_int64_srcptr A,
                                              int64_vector_srcptr c)
{
  ASSERT(A->NumCols == c->dim);
  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    int64_t tmp = 0;
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
  int64_vector_t v;
  int64_vector_init(v, A->NumRows);
  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    int64_t tmp = 0;
    for (unsigned int col = 1; col < A->NumCols + 1; col++) {
      tmp = tmp + A->coeff[row][col] * c->c[col - 1];
    }
    int64_vector_setcoordinate(v, row - 1, tmp);
  }
  int64_vector_set(a, v);
}

void mat_int64_mul_int64_poly(int64_poly_ptr a, mat_int64_srcptr A, int64_poly_srcptr c)
{
  ASSERT(A->NumCols >= (unsigned int) (c->deg + 1));
  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    int64_t tmp = 0;
    for (unsigned int col = 1; col < (unsigned int) (c->deg + 2); col++) {
      tmp = tmp + A->coeff[row][col] * c->coeff[col - 1];
    }
    int64_poly_setcoeff(a, row - 1, tmp);
  }
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

