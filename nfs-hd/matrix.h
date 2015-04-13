#ifndef MATRIX_H
#define MATRIX_H 

typedef struct {
  mpz_t **coeff;
  unsigned int NumRows;
  unsigned int NumCols;
} s_mat_Z_t;

typedef s_mat_Z_t mat_Z_t[1];
typedef s_mat_Z_t * mat_Z_ptr;
typedef const s_mat_Z_t * mat_Z_srcptr;

typedef struct {
  int64_t **coeff;
  unsigned int NumRows;
  unsigned int NumCols;
} s_mat_int64_t;

typedef s_mat_int64_t mat_int64_t[1];
typedef s_mat_int64_t * mat_int64_ptr;
typedef const s_mat_int64_t * mat_int64_srcptr;

#endif /* MATRIX_H */
