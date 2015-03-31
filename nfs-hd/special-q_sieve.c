#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include "utils.h"
#include "makefb.h"
#include "int64_vector.h"
#include "int64_poly.h"
#include "mat_int64.h"
#include "array.h"
#include "uint64_array.h"
#include "utils_mpz.h"
#include "utils_int64.h"
#include "utils_lattice.h"

//TODO: change % in test + add or substract.

/*
 * p = 10822639589
 * n = 6
 * f0 = 4598*x^6 + 4598*x^5 + 4597*x^4 + 4598*x^3 + 4595*x^2 + 4598*x + 4596
 * 4596,4598,4585,4598,4597,4598,4598
 * f1 = -2353645*x^6 - 2353645*x^5 + 638*x^4 - 2353645*x^3 + 4709204*x^2 -
 *   2353645*x + 2354921
 * 2354921,-2353645,4709204,-2353645,638,-2353645,-2353645
 * t = 4
 * lpb = 28
 * N_a = 2^40
 */

/*
 * Different mode:
   - TRACE_POS: follow all the operation make in a case of the array;
   - NUMBER_HIT: count the number of hit;
   - MEAN_NORM_BOUND: mean of the initialized norms;
   - MEAN_NORM_BOUND_SIEVE: mean of the log2 of the residual norms;
   - MEAN_NORM: mean of the norm we can expected if the initialisation compute
      the true norms;
   - …
 */

#ifdef TRACE_POS
FILE * file;
#endif // TRACE_POS

#ifdef NUMBER_HIT
uint64_t number_of_hit = 0;
#endif // NUMBER_HIT

#ifdef MEAN_NORM_BOUND
double norm_bound = 0;
#endif // MEAN_NORM_BOUND

#ifdef MEAN_NORM_BOUND_SIEVE
double norm_bound_sieve = 0;
#endif // MEAN_NORM_BOUND_SIEVE

#ifdef MEAN_NORM
double norm = 0;
#endif // MEAN_NORM

#ifdef ASSERT_SIEVE
uint64_t index_old = 0;
#endif // ASSERT_SIEVE

#ifdef NUMBER_SURVIVALS
uint64_t number_survivals = 0;
uint64_t number_survivals_facto = 0;
#endif // NUMBER_SURVIVALS

/*
 * Build the matrix Mq for an ideal (q, g) with deg(g) = 1. The matrix need Tq
 *  with the good coefficients, i.e. t the coeffecients not reduced modulo q.
 *  Suppose that the matrix is good initialized. Can not verify if ideal->Tq
 *  and matrix has the good size.
 *
 * matrix: the matrix with the good initialization for the number of rows and
 * columns. We can not assert here if the numbers of values for Tq and matrix is
 *  good or not.
 * ideal: the ideal (q, g) with deg(g) = 1.
 */
void build_Mq_ideal_1(mat_Z_ptr matrix, ideal_1_srcptr ideal)
{
  ASSERT(matrix->NumRows == matrix->NumCols);

  //Reinitialize the coefficients.
  for (unsigned int row = 1; row <= matrix->NumRows; row++) {
    for (unsigned int col = 1; col <= matrix->NumCols; col++) {
      mat_Z_set_coeff_uint64(matrix, 0, row, col);
    }
  }
  //Coeffecients of the diagonal.
  mat_Z_set_coeff_uint64(matrix, ideal->ideal->r, 1, 1);
  for (unsigned int row = 2; row <= matrix->NumRows; row++) {
    mat_Z_set_coeff_uint64(matrix, 1, row, row);
  }
  //Coefficients of the first line.
  for (unsigned int col = 2; col <= matrix->NumCols; col++) {
    mat_Z_set_coeff(matrix, ideal->Tr[col - 2], 1, col);
  }
}

/*
 * Build the matrix Mq for an ideal (q, g) with deg(g) > 1. This is probably
 *  non-sens.
 */
void build_Mq_ideal_u(mat_Z_ptr matrix, ideal_u_srcptr ideal)
{
  ASSERT(matrix->NumRows == matrix->NumCols);

  for (unsigned int row = 1; row <= matrix->NumRows; row++) {
    for (unsigned int col = 1; col <= matrix->NumCols; col++) {
      mat_Z_set_coeff_uint64(matrix, 0, row, col);
    }
  }
  for (int row = 1; row <= ideal->ideal->h->deg; row++) {
    mat_Z_set_coeff_uint64(matrix, ideal->ideal->r, row, row);
  }
  for (int row = ideal->ideal->h->deg + 1; row <= (int)matrix->NumRows;
       row++) {
    mat_Z_set_coeff_uint64(matrix, 1, row, row);
  }
  for (int col = ideal->ideal->h->deg + 1; col <= (int)matrix->NumCols;
       col++) {
    for (int row = 1; row <= ideal->ideal->h->deg + 1; row++) {
      mat_Z_set_coeff(matrix,
                      ideal->Tr[row - 1]
                      [col - ideal->ideal->h->deg + 1], row, col);
    }
  }
}

/*
 * Precompute some part of the computation of the bound on the resultant between
 *  f and ha. We precompute (deg(f) + 1) ^ (i/2) * (i + 1) ^ (deg(f) / 2) *
 *  infinite_norm(f), with i the degree of ha.
 *
 * pre_compute: an array in which we store the precompute elements. Need to be
 *  alloc for #(dimension of the lattice) coefficients.
 * f: a polynomial that defines the side.
 * d: degree of ha.
 */
void pre_computation(double * pre_compute, mpz_poly_srcptr f, unsigned int d)
{
  pre_compute[0] = 1;
  mpz_t infinite_norm_f_tmp;
  mpz_init(infinite_norm_f_tmp);
  //Infinite norm of f.
  mpz_poly_infinite_norm(infinite_norm_f_tmp, f);
  double infinite_norm_f = mpz_get_d(infinite_norm_f_tmp);
  mpz_clear(infinite_norm_f_tmp);
  //Do the computation for each i.
  for (unsigned int i = 1; i < d; i++) {
    pre_compute[i] = pow(infinite_norm_f, (double)i) *
      pow((double)(i + 1), ((double)f->deg / 2)) * pow((double)(f->deg + 1),
                                                       ((double)i / 2));
  }
}

/*
 * Compute the norm of a in the number field defined by f.
 *
 * res: the result.
 * f: a polynomial.
 * a: the polynomial for which we want to compute the norm.
 */
void norm_poly(mpz_ptr res, mpz_poly_srcptr f, mpz_poly_srcptr a)
{
  mpz_poly_resultant(res, f, a);
  mpz_abs(res, res);
}

#ifdef MEAN_NORM
/*
 * Compute the norm of a polynomial a and add an approximation of the norm to
 *  the mean of all the polynomial a. If the special-q is set in this number
 *  field, divide the norm by the value q.
 *
 * f: the polynomial that defines the number field.
 * a: the polynomial for which we compute the norm.
 * q: value of the special-q, 1.0 if no special-q.
 */
void mean_norm(mpz_poly_srcptr f, int64_poly_srcptr a, double q)
{
  mpz_poly_t b;
  mpz_poly_init(b, -1);
  mpz_t res;
  mpz_init(res);
  int64_poly_to_mpz_poly(b, a);
  norm_poly(res, f, b);
  norm = norm + mpz_get_d(res) / q;
  mpz_poly_clear(b);
  mpz_clear(res);
}
#endif // MEAN_NORM

#ifdef TRACE_POS
/*
 * Do the first step for the TRACE_POS mode, i.e. print c, a, f, the resultant
 * and the norm we store.
 *
 * i: index of the array we want to follow.
 * vector: the vector corresponding to the ith position in array.
 * a: the polynomial in the original lattice.
 * f: the polynomial that defines the number field.
 * bound_resultant: an approximation of the resultant.
 * special_q: 1 if the special-q is set in this number field, 0 otherwise.
 * array: the array in which we store the norms.
 */
void trace_pos_init(uint64_t i, int64_vector_srcptr vector, int64_poly_srcptr a,
    mpz_poly_srcptr f, double bound_resultant, int special_q,
    array_srcptr array)
{
  if (i == TRACE_POS) {
    fprintf(file, "c = ");
    int64_vector_fprintf(file, vector);
    fprintf(file, "a = ");
    int64_poly_fprintf(file, a);
    fprintf(file, "f = ");
    mpz_poly_fprintf(file, f);
    mpz_t res;
    mpz_init(res);
    mpz_poly_t tmp;
    mpz_poly_init(tmp, -1);
    int64_poly_to_mpz_poly(tmp, a);
    norm_poly(res, f, tmp);
    factor_t factor;
    gmp_factorize(factor, res);
    gmp_fprintf(file, "Resultant: %Zd\n", res);
    fprintf(file, "Factorization: ");
    factor_fprintf(file, factor);
    factor_clear(factor);
    mpz_clear(res);
    if (special_q) {
      fprintf(file, "Initialization (without special-q): %u\n",
              (unsigned char) log2(bound_resultant));
      fprintf(file, "Initialization (with special-q): %u\n", array->array[i]);
    } else {
      fprintf(file, "Initialization: %u\n", array->array[i]);
    }
    mpz_poly_clear(tmp);
  }
}
#endif // TRACE_POS

/*
 * Contains all the mode we can activate during the initialisation of the norms. *
 * special_q: 1 if there is a special_q, 0 otherwise.
 * f: the polynomial that defines the number field.
 * a: the polynomial.
 * q: the q of the special-q.
 * bound_resulant: an approximation of the resultant between a and f.
 * i: position in the array we want to follow.
 * vector: the vector corresponding to the ith position.
 * array: the array in which we store the norm.
 */
void mode_init_norm(MAYBE_UNUSED int special_q, MAYBE_UNUSED mpz_poly_srcptr f,
    MAYBE_UNUSED int64_poly_srcptr a, MAYBE_UNUSED double q,
    MAYBE_UNUSED double bound_resultant, MAYBE_UNUSED uint64_t i,
    MAYBE_UNUSED int64_vector_srcptr vector, MAYBE_UNUSED array_srcptr array)
{
#ifdef MEAN_NORM
  mean_norm(f, a, q);
#endif // MEAN_NORM

#ifdef MEAN_NORM_BOUND
  norm_bound = norm_bound + bound_resultant / q;
#endif // MEAN_NORM_BOUND

#ifdef TRACE_POS
  trace_pos_init(i, vector, a, f, bound_resultant, special_q, array);
#endif // TRACE_POS
}

/*
 * Compute the norm of a case in the array.
 *
 * array: array in which we store the norm.
 * i: index of the array in which we store the computed norm.
 * a: polynomial equal to matrix * vector.
 * pre_compute: array with precomputed value to compute upper bound of the norm.
 * H: sieving bound.
 * matrix: to compute the current value of a.
 * beg: index of the upper changed coordinate during the addition of 1 in the
 *  sieving region for vector.
 * f: function that defines the number field.
 * spq: the special-q we set.
 * special_q: if the special-q is set in this side, special-q is equal to 1.
 * vector: only for TRACE_POS mode.
 */
void init_each_case(array_ptr array, uint64_t i, int64_poly_ptr a,
    double * pre_compute, sieving_bound_srcptr H, mat_int64_srcptr matrix,
    unsigned int beg, mpz_poly_srcptr f, ideal_1_srcptr spq, int special_q,
    MAYBE_UNUSED int64_vector_srcptr vector)
{
  double bound_resultant;
  //a = a + sum(matrix[j][beg + 1], j, 1, H->t + 1) * x^beg.
  for (unsigned int j = 0; j < H->t; j++) {
    int64_poly_setcoeff(a, j, a->coeff[j] + matrix->coeff[j + 1][beg + 1]);
  }
  //Substract -(2 * H[i] - 1) when you change an other coordinate than c0.
  for (unsigned int k = 0; k < beg; k++) {
    for (unsigned int j = 0; j < H->t; j++) {
      int64_t tmp = a->coeff[j];
      tmp = tmp - matrix->coeff[j + 1][k + 1] * (2 * H->h[j] - 1);
      int64_poly_setcoeff(a, j, tmp);
    }
  }

  ASSERT(a->deg >= -1);
  //Verify if the computation is good.
#ifndef NDEBUG
  int64_poly_t a_tmp;
  int64_poly_init(a_tmp, -1);
  mat_int64_mul_int64_vector_to_int64_poly(a_tmp, matrix, vector);
  ASSERT(int64_poly_equal(a, a_tmp) == 0);
  int64_poly_clear(a_tmp);
#endif // NDEBUG

  //Compute an approximation of the norm.
  if (a->deg > 0) {
    uint64_t tmp = int64_poly_infinite_norm(a);
    bound_resultant = pow((double)tmp, (double)f->deg) * pre_compute[a->deg];
  } else {
    bound_resultant = 1;
  }
  /*
   * Store in array the value of the norm, without the contribution of the
   *  special-q if it is set in this number field.
   */
  if (special_q) {
    array->array[i] = (unsigned char)log2(bound_resultant) - spq->log;
    mode_init_norm(special_q, f, a, (double)spq->ideal->r, bound_resultant, i,
        vector, array);
  } else {
    array->array[i] = (unsigned char) log2(bound_resultant);
    mode_init_norm(special_q, f, a, 1.0, bound_resultant, i, vector, array);
  }
}

/*
 * Init the norm for a special-q (q, g) with deg(g) = 1.
 *
 * array: the array in which the norms are initialized.
 * pre_compute: the array of precomputed value obtained with pre_computation.
 * H: the sieving bound which gives the sieving region.
 * matrix: the MqLLL matrix.
 * f: the f which defines the side.
 * spq: the special-q.
 * special_q: 0 if there is no special-q in this side, else 1.
 */
void init_norm_1(array_ptr array, double * pre_compute,
    sieving_bound_srcptr H, mat_Z_srcptr matrix, mpz_poly_srcptr f,
    ideal_1_srcptr spq, int special_q)
{
  ASSERT(special_q == 0 || special_q == 1);

  int64_vector_t vector;
  int64_vector_init(vector, H->t);
  int64_vector_setcoordinate(vector, 0, -(int64_t)H->h[0] - 1);
  for (unsigned int i = 1; i < H->t - 1; i++) {
    int64_vector_setcoordinate(vector, i, -(int64_t)H->h[i]);
  }
  int64_vector_setcoordinate(vector, (int64_t)H->t - 1, 0);
  /*
   * Store the index of the upper changed coordinate during the addition of 1
   *  in the sieving region for vector.
   */
  unsigned int beg = 0;

  mat_int64_t matrix_int;
  mat_int64_init(matrix_int, matrix->NumRows, matrix->NumCols);
  mat_Z_to_mat_int64(matrix_int, matrix);

  //Comupte a from vector and MqLLL.
  int64_poly_t a;
  int64_poly_init(a, H->t - 1);
  mat_int64_mul_int64_vector_to_int64_poly(a, matrix_int, vector);

  //Reinitialise the value that store the mean of the norms.
#ifdef MEAN_NORM_BOUND
    norm_bound = 0;
#endif // MEAN_NORM_BOUND

#ifdef MEAN_NORM
  norm = 0;
#endif // MEAN_NORM

  int64_vector_setcoordinate(vector, 0, -(int64_t)H->h[0]);
  init_each_case(array, 0, a, pre_compute, H, matrix_int, beg, f, spq,
      special_q, vector);

  for (uint64_t i = 1; i < array->number_element; i++) {
    beg = int64_vector_add_one_i(vector, 0, H);
    init_each_case(array, i, a, pre_compute, H, matrix_int, beg, f, spq,
        special_q, vector);
  }

  int64_poly_clear(a);
  int64_vector_clear(vector);
  mat_int64_clear(matrix_int);
}

/*
 * Tqr is the normalised Tqr obtained by compute_Tqr_1. If Tqr[0] != 0,
 *  pseudo_Tqr = [(-Tqr[0])^-1 mod r = a, a * Tqr[1], …].
 *
 * pseudo_Tqr: a matrix (a line here) obtained as describe above.
 * Tqr: the Tqr matrix.
 * t: dimension of the lattice.
 * ideal: the ideal r.
 */
void compute_pseudo_Tqr_1(uint64_t * pseudo_Tqr, uint64_t * Tqr,
    unsigned int t, ideal_1_srcptr ideal)
{
  unsigned int i = 0;
  while (Tqr[i] == 0 && i < t) {
    pseudo_Tqr[i] = 0;
    i++;
  }
  uint64_t inverse = ideal->ideal->r - 1;

  ASSERT(Tqr[i] == 1);
  ASSERT(inverse ==
         invmod_uint64((uint64_t)(-Tqr[i] + (int64_t)ideal->ideal->r),
                       ideal->ideal->r));
  pseudo_Tqr[i] = inverse;
  for (unsigned int j = i + 1; j < t; j++) {
    pseudo_Tqr[j] = (inverse * Tqr[j]) % ((int64_t)ideal->ideal->r);
  }
}

/*
 * Mqr is the a matrix that can generate vector c such that Tqr * c = 0 mod r.
 * Mqr = [[r a 0 … 0]
 *       [0 1 b … 0]
 *       [| | | | |]
 *       [0 0 0 0 1]]
 *
 * Mqr: the Mqr matrix.
 * Tqr: the Tqr matrix.
 * t: dimension of the lattice.
 * ideal: the ideal r.
 */
void compute_Mqr_1(mat_int64_ptr Mqr, uint64_t * Tqr, unsigned int t,
    ideal_1_srcptr ideal)
{
  mat_int64_init(Mqr, t, t);
  mat_int64_set_zero(Mqr);
  Mqr->coeff[1][1] = ideal->ideal->r;
  for (unsigned int i = 2; i <= t; i++) {
    Mqr->coeff[i][i] = 1;
  }
  for (unsigned int col = 2; col <= t; col++) {
    if (Tqr[col-1] != 0) {
      Mqr->coeff[1][col] = (-(int64_t)Tqr[col-1]) + (int64_t)ideal->ideal->r;
    } else {
      Mqr->coeff[1][col] = 0;
    }
  }
}

/*
 * Compute Tqr for r an ideal of degree 1 and normalised it.
 *
 * Tqr: the Tqr matrix (Tqr is a line).
 * matrix: the MqLLL matrix.
 * t: dimension of the lattice.
 * ideal: the ideal r.
 */
void compute_Tqr_1(uint64_t * Tqr, mat_Z_srcptr matrix,
    unsigned int t, ideal_1_srcptr ideal)
{
  mpz_t tmp;
  mpz_init(tmp);
  unsigned int i = 0;
  Tqr[i] = 0;
  mpz_t invert;
  mpz_init(invert);

  //Tqr = Mq,1 - Tr * Mq,2.
  for (unsigned int j = 0; j < t; j++) {
    mpz_set(tmp, matrix->coeff[1][j + 1]);
    for (unsigned int k = 0; k < t - 1; k++) {
      mpz_submul(tmp, ideal->Tr[k], matrix->coeff[k + 2][j + 1]);
    }
    if (Tqr[i] == 0) {
      mpz_mod_ui(tmp, tmp, ideal->ideal->r);
      if (mpz_cmp_ui(tmp, 0) != 0) {
        mpz_invert_ui(tmp, tmp, ideal->ideal->r);
        mpz_set(invert, tmp);
        Tqr[j] = 1;
        i = j;
      } else {
        Tqr[j] = 0;
      }
    } else {
      mpz_mul(tmp, tmp, invert);
      mpz_mod_ui(tmp, tmp, ideal->ideal->r);
      Tqr[j] = mpz_get_ui(tmp);
    }
  }

  mpz_clear(tmp);
  mpz_clear(invert);
}

/*
 * Generate the Mqr matrix, r is an ideal of degree 1.
 *
 * Mqr: the matrix.
 * Tqr: the Tqr matrix (Tqr is a line).
 * t: dimension of the lattice.
 * ideal: the ideal r.
 */
void generate_Mqr_1(mat_int64_ptr Mqr, uint64_t * Tqr,
    unsigned int t, ideal_1_srcptr ideal)
{
  ASSERT(Tqr[0] != 0);
  ASSERT(t == Mqr->NumCols);
  ASSERT(t == Mqr->NumRows);

#ifndef NDEBUG
  for (unsigned int row = 1; row <= Mqr->NumRows; row++) {
    for (unsigned int col = 1; col <= Mqr->NumCols; col++) {
      ASSERT(Mqr->coeff[row][col] == 0);
    }
  }
#endif // NDEBUG

  Mqr->coeff[1][1] = ideal->ideal->r;
  for (unsigned int col = 2; col <= Mqr->NumRows; col++) {
    if (Tqr[col] == 0) {
      Mqr->coeff[1][col] = 0;
    } else {
      Mqr->coeff[1][col] = (-Tqr[col - 1]) + ideal->ideal->r;
    }
    ASSERT(Mqr->coeff[col][1] < (int64_t)ideal->ideal->r);
  }

  for (unsigned int row = 2; row <= Mqr->NumRows; row++) {
    Mqr->coeff[row][row] = 1;
  }
}

/*
 * WARNING: Need to assert if the function does what it is supposed to do.
 * Compute Tqr for an ideal (r, h) with deg(h) > 1.
 *
 * Tqr: Tqr is a matrix in this case.
 * matrix: MqLLL.
 * t: dimension of the lattice.
 * ideal: (r, h).
 */
void compute_Tqr_u(mpz_t ** Tqr, mat_Z_srcptr matrix, unsigned int t,
    ideal_u_srcptr ideal)
{
  /* Tqr = Mq,1 - Tr * Mq,2. */
  for (int row = 0; row < ideal->ideal->h->deg; row++) {
    for (unsigned int col = 0; col < t; col++) {
      mpz_init(Tqr[row][col]);
      mpz_set(Tqr[row][col], matrix->coeff[row + 1][col + 1]);
      for (int k = 0; k < (int)t - ideal->ideal->h->deg; k++) {
        mpz_submul(Tqr[row][col], ideal->Tr[row][k],
            matrix->coeff[k + ideal->ideal->h->deg + 1][col + 1]);
      }
      mpz_mod_ui(Tqr[row][col], Tqr[row][col], ideal->ideal->r);
    }
  }
}

#ifdef ASSERT_SIEVE
/*
 * To assert that the ideal is a factor of a = matrix * c mapped in the number
 *  field defined by f.
 *
 * H: the sieving bound.
 * index: index of c in the array in which we store the norm.
 * number_element: number of element contains in the sieving region.
 * matrix: MqLLL.
 * f: the polynomial that defines the number field.
 * ideal: the ideal involved in the factorisation in ideals.
 * c: the vector corresponding to the ith value of the array in which we store
 *  the norm.
 */
void assert_sieve(sieving_bound_srcptr H, uint64_t index,
    uint64_t number_element, mat_Z_srcptr matrix, mpz_poly_srcptr f,
    ideal_1_srcptr ideal, int64_vector_srcptr c)
{
  /*
   * Verify if we do not sieve the same index, because we do not sieve the power
   *  of an ideal.
   */
  if (index_old == 0) {
    ASSERT(index_old <= index);
  } else {
    ASSERT(index_old < index);
  }
  index_old = index;

  mpz_vector_t v;
  mpz_vector_init(v, H->t);
  mpz_poly_t a;
  mpz_poly_init(a, -1);
  mpz_t res;
  mpz_init(res);

  array_index_mpz_vector(v, index, H, number_element);
  mat_Z_mul_mpz_vector_to_mpz_poly(a, matrix, v);
  norm_poly(res, a, f);
  //Verify if res = r * …
  ASSERT(mpz_congruent_ui_p(res, 0, ideal->ideal->r));

  //Verify c and v are the same.
  if (c != NULL) {
    mpz_vector_t v_tmp;
    mpz_vector_init(v_tmp, c->dim);
    int64_vector_to_mpz_vector(v_tmp, c);
    ASSERT(mpz_vector_cmp(v, v_tmp) == 0);
    mpz_vector_clear(v_tmp);
  }

  //Verify if h divides a mod r.
  if (a->deg != -1) {
    mpz_poly_t r;
    mpz_poly_init(r, -1);
    mpz_poly_t q;
    mpz_poly_init(q, -1);
    mpz_t p;
    mpz_init(p);
    mpz_set_uint64(p, ideal->ideal->r);

    ASSERT(mpz_poly_div_qr(q, r, a, ideal->ideal->h, p));
    ASSERT(r->deg == -1);
    mpz_poly_clear(r);
    mpz_poly_clear(q);
    mpz_clear(p);
  }

  mpz_clear(res);
  mpz_poly_clear(a);
  mpz_vector_clear(v);
}
#endif // ASSERT_SIEVE

#ifdef TRACE_POS
/*
 * Say all the ideal we delete in the ith position.
 *
 * index: index in array we want to follow.
 * ideal: the ideal we delete to the norm.
 * array: the array in which we store the resulting norm.
 */
void trace_pos_sieve(uint64_t index, ideal_1_srcptr ideal, array_srcptr array)
{
  if (index == TRACE_POS) {
    fprintf(file, "The ideal is: ");
    ideal_fprintf(file, ideal->ideal);
    fprintf(file, "log: %u\n", ideal->log);
    fprintf(file, "The new value of the norm is %u.\n", array->array[index]);
  }
}
#endif // TRACE_POS

/*
 * All the mode we can active during the sieving step.
 *
 * H: the sieving bound.
 * index: index of the array in which we store the resulting norm.
 * array: the array in which we store the resulting norm.
 * matrix: MqLLL.
 * f: the polynomial that defines the number field.
 * ideal: the ideal we want to delete.
 * c: the vector corresponding to indexth position in array.
 * number_c_l: number of possible c with the same ci, ci+1, …, ct. (cf line_sieve_ci)
 * nbint: say if we have already count the number_c_l or not.
 */
void mode_sieve(MAYBE_UNUSED sieving_bound_srcptr H,
    MAYBE_UNUSED uint64_t index, MAYBE_UNUSED array_srcptr array,
    MAYBE_UNUSED mat_Z_srcptr matrix, MAYBE_UNUSED mpz_poly_srcptr f,
    MAYBE_UNUSED ideal_1_srcptr ideal, MAYBE_UNUSED int64_vector_srcptr c,
    MAYBE_UNUSED uint64_t number_c_l, MAYBE_UNUSED unsigned int nbint)
{
#ifdef ASSERT_SIEVE
  assert_sieve(H, index, array->number_element, matrix, f, ideal, c);
#endif // ASSERT_SIEVE

#ifdef NUMBER_HIT
  if (nbint) {
    number_of_hit = number_of_hit + number_c_l;
  }
#endif // NUMBER_HIT

#ifdef TRACE_POS
  trace_pos_sieve(index, ideal, array);
#endif // TRACE_POS
}

/*
 * Sieve for a special-q of degree 1 and Tqr with a zero coefficient at the
 *  first place. This function sieve a c in the q lattice.
 *
 * array: in which we store the norms.
 * c: element of the q lattice.
 * ideal: an ideal with r < q.
 * ci: the possible first coordinate of c to have c in the sieving region.
 * H: the sieving bound.
 * i: index of the first non-zero coefficient in pseudo_Tqr.
 * number_c_l: number of possible c with the same ci, ci+1, …, ct.
 */
void line_sieve_ci(array_ptr array, int64_vector_ptr c, ideal_1_srcptr ideal,
    int64_t ci, sieving_bound_srcptr H, unsigned int i, uint64_t number_c_l,
    MAYBE_UNUSED mat_Z_srcptr matrix, MAYBE_UNUSED mpz_poly_srcptr f)
{
#ifndef NDEBUG
  if (i == 0) {
    ASSERT(number_c_l == 1);
  } else {
    ASSERT(number_c_l % 2 == 0);
  }
#endif // NDEBUG

  uint64_t index = 0;

#ifdef ASSERT_SIEVE
  index_old = 0;
#endif // ASSERT_SIEVE

  if (ci < (int64_t)H->h[i]) {
    //Change the ith coordinate of c.
    int64_vector_setcoordinate(c, i, ci);
    index = array_int64_vector_index(c, H, array->number_element);
    array->array[index] = array->array[index] - ideal->log;

    mode_sieve(H, index, array, matrix, f, ideal, c, number_c_l, 1);

    for (uint64_t k = 1; k < number_c_l; k++) {
      array->array[index + k] = array->array[index + k] - ideal->log;
    
      mode_sieve(H, index + k, array, matrix, f, ideal, NULL, number_c_l, 0);
    }

    int64_t tmp = ci;
    tmp = tmp + (int64_t)ideal->ideal->r;

    while(tmp < (int64_t)H->h[i]) {
      index = index + ideal->ideal->r * number_c_l;

      array->array[index] = array->array[index] - ideal->log;
 
      //TODO: nbint = 0 is strange.
      mode_sieve(H, index, array, matrix, f, ideal, NULL, number_c_l, 0);

      for (uint64_t k = 1; k < number_c_l; k++) {
        array->array[index + k] = array->array[index + k] - ideal->log;
 
        mode_sieve(H, index + k, array, matrix, f, ideal, NULL, number_c_l, 0);
      }

      tmp = tmp + (int64_t)ideal->ideal->r;
    }
  }
}

/*
 *
 *
 * array: the array in which we store the resulting norms.
 * c: the 
 */
void line_sieve_1(array_ptr array, int64_vector_ptr c, uint64_t * pseudo_Tqr,
    ideal_1_srcptr ideal, sieving_bound_srcptr H, unsigned int i,
    uint64_t number_c_l, int64_t * ci, unsigned int pos,
    MAYBE_UNUSED mat_Z_srcptr matrix, MAYBE_UNUSED mpz_poly_srcptr f)
{
  ASSERT(pos >= i);

  //Recompute ci.
  for (unsigned int j = i + 1; j < pos; j++) {
    * ci = * ci - ((int64_t)pseudo_Tqr[j] * (2 * (int64_t)H->h[j] - 1));
    if (* ci >= (int64_t)ideal->ideal->r || * ci < 0) {
      * ci = * ci % (int64_t)ideal->ideal->r;
    }
  }
  * ci = * ci + pseudo_Tqr[pos];
  if (* ci >= (int64_t)ideal->ideal->r) {
    * ci = * ci - (int64_t)ideal->ideal->r;
  }
  if (* ci < 0) {
    * ci = * ci + (int64_t)ideal->ideal->r;
  }
  ASSERT(* ci >= 0 && * ci < (int64_t)ideal->ideal->r);

  //Compute the lowest value such that Tqr*(c[0], …, c[i], …, c[t-1]) = 0 mod r.
  //lb: the minimal accessible value in the sieving region for this coordinate.
  int64_t lb = 0;
  if (i < H->t - 1) {
    lb =  -(int64_t)H->h[i];
  }
  int64_t k = (lb - * ci) / (-(int64_t) ideal->ideal->r);
  if ((lb - * ci) > 0) {
    k--;
  }
  int64_t ci_tmp = -k * (int64_t)ideal->ideal->r + * ci;
  ASSERT(ci_tmp >= lb);

#ifndef NDEBUG
  int64_t tmp_ci = ci_tmp;
  tmp_ci = tmp_ci % (int64_t) ideal->ideal->r;
  if (tmp_ci < 0) {
    tmp_ci = tmp_ci + (int64_t) ideal->ideal->r;
  }
  int64_t tmp = 0;

  for (unsigned int j = i + 1; j < H->t; j++) {
    tmp = tmp + ((int64_t)pseudo_Tqr[j] * c->c[j]);
    tmp = tmp % (int64_t) ideal->ideal->r;
    if (tmp < 0) {
      tmp = tmp + (int64_t) ideal->ideal->r;
    }
  }
  ASSERT(tmp == tmp_ci);
#endif // NDEBUG

  line_sieve_ci(array, c, ideal, ci_tmp, H, i, number_c_l, matrix, f);
}

/*
 * Plane sieve.
 *
 * array: the array in which we store the norms.
 * r: the ideal we consider.
 * Mqr: the Mqr matrix.
 * H: the sieving bound that defines the sieving region.
 */
void plane_sieve_1(array_ptr array, ideal_1_srcptr r,
    mat_int64_srcptr Mqr, sieving_bound_srcptr H,
    MAYBE_UNUSED mat_Z_srcptr matrix, MAYBE_UNUSED mpz_poly_srcptr f)
{
  ASSERT(Mqr->NumRows == Mqr->NumCols);
  ASSERT(Mqr->NumRows == 3);

  //Perform the Franke-Kleinjung algorithm.
  int64_vector_t * vec = malloc(sizeof(int64_vector_t) * Mqr->NumRows);
  for (unsigned int i = 0; i < Mqr->NumCols; i++) {
    int64_vector_init(vec[i], Mqr->NumRows);
    mat_int64_extract_vector(vec[i], Mqr, i);
  }
  int64_vector_t e0;
  int64_vector_t e1;
  int64_vector_init(e0, vec[0]->dim);
  int64_vector_init(e1, vec[1]->dim);
  int boolean = reduce_qlattice(e0, e1, vec[0], vec[1], 2*H->h[0]);
  
  //Find some short vectors to go from z = d to z = d + 1.
  list_vector_t SV;
  list_vector_init(SV);
  SV4(SV, vec[0], vec[1], vec[2]);

  //Reduce q-lattice is not possible.
  if (boolean == 0) {
    mat_int64_fprintf_comment(stdout, Mqr);
 
    //plane_sieve_whithout_FK(SV, vec);

    int64_vector_clear(e0);  
    int64_vector_clear(e1);
    for (unsigned int i = 0; i < Mqr->NumCols; i++) {
      int64_vector_clear(vec[i]);
    }
    free(vec);
    list_vector_clear(SV);

    return;
  }

  for (unsigned int i = 0; i < Mqr->NumCols; i++) {
    int64_vector_clear(vec[i]);
  }
  free(vec);

  int64_vector_t vs;
  int64_vector_init(vs, e0->dim);
  int64_vector_set_zero(vs);
  int64_vector_t v;
  int64_vector_init(v, vs->dim);

  int64_vector_t vs_tmp;
  int64_vector_init(vs_tmp, vs->dim);

  uint64_t coord_e0 = 0, coord_e1 = 0;
  coordinate_FK_vector(&coord_e0, &coord_e1, e0, e1, H, array->number_element);

  //Enumerate the element of the sieving region.
  for (unsigned int d = 0; d < H->h[2]; d++) {
    int64_vector_set(v, vs);
    int64_vector_set(vs_tmp, vs);

    unsigned int FK_value = 0;
    //Say if we have a vector in the sieving region.
    unsigned int flag_sr = 0;
    uint64_t index_v = 0;

    //Perform Franke-Kleinjung enumeration.
    while (v->c[1] < (int64_t)H->h[1]) {
      if (v->c[1] >= -(int64_t)H->h[1]) {
        if (!flag_sr) {
          index_v = array_int64_vector_index(v, H, array->number_element);
          flag_sr = 1;
        } else {
          if (FK_value == 0) {
            index_v = index_v + coord_e0;
          } else if (FK_value == 1) {
            index_v = index_v + coord_e1;
          } else {
            ASSERT(FK_value == 2);
            index_v = index_v + coord_e0 + coord_e1;
          }

          mode_sieve(H, index_v, array, matrix, f, r, v, 1, 1);

        }
        array->array[index_v] = array->array[index_v] - r->log;
      }
      if (ABS(v->c[0]) < ABS(vs_tmp->c[0])) {
        int64_vector_set(vs_tmp, v);
      }
      FK_value = enum_pos_with_FK(v, v, e0, e1, -(int64_t)H->h[0], 2 * (int64_t)
          H->h[0]);
    }

    flag_sr = 0;
    int64_vector_set(v, vs);
    FK_value = enum_neg_with_FK(v, v, e0, e1, -(int64_t)H->h[0] + 1, 2 *
        (int64_t)H->h[0]);
    while (v->c[1] >= -(int64_t)H->h[1]) {
      if (v->c[1] < (int64_t)H->h[1]) {
        if (!flag_sr) {
          index_v = array_int64_vector_index(v, H, array->number_element);
          flag_sr = 1;
        } else {
          if (FK_value == 0) {
            index_v = index_v - coord_e0;
          } else if (FK_value == 1) {
            index_v = index_v - coord_e1;
          } else {
            ASSERT(FK_value == 2);
            index_v = index_v - coord_e0 - coord_e1;
          }

          mode_sieve(H, index_v, array, matrix, f, r, v, 1, 1);

        }
        array->array[index_v] = array->array[index_v] - r->log;
      }
      if (ABS(v->c[0]) < ABS(vs_tmp->c[0])) {
        int64_vector_set(vs_tmp, v);
      }
      FK_value = enum_neg_with_FK(v, v, e0, e1, -(int64_t)H->h[0] + 1, 2 *
          (int64_t) H->h[0]);
    }
    int64_vector_set(vs, vs_tmp);

    //Jump in the next plane.
    list_vector_t list_vector_tmp;
    list_vector_init(list_vector_tmp);
    int64_vector_t v_tmp;
    int64_vector_init(v_tmp, SV->v[0]->dim);
    /*
     * 0: in the sieving region.
     * 1: in [-H_0, H_0[
     * 2: outside
     */
    unsigned char * assert = malloc(sizeof(unsigned char) * SV->length);
    unsigned char found = 2;

    for (unsigned int i = 0; i < SV->length; i++) {
      int64_vector_add(v_tmp, vs, SV->v[i]);
      list_vector_add_int64_vector(list_vector_tmp, v_tmp);

      if (-(int64_t) H->h[0] <= v_tmp->c[0] && (int64_t)H->h[0] > v_tmp->c[0]) {
        if (-(int64_t) H->h[1] <= v_tmp->c[1] && (int64_t)H->h[1] > v_tmp->c[1])
        {
          assert[i] = 0;
          found = 0;
        } else {
          assert[i] = 1;
          if (found > 1) {
            found = 1;
          }
        }
      } else {
        assert[i] = 2;
      }
    }
    int64_vector_clear(v_tmp);

    unsigned int pos = 0;
    if (found == 0) {
      int64_t min = -1;
      for (unsigned int i = 1; i < list_vector_tmp->length; i++) {
        if (assert[i] == 0) {
          int64_t tmp= ABS(vs->c[0] - list_vector_tmp->v[i]->c[0]);
          if (min == -1) {
            min = tmp;
          }
          if (min >= tmp) {
            min = tmp;
            pos = i;
          }
        }
      }
      int64_vector_set(vs, list_vector_tmp->v[pos]);
      ASSERT(vs->c[0] < (int64_t)H->h[0]);
      ASSERT(vs->c[0] >= -(int64_t)H->h[0]);
    } else if (found == 1) {
      pos = find_min_x(list_vector_tmp, H, assert, 1);
      int64_vector_set(vs, list_vector_tmp->v[pos]);
      ASSERT(vs->c[0] < (int64_t)H->h[0]);
      ASSERT(vs->c[0] >= -(int64_t)H->h[0]);
    } else {
      ASSERT(found == 2);
      pos = find_min_x(list_vector_tmp, H, assert, 2);
      int64_vector_set(vs, list_vector_tmp->v[pos]);
      add_FK_vector(vs, e0, e1, H);
      ASSERT(vs->c[0] < (int64_t)H->h[0]);
      ASSERT(vs->c[0] >= -(int64_t)H->h[0]);
    }

    free(assert);
    list_vector_clear(list_vector_tmp);
  }

  list_vector_clear(SV);
  int64_vector_clear(v);
  int64_vector_clear(vs);
  int64_vector_clear(vs_tmp);
  int64_vector_clear(e0);
  int64_vector_clear(e1);
}

#ifdef SIEVE_U
void sieve_u(array_ptr array, mpz_t ** Tqr, ideal_u_srcptr ideal,
    mpz_vector_srcptr c, sieving_bound_srcptr H)
{
  //Not optimal
  int nul = 1;
  mpz_t * tmp = (mpz_t * )
    malloc(sizeof(mpz_t) * ideal->ideal->h->deg);
  for (int row = 0; row < ideal->ideal->h->deg;
       row++) {
    mpz_init(tmp[row]);
    for (unsigned int col = 0; col < c->dim; col++) {
      mpz_addmul(tmp[row], Tqr[row][col], c->c[col]);
    }
    mpz_mod_ui(tmp[row], tmp[row], ideal->ideal->r);
    if (mpz_cmp_ui(tmp[row], 0) != 0) {
      nul = 0;
    }
  }

  if (nul == 1) {
    uint64_t index = 0;
    index = array_mpz_vector_index(c, H, array->number_element);
    array->array[index] = array->array[index] -
      ideal->log;

#ifdef TRACE_POS
    if (index == TRACE_POS) {
      fprintf(file, "The ideal is: ");
      ideal_u_fprintf(file, ideal, H->t);
      fprintf(file, "The new value of the norm is %u.\n", array->array[index]);
    }
#endif // TRACE_POS
  }

  for (int row = 0; row < ideal->ideal->h->deg;
       row++) {
    mpz_clear(tmp[row]);
  }
  free(tmp);
}
#endif // SIEVE_U

/*
 *
 *
 * array: the array in which we store the resulting norms.
 * matrix: MqLLL.
 * fb: factor base of this side.
 * H: sieving bound.
 * f: the polynomial that defines the number field.
 */
void special_q_sieve(array_ptr array, mat_Z_srcptr matrix,
    factor_base_srcptr fb, sieving_bound_srcptr H,
    MAYBE_UNUSED mpz_poly_srcptr f)
{
#ifdef TIMER_SIEVE
  double time_sieve = 0;
  double sec = 0;
#endif // TIMER_SIEVE

  //For all the ideal_1 of the factor base.
  for (uint64_t i = 0; i < fb->number_element_1; i++) {

    ideal_1_t r;
    ideal_1_init(r);
    ideal_1_set(r, fb->factor_base_1[i], H->t);

#ifdef TIMER_SIEVE
    if (r->ideal->r > TIMER_SIEVE) {
      sec = seconds();
    }
#endif // TIMER_SIEVE

    uint64_t * Tqr = (uint64_t *) malloc(sizeof(uint64_t) * (H->t));

    //Compute the true Tqr
    compute_Tqr_1(Tqr, matrix, H->t, r);

    if (r->ideal->r < 2 * (uint64_t) H->h[0]) {
      uint64_t * pseudo_Tqr = (uint64_t *) malloc(sizeof(uint64_t) * (H->t));
      //Compute a usefull pseudo_Tqr for line_sieve.
      compute_pseudo_Tqr_1(pseudo_Tqr, Tqr, H->t, r);

      //Initialise c = (-H[0], -H[1], …, 0).
      int64_vector_t c;
      int64_vector_init(c, H->t);
      for (unsigned int j = 0; j < H->t - 1; j++) {
        int64_vector_setcoordinate(c, j, -(int64_t)H->h[j]);
      }
      int64_vector_setcoordinate(c, H->t - 1, 0);

      if (pseudo_Tqr[0] != 0) {
        unsigned int pos = 0;

        //Set c = (-H[0] - 1, -H[1], …, 0).
        int64_t c0 = 0;
        int64_vector_setcoordinate(c, 1, c->c[1] - 1);
        //Compute c0 = (-Tqr[0])^(-1) * (Tqr[1]*c[1] + …) mod r.
        for (unsigned int j = 1; j < H->t; j++) {
          c0 = c0 + (int64_t)pseudo_Tqr[j] * c->c[j];
          if (c0 >= (int64_t)r->ideal->r || c0 < 0) {
            c0 = c0 % (int64_t)r->ideal->r;
          }
        }
        if (c0 < 0) {
          c0 = c0 + (int64_t)r->ideal->r;
        }
        ASSERT(c0 >= 0);
        //Number of c with the same c[1], …, c[t-1].
        uint64_t number_c = array->number_element / (2 * H->h[0]);

#ifdef TIMER_SIEVE
        if (r->ideal->r > TIMER_SIEVE) {
          sec = seconds();
        }
#endif // TIMER_SIEVE

        for (uint64_t j = 0; j < number_c; j++) {
          pos = int64_vector_add_one_i(c, 1, H);
          line_sieve_1(array, c, pseudo_Tqr, r, H, 0, 1, &c0, pos,
              matrix, f);
        }

#ifdef TIMER_SIEVE
        if (r->ideal->r > TIMER_SIEVE) {
          time_sieve = time_sieve + seconds() - sec;
        }
#endif // TIMER_SIEVE

#ifdef NUMBER_HIT
        printf("Number of hits: %" PRIu64 " for r: %" PRIu64 "\n", number_of_hit,
            r->ideal->r);
        printf("Estimated number of hits: %u.\n",
            (unsigned int) nearbyint((double) array->number_element /
              (double) r->ideal->r));
        number_of_hit = 0;
#endif // NUMBER_HIT

        free(pseudo_Tqr);
      }

      int64_vector_clear(c);
    } else if (r->ideal->r < 4 * (int64_t)(H->h[0] * H->h[1])) {

#ifdef TEST_MQR
      printf("# Tqr = [");
      for (unsigned int i = 0; i < H->t - 1; i++) {
        printf("%" PRIu64 ", ", Tqr[i]);
      }
      printf("%" PRIu64 "]\n# Mqr =\n", Tqr[H->t - 1]);
      mat_int64_t Mqr_test;
      compute_Mqr_1(Mqr_test, Tqr, H->t, r);
      mat_int64_fprintf_comment(stdout, Mqr_test);
      mat_int64_clear(Mqr_test);
#endif // TEST_MQR

      ASSERT(H->t == 3);

      mat_int64_t Mqr;
      mat_int64_init(Mqr, H->t, H->t);

      compute_Mqr_1(Mqr, Tqr, H->t, r);

#ifdef TIMER_SIEVE
      if (r->ideal->r > TIMER_SIEVE) {
        sec = seconds();
      }
#endif // TIMER_SIEVE

      plane_sieve_1(array, r, Mqr, H, matrix, f);

#ifdef TIMER_SIEVE
      if (r->ideal->r > TIMER_SIEVE) {
        time_sieve = time_sieve + seconds() - sec;
      }
#endif // TIMER_SIEVE

      mat_int64_clear(Mqr);
    }
    ideal_1_clear(r, H->t);
    free(Tqr);
  }

#ifdef TIMER_SIEVE
  printf("Bound: %d\n", TIMER_SIEVE);
  printf("Perform sieve: %f\n", time_sieve);
#endif // TIMER_SIEVE

}

#if 0
#ifdef SIEVE_TQR
    else {
      unsigned int index = 1;
      while (pseudo_Tqr[index] == 0 && index < H->t) {
        index++;
      }

      if (index == H->t - 1) {
        /* uint64_t number_c_u = array->number_element; */
        /* uint64_t number_c_l = 1; */
        /* for (uint64_t j = 0; j < index; j++) { */
        /*   number_c_u = number_c_u / (2 * H->h[j] + 1); */
        /*   number_c_l = number_c_l * (2 * H->h[j] + 1); */
        /* } */

        /* line_sieve(array, c, r, 0, H, index, number_c_l, */
        /*          matrix, f); */
      } else {
        int64_t ci = 0;
        if (index + 1 < H->t - 1) {
          int64_vector_setcoordinate(c, index + 1, c->c[index + 1] - 1);
        } else if (index + 1 == H->t - 1) {
          int64_vector_setcoordinate(c, index + 1, -1);
        }

        if (index < H->t - 1) {
          for (unsigned int j = index + 1; j < H->t; j++) {
            ci = ci + (int64_t)pseudo_Tqr[j] * c->c[j];
            if (ci >= (int64_t)r->ideal->r || ci < 0) {
              ci = ci % (int64_t)r->ideal->r;
            }
          }
          if (ci < 0) {
            ci = ci + (int64_t)r->ideal->r;
          }
          ASSERT(ci >= 0);
        }

        uint64_t number_c_u = array->number_element;
        uint64_t number_c_l = 1;
        for (uint64_t j = 0; j < index; j++) {
          number_c_u = number_c_u / (2 * H->h[j]);
          number_c_l = number_c_l * (2 * H->h[j]);
        }
        number_c_u = number_c_u / (2 * H->h[index]);

        unsigned int pos = 0;
        pos = int64_vector_add_one_i(c, index + 1, H);
        line_sieve_1(array, c, pseudo_Tqr, r, H, index,
                number_c_l, &ci, pos, matrix, f);
        for (uint64_t j = 1; j < number_c_u; j++) {
          pos = int64_vector_add_one_i(c, index + 1, H);
          line_sieve_1(array, c, pseudo_Tqr, r, H, index,
                  number_c_l, &ci, pos, matrix, f);
        }
      }
    }
#endif // SIEVE_TQR

#ifdef SIEVE_U
  for (uint64_t i = 0; i < fb->number_element_u; i++) {
    mpz_t ** Tqr = (mpz_t **)
      malloc(sizeof(mpz_t *) * (fb->factor_base_u[i]->ideal->h->deg));
    for (int j = 0; j < fb->factor_base_u[i]->ideal->h->deg; j++) {
      Tqr[j] = (mpz_t *) malloc(sizeof(mpz_t) * (H->t));
    }

    compute_Tqr_u(Tqr, matrix, H->t, fb->factor_base_u[i]);

    mpz_vector_t c;
    mpz_vector_init(c, H->t);
    for (unsigned int j = 0; j < H->t - 1; j++) {
      mpz_vector_setcoordinate_si(c, j, -H->h[j]);
    }
    mpz_vector_setcoordinate_si(c, H->t - 1, 0);

    sieve_u(array, Tqr, fb->factor_base_u[i], c, H);

    for (uint64_t j = 1; j < array->number_element; j++) {
      mpz_vector_add_one(c, H);
      sieve_u(array, Tqr, fb->factor_base_u[i], c, H);
    }
    mpz_vector_clear(c);
    for (unsigned int col = 0; col < H->t; col++) {
      for (int row = 0; row < fb->factor_base_u[i]->ideal->h->deg;
           row++) {
        mpz_clear(Tqr[row][col]);
      }
    }
    for (int j = 0; j < fb->factor_base_u[i]->ideal->h->deg; j++) {
      free(Tqr[j]);
    }
    free(Tqr);
  }
#endif // SIEVE_U
#endif // Draft for special-q_sieve

void find_index(uint64_array_ptr indexes, array_srcptr array,
    unsigned char thresh)
{
#ifdef MEAN_NORM_BOUND_SIEVE
  norm_bound_sieve = 0;
#endif // MEAN_NORM_BOUND_SIEVE

  uint64_t ind = 0;
  for (uint64_t i = 0; i < array->number_element; i++) {

#ifdef MEAN_NORM_BOUND_SIEVE
    norm_bound_sieve = norm_bound_sieve + (double)array->array[i];
#endif // MEAN_NORM_BOUND_SIEVE

    if (array->array[i] <= thresh) {
      ASSERT(ind < indexes->length);
      indexes->array[ind] = i;
      ind++;
    }
  }
  uint64_array_realloc(indexes, ind);
}

/*
 * Print a relation.
 *
 * factor: factorisation of a if the norm is smooth is the corresponding number
 *  field.
 * I:
 * L:
 * a: the polynomial in the original lattice.
 * t: dimension of the lattice.
 * V: number of number fields.
 * size:
 * assert_facto: 1 if the factorisation is correct, 0 otherwise.
 */
#ifdef ASSERT_FACTO
void printf_relation(factor_t * factor, unsigned int * I, unsigned int * L,
    mpz_poly_srcptr a, unsigned int t, unsigned int V, unsigned int size,
    unsigned int * assert_facto)
#else
void printf_relation(factor_t * factor, unsigned int * I, unsigned int * L,
    mpz_poly_srcptr a, unsigned int t, unsigned int V, unsigned int size)
#endif
{
  //Print a.
  printf("# ");
  mpz_poly_fprintf(stdout, a);

  //Print if the factorisation is good.
#ifdef ASSERT_FACTO
  unsigned int index = 0;
  printf("# ");
  for (unsigned int i = 0; i < V - 1; i++) {
    if (index < size) {
      if (i == L[index]) {
        if (I[index]) {
          printf("%u:", assert_facto[index]);
        } else {
          printf(":");
        }
        index++;
      } else {
        printf(":");
      }
    } else {
      printf(":");
    }
  }

  if (index < size) {
    if (V - 1 == L[index]) {
      if (I[index]) {
        printf("%u", assert_facto[index]);
      }
      index++;
    }
  }

  printf("\n");
#else
  unsigned int index = 0;
#endif

  //Print the coefficient of a.
  for (int i = 0; i < a->deg; i++) {
    gmp_printf("%Zd,", a->coeff[i]);
  }
  if ((int)t - 1 == a->deg) {
    gmp_printf("%Zd:", a->coeff[a->deg]);
  } else {
    gmp_printf("%Zd,", a->coeff[a->deg]);
    for (int i = a->deg + 1; i < (int)t - 1; i++) {
      printf("0,");
    }
    printf("0:");
  }

  //Print the factorisation in the different number field.
  index = 0;
  for (unsigned int i = 0; i < V - 1; i++) {
    if (index < size) {
      if (i == L[index]) {
        if (I[index]) {
          for (unsigned int j = 0; j < factor[index]->number - 1; j++) {
            gmp_printf("%Zd,", factor[index]->factorization[j]);
          }
          gmp_printf("%Zd:", factor[index]->factorization[
                       factor[index]->number - 1]);
        } else {
          printf(":");
        }
        index++;
      } else {
        printf(":");
      }
    } else {
      printf(":");
    }
  }

  if (index < size) {
    if (V - 1 == L[index]) {
      if (I[index]) {
        for (unsigned int j = 0; j < factor[index]->number - 1; j++) {
          gmp_printf("%Zd,", factor[index]->factorization[j]);
        }
        gmp_printf("%Zd", factor[index]->factorization[
                     factor[index]->number - 1]);
      }
      index++;
    }
  }

  printf("\n");
}

/*
 * A potential polynomial is found. Factorise it to verify if it gives a
 *  relation.
 *
 * a: the polynomial in the original lattice.
 * f: the polynomials that define the number fields.
 * lpb: the large prime bands.
 * L:
 * size:
 * t: dimension of the lattice.
 * V: number of number fields.
 */
void good_polynomial(mpz_poly_srcptr a, mpz_poly_t * f, mpz_t * lpb,
    unsigned int * L, unsigned int size, unsigned int t, unsigned int V)
{
  mpz_t * res = (mpz_t * ) malloc(sizeof(mpz_t) * size);
  factor_t * factor = (factor_t * ) malloc(sizeof(factor_t) * size);
  unsigned int * I = (unsigned int * ) malloc(sizeof(unsigned int) * size);

#ifdef ASSERT_FACTO
  unsigned int * assert_facto = (unsigned int * ) malloc(sizeof(unsigned int) * size);
#endif

  unsigned int find = 0;
  /* Not optimal */
  for (unsigned int i = 0; i < size; i++) {
    mpz_init(res[i]);
    norm_poly(res[i], f[L[i]], a);

#ifdef ASSERT_FACTO
    assert_facto[i] = gmp_factorize(factor[i], res[i]);
#else
    gmp_factorize(factor[i], res[i]);
#endif

    if (factor_is_smooth(factor[i], lpb[L[i]], 0)) {
      find++;
      I[i] = 1;
    } else {
      I[i] = 0;
    }
  }
  if (find >= 2) {
#ifdef ASSERT_FACTO
    printf_relation(factor, I, L, a, t, V, size, assert_facto);
#else
    printf_relation(factor, I, L, a, t, V, size);
#endif
  }

  for (unsigned int i = 0; i < size; i++) {
    factor_clear(factor[i]);
    mpz_clear(res[i]);
  }
  free(res);
  free(factor);
  free(I);
#ifdef ASSERT_FACTO
  free(assert_facto);
#endif
}

/*
 * Return the sum of the element in index. index has size V.
 *
 * index: array of index.
 * V: number of element, number of number fields.
 */
uint64_t sum_index(uint64_t * index, unsigned int V)
{
  uint64_t sum = 0;
  for (unsigned int i = 0; i < V; i++) {
    sum += index[i];
  }
  return sum;
}

/*
 * Find in all the indexes which one has the smallest value.
 *
 *
 *
 */
unsigned int find_indexes_min(unsigned int ** L,
    uint64_array_t * indexes, uint64_t * index, unsigned int V)
{
  * L = (unsigned int * ) malloc(sizeof(unsigned int) * (V));
  unsigned int i = 0;
  unsigned int size = 0;
  while(indexes[i]->length == 0)
  {
    i++;
  }
  uint64_t min = indexes[i]->array[index[i]];

  for (; i < V; i++) {
    if (index[i] < indexes[i]->length) {
      min = MIN(min, indexes[i]->array[index[i]]);
    }
  }
  for (i = 0; i < V; i++) {
    if (indexes[i]->length != 0) {
      if (min == indexes[i]->array[index[i]]) {
        (*L)[size] = i;
        size++;
      }
    }
  }
  * L = realloc(* L, size * sizeof(unsigned int));
  return size;
}

/*
 * For 
 */
void find_relation(uint64_array_t * indexes, uint64_t * index,
    uint64_t number_element, mpz_t * lpb, mat_Z_srcptr matrix, mpz_poly_t * f,
    sieving_bound_srcptr H, unsigned int V)
{
  unsigned int * L;
  unsigned int size = find_indexes_min(&L, indexes, index, V);
  if (size >= 2) {
    mpz_vector_t c;
    mpz_t gcd;
    mpz_init(gcd);
    mpz_vector_init(c, H->t);
    array_index_mpz_vector(c, indexes[L[0]]->array[index[L[0]]], H,
        number_element);

    mpz_poly_t a;
    mpz_poly_init(a, 0);
    mat_Z_mul_mpz_vector_to_mpz_poly(a, matrix, c);
    mpz_poly_content(gcd, a);

#ifdef NUMBER_SURVIVALS
      number_survivals++;
#endif // NUMBER_SURVIVALS

    //a must be irreducible.
    if (mpz_cmp_ui(gcd, 1) == 0 && a->deg > 0 &&
        mpz_cmp_ui(mpz_poly_lc_const(a), 0) > 0) {

#ifdef NUMBER_SURVIVALS
      number_survivals_facto++;
#endif // NUMBER_SURVIVALS

      good_polynomial(a, f, lpb, L, size, H->t, V);
    }

    mpz_poly_clear(a);

    mpz_clear(gcd);
    mpz_vector_clear(c);

    for (unsigned int i = 0; i < size; i++) {
      index[L[i]] = index[L[i]] + 1;
    }
  } else {
    index[L[0]] = index[L[0]] + 1;
  }
  free(L);
}

/*
 * To find the relations.
 *
 * indexes: for each number field, position where the norm is less as the
 *  threshold.
 * number_element: number of element in the sieving region.
 * lpb: large prime bands.
 * matrix: MqLLL.
 * f: polynomials that define the number fields.
 * H: sieving bound.
 * V: number of number fields.
 */
void find_relations(uint64_array_t * indexes, uint64_t number_element,
    mpz_t * lpb, mat_Z_srcptr matrix, mpz_poly_t * f, sieving_bound_srcptr H,
    unsigned int V)
{
  //Index is the current index to move in the indexes array.
  uint64_t * index = (uint64_t * ) malloc(sizeof(uint64_t) * V);
  /* Compute sum of the length of all the uint64_array. */
  uint64_t length_tot = 0;
  for (unsigned int i = 0; i < V; i++) {
    index[i] = 0;
    if (indexes[i]->length != 0) {
      length_tot += (indexes[i]->length - 1);
    }
  }

  if (0 != length_tot) {
    while(sum_index(index, V) < length_tot) {
      find_relation(indexes, index, number_element, lpb, matrix, f, H, V);
    }
    find_relation(indexes, index, number_element, lpb, matrix, f, H, V);
    free(index);
  } else {
    printf("# No relations\n");
  }
}

void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "H", "the sieving region");
  param_list_decl_usage(pl, "V", "number of number field");
  param_list_decl_usage(pl, "fbb0", "factor base bound on the number field 0");
  param_list_decl_usage(pl, "fbb1", "factor base bound on the number field 1");
  param_list_decl_usage(pl, "thresh0", "threshold on the number field 0");
  param_list_decl_usage(pl, "thresh1", "threshold on the number field 1");
  param_list_decl_usage(pl, "lpb0", "threshold on the number field 0");
  param_list_decl_usage(pl, "lpb1", "threshold on the number field 1");
  param_list_decl_usage(pl, "f0", "polynomial that defines the number field 0");
  param_list_decl_usage(pl, "f1", "polynomial that defines the number field 1");
  param_list_decl_usage(pl, "q_min", "minimum of the special-q");
  param_list_decl_usage(pl, "q_max", "maximum of the special-q");
  param_list_decl_usage(pl, "q_side", "side of the special-q");

  /* MNFS */

  param_list_decl_usage(pl, "fbb2", "factor base bound on the number field 2");
  param_list_decl_usage(pl, "fbb3", "factor base bound on the number field 3");
  param_list_decl_usage(pl, "fbb4", "factor base bound on the number field 4");
  param_list_decl_usage(pl, "fbb5", "factor base bound on the number field 5");
  param_list_decl_usage(pl, "fbb6", "factor base bound on the number field 6");
  param_list_decl_usage(pl, "fbb7", "factor base bound on the number field 7");
  param_list_decl_usage(pl, "fbb8", "factor base bound on the number field 8");
  param_list_decl_usage(pl, "fbb9", "factor base bound on the number field 9");
  param_list_decl_usage(pl, "thresh2", "threshold on the number field 2");
  param_list_decl_usage(pl, "thresh3", "threshold on the number field 3");
  param_list_decl_usage(pl, "thresh4", "threshold on the number field 4");
  param_list_decl_usage(pl, "thresh5", "threshold on the number field 5");
  param_list_decl_usage(pl, "thresh6", "threshold on the number field 6");
  param_list_decl_usage(pl, "thresh7", "threshold on the number field 7");
  param_list_decl_usage(pl, "thresh8", "threshold on the number field 8");
  param_list_decl_usage(pl, "thresh9", "threshold on the number field 9");
  param_list_decl_usage(pl, "lpb2", "threshold on the number field 2");
  param_list_decl_usage(pl, "lpb3", "threshold on the number field 3");
  param_list_decl_usage(pl, "lpb4", "threshold on the number field 4");
  param_list_decl_usage(pl, "lpb5", "threshold on the number field 5");
  param_list_decl_usage(pl, "lpb6", "threshold on the number field 6");
  param_list_decl_usage(pl, "lpb7", "threshold on the number field 7");
  param_list_decl_usage(pl, "lpb8", "threshold on the number field 8");
  param_list_decl_usage(pl, "lpb9", "threshold on the number field 9");
  param_list_decl_usage(pl, "f2", "polynomial that defines the number field 2");
  param_list_decl_usage(pl, "f3", "polynomial that defines the number field 3");
  param_list_decl_usage(pl, "f4", "polynomial that defines the number field 4");
  param_list_decl_usage(pl, "f5", "polynomial that defines the number field 5");
  param_list_decl_usage(pl, "f6", "polynomial that defines the number field 6");
  param_list_decl_usage(pl, "f7", "polynomial that defines the number field 7");
  param_list_decl_usage(pl, "f8", "polynomial that defines the number field 8");
  param_list_decl_usage(pl, "f9", "polynomial that defines the number field 9");
}

/*
 * Initialise the parameters of the special-q sieve.
 *
 * f: the V functions to define the number fields.
 * fbb: the V factor base bounds.
 * t: dimension of the lattice.
 * H: sieving bound.
 * q_min: lower bound of the special-q range.
 * q_max: upper bound of the special-q range.
 * thresh: the V threshold.
 * lpb: the V large prime bounds.
 * array: array in which the norms are stored.
 * matrix: the Mq matrix (set with zero coefficients).
 * q_side: side of the special-q.
 * V: number of number fields.
 */
void initialise_parameters(int argc, char * argv[], mpz_poly_t ** f,
    uint64_t ** fbb, factor_base_t ** fb, sieving_bound_ptr H,
    uint64_t * q_min, uint64_t * q_max, unsigned char ** thresh, mpz_t ** lpb,
    array_ptr array, mat_Z_ptr matrix, unsigned int * q_side, unsigned int * V)
{
  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  FILE * fpl;
  char * argv0 = argv[0];

  argv++, argc--;
  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

    /* Could also be a file */
    if ((fpl = fopen(argv[0], "r")) != NULL) {
      param_list_read_stream(pl, fpl, 0);
      fclose(fpl);
      argv++,argc--;
      continue;
    }

    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    param_list_print_usage(pl, argv0, stderr);
    exit (EXIT_FAILURE);
  }

  param_list_parse_uint(pl, "V", V);
  ASSERT(* V >= 2 && * V < 11);

  unsigned int t;
  int * r;
  param_list_parse_int_list_size(pl, "H", &r, &t, ".,");
  ASSERT(t > 2);
  ASSERT(t == 3); // TODO: remove as soon as possible.
  sieving_bound_init(H, t);
  for (unsigned int i = 0; i < t; i++) {
    sieving_bound_set_hi(H, i, (unsigned int) r[i]);
  }
  free(r);

  //TODO: something strange here.
  * fbb = malloc(sizeof(uint64_t) * (* V));
  * thresh = malloc(sizeof(unsigned char) * (* V));
  * fb = malloc(sizeof(factor_base_t) * (* V));
  * f = malloc(sizeof(mpz_poly_t) * (* V));
  * lpb = malloc(sizeof(mpz_t) * (* V));

  for (unsigned int i = 0; i < * V; i++) {
    char str [4];
    sprintf(str, "fbb%u", i);
    param_list_parse_uint64(pl, str, &fbb[0][i]);
    factor_base_init((*fb)[i], fbb[0][i], fbb[0][i], fbb[0][i]);
  }

  for (unsigned int i = 0; i < * V; i++) {
    char str [7];
    sprintf(str, "thresh%u", i);
    param_list_parse_uchar(pl, str, (* thresh) + i);
  }

  for (unsigned int i = 0; i < * V; i++) {
    char str [2];
    sprintf(str, "f%u", i);
    param_list_parse_mpz_poly(pl, str, (*f)[i], ".,");
  }

  for (unsigned int i = 0; i < * V; i++) {
    char str [4];
    sprintf(str, "lpb%u", i);
    mpz_init((*lpb)[i]);
    param_list_parse_mpz(pl, str,(*lpb)[i]);
    ASSERT(mpz_cmp_ui((*lpb)[i], fbb[0][i]) >= 0);
  }

  uint64_t number_element = sieving_bound_number_element(H);
  ASSERT(number_element >= 6);

  //TODO: q_side is an unsigned int.
  param_list_parse_int(pl, "q_side", (int *)q_side);
  ASSERT(* q_side < * V);
  param_list_parse_uint64(pl, "q_min", q_min);
  ASSERT(* q_min > fbb[0][* q_side]);
  param_list_parse_uint64(pl, "q_max", q_max);
  ASSERT(* q_min < * q_max);

  param_list_clear(pl);

  array_init(array, number_element);

  mat_Z_init(matrix, t, t);
}

/*
  The main.
*/
int main(int argc, char * argv[])
{
  unsigned int V;
  mpz_poly_t * f;
  uint64_t * fbb;
  sieving_bound_t H;
  uint64_t q_min;
  uint64_t q_max;
  unsigned int q_side;
  unsigned char * thresh;
  mpz_t * lpb;
  array_t array;
  mat_Z_t matrix;
  factor_base_t * fb;
  uint64_t q;

  initialise_parameters(argc, argv, &f, &fbb, &fb, H, &q_min, &q_max,
                        &thresh, &lpb, array, matrix, &q_side, &V);

#ifdef PRINT_PARAMETERS
  printf("# H =\n");
  sieving_bound_fprintf_comment(stdout, H);
  printf("# V = %u\n", V);
  for (unsigned int i = 0; i < V; i++) {
    printf("# fbb%u = %" PRIu64 "\n", i, fbb[i]);
  }
  for (unsigned int i = 0; i < V; i++) {
    printf("# thresh%u = %u\n", i, (unsigned int)thresh[i]);
  }
  for (unsigned int i = 0; i < V; i++) {
    gmp_printf("# lpb%u = %Zd\n", i, lpb[i]);
  }
  for (unsigned int i = 0; i < V; i++) {
    printf("# f%u = ", i);
    mpz_poly_fprintf(stdout, f[i]);
  }
  printf("# q_min = %" PRIu64 "\n", q_min);
  printf("# q_max = %" PRIu64 "\n", q_min);
  printf("# q_side = %u\n", q_side);
#endif // PRINT_PARAMETERS

  uint64_array_t * indexes =
    (uint64_array_t * ) malloc(sizeof(uint64_array_t) * V);

  double ** pre_compute = (double ** ) malloc(sizeof(double * ) * V);
  for (unsigned int i = 0; i < V; i++) {
    pre_compute[i] = (double * ) malloc((H->t) * sizeof(double));
    pre_computation(pre_compute[i], f[i], H->t);
  }

  double ** time = (double ** ) malloc(sizeof(double * ) * V);
  for (unsigned int i = 0; i < V; i++) {
    time[i] = (double * ) malloc(sizeof(double) * 3);
  }
  double sec_tot;
  double sec_cofact;

#ifdef NUMBER_SURVIVALS
  uint64_t * numbers_survivals = (uint64_t * ) malloc(sizeof(uint64_t) * V);
#endif // NUMBER_SURVIVALS

#ifdef TRACE_POS
  file = fopen("TRACE_POS.txt", "w+");
  fprintf(file, "TRACE_POS: %d\n", TRACE_POS);
#endif // TRACE_POS

  double sec = seconds();
  makefb(fb, f, fbb, H->t, lpb, V);

  printf("# Time for makefb: %f.\n", seconds() - sec);

  ASSERT(q_min >= fbb[q_side]);
  ideal_1_t special_q;
  gmp_randstate_t state;
  mpz_t a;
  mpz_poly_factor_list l;

  mpz_poly_factor_list_init(l);
  gmp_randinit_default(state);
  mpz_init(a);
  ideal_1_init(special_q);

  //Pass all the prime less than q_min.
  for (q = 2; q < q_min; q = getprime(q)) {}

  for ( ; q <= q_max; q = getprime(q)) {
    mpz_set_si(a, q);
    mpz_poly_factor(l, f[q_side], a, state);
    for (int i = 0; i < l->size ; i++) {
      //Only deal with special-q of degre
      if (l->factors[i]->f->deg == 1) {
        ideal_1_set_part(special_q, q, l->factors[i]->f, H->t);
        printf("# Special-q: q: %" PRIu64 ", g: ", q);
        mpz_poly_fprintf(stdout, l->factors[i]->f);

#ifdef TRACE_POS
        fprintf(file, "Special-q: q: %" PRIu64 ", g: ", q);
        mpz_poly_fprintf(file, l->factors[i]->f);
#endif // TRACE_POS

        sec = seconds();
        sec_tot = sec;

        /* LLL part */
        build_Mq_ideal_1(matrix, special_q);

#ifdef TRACE_POS
        fprintf(file, "Mq:\n");
        mat_Z_fprintf(file, matrix);
#endif // TRACE_POS

        mat_Z_LLL_transpose(matrix);

        /* /\* */
        /*   TODO: continue here. */
        /*   If the last coefficient of the last vector is negative, we do not */
        /*    sieve in the good direction. We therefore take the opposite vector */
        /*    because we want that a = MqLLL * c, with c in the sieving region */
        /*    has the last coordinate positive. */
        /* *\/ */
        /* if (mpz_cmp_ui(matrix->coeff[matrix->NumRows] */
        /*                [matrix->NumCols], 0) < 0) { */
        /*   for (unsigned int rows = 1; rows <= matrix->NumRows; rows++) { */
        /*     mpz_mul_si(matrix->coeff[rows][matrix->NumCols], */
        /*                matrix->coeff[rows][matrix->NumCols], -1); */
        /*   } */
        /* } */

#ifdef TRACE_POS
        fprintf(file, "MqLLL:\n");
        mat_Z_fprintf(file, matrix);
#endif // TRACE_POS

#ifdef MEAN_NORM_BOUND
        double * norms_bound = (double * ) malloc(sizeof(double) * V);
#endif // MEAN_NORM_BOUND

#ifdef MEAN_NORM
        double * norms = (double * ) malloc(sizeof(double) * V);
#endif // MEAN_NORM

#ifdef MEAN_NORM_BOUND_SIEVE
        double * norms_bound_sieve = (double * ) malloc(sizeof(double) * V);
#endif // MEAN_NORM_BOUND_SIEVE

        for (unsigned int j = 0; j < V; j++) {
          sec = seconds();
          uint64_array_init(indexes[j], array->number_element);
          init_norm_1(array, pre_compute[j], H, matrix, f[j], special_q,
                      !(j ^ q_side));
          time[j][0] = seconds() - sec;

#ifdef MEAN_NORM_BOUND
          norms_bound[j] = norm_bound;
#endif // MEAN_NORM_BOUND

#ifdef MEAN_NORM
          norms[j] = norm;
#endif // MEAN_NORM

          sec = seconds();
          special_q_sieve(array, matrix, fb[j], H, f[j]);
          time[j][1] = seconds() - sec;
          sec = seconds();
          find_index(indexes[j], array, thresh[j]);
          time[j][2] = seconds() - sec;

#ifdef MEAN_NORM_BOUND_SIEVE
        norms_bound_sieve[j] = norm_bound_sieve;
#endif // MEAN_NORM_BOUND_SIEVE

#ifdef TRACE_POS
        fprintf(file, "********************\n");
#endif // TRACE_POS
        }

        sec = seconds();
        find_relations(indexes, array->number_element, lpb, matrix, f, H, V);
        sec_cofact = seconds() - sec;

        for (unsigned j = 0; j < V; j++) {

#ifdef NUMBER_SURVIVALS
          numbers_survivals[j] = indexes[j]->length;
#endif // NUMBER_SURVIVALS

          uint64_array_clear(indexes[j]);
        }

        printf("# Time for this special-q: %fs.\n", seconds() - sec_tot);
        for (unsigned int j = 0; j < V; j++) {
          printf("# Time to init norm %d: %fs.\n", j, time[j][0]);

#ifdef MEAN_NORM
          printf("# Mean of the norm (bit size) %d: %f.\n",
                 j, log2(norms[j] / (double) array->number_element));
#endif // MEAN_NORM

#ifdef MEAN_NORM_BOUND
          printf("# Mean of the bound of the norm (bit size) %d: %f.\n",
                 j, log2(norms_bound[j] / (double) array->number_element));
#endif // MEAN_NORM_BOUND

#ifdef MEAN_NORM_BOUND_SIEVE
          printf("# Mean of the bit size of the bound of the residual norm %d: %f.\n",
                 j, norms_bound_sieve[j] / (double) array->number_element);
#endif // MEAN_NORM_BOUND_SIEVE

          printf("# Time to sieve %d: %fs.\n", j, time[j][1]);
          printf("# Time to find indexes %d: %fs.\n", j, time [j][2]);

#ifdef NUMBER_SURVIVALS
          printf("# Number of survivals %d: %" PRIu64 ".\n", j,
                 numbers_survivals[j]);
#endif // NUMBER_SURVIVALS
        }

#ifdef NUMBER_SURVIVALS
        free(numbers_survivals);
#endif // NUMBER_SURVIVALS

        printf("# Time to factorize: %fs.\n", sec_cofact);

#ifdef NUMBER_SURVIVALS
        printf("# Number total of survivals: %" PRIu64 ".\n",
               number_survivals);
        printf("# Number total of polynomial a survivals: %" PRIu64 ".\n",
               number_survivals_facto);
#endif // NUMBER_SURVIVALS

        printf("# ----------------------------------------\n");

#ifdef MEAN_NORM
        free(norms);
#endif // MEAN_NORM

#ifdef MEAN_NORM_BOUND
        free(norms_bound);
#endif // MEAN_NORM_BOUND

#ifdef MEAN_NORM_BOUND_SIEVE
        free(norms_bound_sieve);
#endif // MEAN_NORM_BOUND_SIEVE

#ifdef TRACE_POS
        fprintf(file, "----------------------------------------\n");
#endif // TRACE_POS

      }
    }
  }

  mpz_poly_factor_list_clear(l);
  ideal_1_clear(special_q, H->t);
  gmp_randclear(state);
  mpz_clear(a);
  getprime(0);

#ifdef TRACE_POS
  fclose(file);
#endif // TRACE_POS

  mat_Z_clear(matrix);
  for (unsigned int i = 0; i < V; i++) {
    mpz_clear(lpb[i]);
    mpz_poly_clear(f[i]);
    factor_base_clear(fb[i], H->t);
    free(pre_compute[i]);
    free(time[i]);
  }
  free(time);
  free(indexes);
  free(pre_compute);
  array_clear(array);
  free(lpb);
  free(f);
  free(fb);
  sieving_bound_clear(H);
  free(fbb);
  free(thresh);

  return 0;
}
