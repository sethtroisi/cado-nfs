#include <cado.h>
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

/*
  p = 10822639589
  n = 6
  f0 = 4598*x^6 + 4598*x^5 + 4597*x^4 + 4598*x^3 + 4595*x^2 + 4598*x + 4596
  4596,4598,4585,4598,4597,4598,4598
  f1 = -2353645*x^6 - 2353645*x^5 + 638*x^4 - 2353645*x^3 + 4709204*x^2 -
    2353645*x + 2354921
  2354921,-2353645,4709204,-2353645,638,-2353645,-2353645
  t = 4
  lpb = 28
  N_a = 2^40
*/

/*
  Different mode:
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
#endif

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

/*
  Build the matrix Mq for an ideal (q, g) with deg(g) = 1. The matrix need Tq
   with the good coefficients, i.e. t the coeffecients not reduced modulo q.

  matrix: the matrix with the good initialization for the number of rows and
   columns. We can not assert here if the numbers of values for Tq and matrix is
   good or not.
  ideal: the ideal (q, g) with deg(g) = 1.
*/
void build_Mq_ideal_1(mat_Z_ptr matrix, ideal_1_srcptr ideal)
{
  for (unsigned int row = 1; row <= matrix->NumRows; row++) {
    for (unsigned int col = 1; col <= matrix->NumCols; col++) {
      mat_Z_set_coeff_uint64(matrix, 0, row, col);
    }
  }
  mat_Z_set_coeff_uint64(matrix, ideal->ideal->r, 1, 1);
  for (unsigned int row = 2; row <= matrix->NumRows; row++) {
    mat_Z_set_coeff_uint64(matrix, 1, row, row);
  }
  for (unsigned int col = 2; col <= matrix->NumCols; col++) {
    mat_Z_set_coeff(matrix, ideal->Tr[col - 2], 1, col);
  }
}

/*
  Build the matrix Mq for an ideal (q, g) with deg(g) > 1. This is probably
   non-sens.
*/
void build_Mq_ideal_u(mat_Z_ptr matrix, ideal_u_srcptr ideal)
{
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
  Precompute some part of the computation of the bound on the resultant between
   f and ha. We precompute (deg(f) + 1) ^ (i/2) * (i + 1) ^ (deg(f) / 2) *
   infinite_norm(f), with i the degree of ha.
  pre_compute: an array in which we store the precompute elements.
  f: a polynomial that defines the side.
  t: t gives the degree of ha.
*/
void pre_computation(double * pre_compute, mpz_poly_srcptr f, unsigned int t)
{
  pre_compute[0] = 1;
  mpz_t infinite_norm_f_tmp;
  mpz_init(infinite_norm_f_tmp);
  mpz_poly_infinite_norm(infinite_norm_f_tmp, f);
  double infinite_norm_f = mpz_get_d(infinite_norm_f_tmp);
  mpz_clear(infinite_norm_f_tmp);
  for (unsigned int i = 1; i < t; i++) {
    pre_compute[i] = pow(infinite_norm_f, (double)i) *
      pow((double)(i + 1), ((double)f->deg / 2)) * pow((double)(f->deg + 1),
                                                       ((double)i / 2));
  }
}

/*
  Compute the norm of a case in the array.

  array: array in which we store the norm.
  i: index of the array in which we store the computed norm.
  a: polynomial equal to matrix * vector.
  pre_compute: array with precomputed value to compute upper bound of the norm.
  H: sieving interval.
  matrix: to compute the current value of a.
  beg: index of the upper changed coordinate during the addition of 1 in the
   sieving region for vector.
  f: function that defines the number field.
  ideal: the special-q we set.
  special_q: if the special-q is set in this side, special-q is equal to 1.
  vector: only for TRACE_POS mode.
*/
void init_each_case(array_ptr array, uint64_t i, int64_poly_ptr a,
                    double * pre_compute, sieving_interval_srcptr H,
                    mat_int64_srcptr matrix, unsigned int beg,
                    mpz_poly_srcptr f, ideal_1_srcptr ideal, int special_q,
                    MAYBE_UNUSED int64_vector_srcptr vector)
{
  double bound_resultant;

  //a = a + sum(matrix[j][beg + 1], j, 1, H->t + 1) * x^beg.
  for (unsigned int j = 0; j < H->t; j++) {
    int64_poly_setcoeff(a, j, a->coeff[j] + matrix->coeff[j + 1][beg + 1]);
  }

  //Substract -2 * H[i] when you change an other coordinate than c0.
  for (unsigned int k = 0; k < beg; k++) {
    for (unsigned int j = 0; j < H->t; j++) {
      int64_t tmp = a->coeff[j];
      tmp = tmp - matrix->coeff[j + 1][k + 1] * 2 * H->h[j];
      int64_poly_setcoeff(a, j, tmp);
    }
  }

  ASSERT(a->deg >= -1);

#ifndef NDEBUG
  int64_poly_t a_tmp;
  int64_poly_init(a_tmp, -1);
  mat_int64_mul_int64_vector_to_int64_poly(a_tmp, matrix, vector);
  ASSERT(int64_poly_equal(a, a_tmp) == 0);
  int64_poly_clear(a_tmp);
#endif // NDEBUG

  if (a->deg > 0) {
    uint64_t tmp = 0;
    int64_poly_infinite_norm(&tmp, a);
    bound_resultant = pow((double)tmp, f->deg) * pre_compute[a->deg];
  } else {
    bound_resultant = 1;
  }

#ifdef MEAN_NORM_BOUND
  norm_bound = norm_bound + bound_resultant;
#endif // MEAN_NORM_BOUND

#ifdef MEAN_NORM
  mpz_poly_t b;
  mpz_poly_init(b, -1);
  mpz_t res;
  mpz_init(res);
  int64_poly_to_mpz_poly(b, a);
  mpz_poly_resultant(res, f, b);
  mpz_abs(res, res);
  norm = norm + mpz_get_d(res);
  mpz_poly_clear(b);
  mpz_clear(res);
#endif // MEAN_NORM

  //WARNING: an assert here is necessary.
  if (special_q) {
    array->array[i] = (unsigned char) (log2(bound_resultant) -
                                       ideal->log);

#ifdef TRACE_POS
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
      mpz_poly_resultant(res, f, tmp);
      factor_t factor;
      mpz_abs(res, res);
      gmp_factorize(factor, res);
      gmp_fprintf(file, "Resultant: %Zd\n", res);
      fprintf(file, "Factorization: ");
      factor_fprintf(file, factor);
      factor_clear(factor);
      mpz_clear(res);
      fprintf(file, "Initialization (without special-q): %u\n",
              (unsigned char) log2(bound_resultant));
      fprintf(file, "Initialization (with special-q): %u\n", array->array[i]);
      mpz_poly_clear(tmp);
    }
#endif

  } else {
    array->array[i] = (unsigned char) log(bound_resultant);

#ifdef TRACE_POS
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
      mpz_poly_resultant(res, f, tmp);
      factor_t factor;
      mpz_abs(res, res);
      gmp_factorize(factor, res);
      gmp_fprintf(file, "Norm: %Zd\n", res);
      fprintf(file, "Factorization: ");
      factor_fprintf(file, factor);
      factor_clear(factor);
      mpz_clear(res);
      mpz_poly_clear(tmp);
      fprintf(file, "Initialization: %u\n", array->array[i]);
    }
#endif

  }
}

/*
  Init the norm for a special-q (q, g) with deg(g) = 1.

  array: the array in which the norms are initialized.
  pre_compute: the array of precomputed value obtained with pre_computation.
  H: the sieving interval which gives the sieving region.
  matrix: the MqLLL matrix.
  f: the f which defines the side.
  ideal: the special-q.
  special_q: 0 if there is no special-q in this side, else 1.
*/
void init_norm_1(array_ptr array, double * pre_compute,
                 sieving_interval_srcptr H, mat_Z_srcptr matrix,
                 mpz_poly_srcptr f, ideal_1_srcptr ideal, int special_q)
{
  ASSERT(special_q == 0 || special_q == 1);

  int64_vector_t vector;
  int64_vector_init(vector, H->t);
  int64_vector_setcoordinate(vector, 0, -(int64_t)H->h[0] - 1);
  for (unsigned int i = 1; i < H->t - 1; i++) {
    int64_vector_setcoordinate(vector, i, -(int64_t)H->h[i]);
  }
  int64_vector_setcoordinate(vector, (int64_t)H->t - 1, 0);
  unsigned int beg = 0;

  mat_int64_t matrix_int;
  mat_int64_init(matrix_int, matrix->NumRows, matrix->NumCols);
  mat_Z_to_mat_int64(matrix_int, matrix);

  int64_poly_t a;
  int64_poly_init(a, H->t - 1);
  mat_int64_mul_int64_vector_to_int64_poly(a, matrix_int, vector);

#ifdef MEAN_NORM_BOUND
    norm_bound = 0;
#endif // MEAN_NORM_BOUND

#ifdef MEAN_NORM
  norm = 0;
#endif // MEAN_NORM

  int64_vector_setcoordinate(vector, 0, -(int64_t)H->h[0]);
  init_each_case(array, 0, a, pre_compute, H, matrix_int, beg, f, ideal, special_q,
                 vector);

  for (uint64_t i = 1; i < array->number_element; i++) {
    beg = int64_vector_add_one_i(vector, 0, H);
    init_each_case(array, i, a, pre_compute, H, matrix_int, beg, f, ideal,
                   special_q, vector);
  }

  int64_poly_clear(a);
  int64_vector_clear(vector);
  mat_int64_clear(matrix_int);
}

/*
  Compute Tqr for an ideal (r, h) with deg(h) = 1.
  In fact, Tqr is not just computed. For example, if Tqr[0] != 0, the output is
   equal to [-Tqr[0]^{-1} mod r = a, a * Tqr[1] mod r, …].
  The PRINT_TQR mode print the true Tqr matrix.

  Tqr: Tqr is just a line in this case.
  matrix: MqLLL.
  t: dimension of the lattice.
  ideal: special-q.
*/
void compute_Tqr_1(uint64_t * Tqr, mat_Z_srcptr matrix, unsigned int t,
                   ideal_1_srcptr ideal)
{
  mpz_t tmp;
  mpz_init(tmp);
  unsigned int i = 0;
  Tqr[i] = 0;
  //Tqr = Mq,1 - Tr * Mq,2.

#ifdef PRINT_TQR
  mpz_t invert;
  mpz_init(invert);
  printf("Tqr = [");
#endif // PRINT_TQR

  for (unsigned int j = 0; j < t; j++) {
    mpz_set(tmp, matrix->coeff[1][j + 1]);
    for (unsigned int k = 0; k < t - 1; k++) {

      mpz_submul(tmp, ideal->Tr[k],
                 matrix->coeff[k + 2][j + 1]);
    }

    if (Tqr[i] == 0) {
      mpz_mod_ui(tmp, tmp, ideal->ideal->r);
      if (mpz_cmp_ui(tmp, 0) != 0) {
        mpz_invert_ui(tmp, tmp, ideal->ideal->r);

#ifdef PRINT_TQR
        mpz_set(invert, tmp);
        gmp_printf("1, ", invert);
#endif // PRINT_TQR

        mpz_mul_si(tmp, tmp, -1);
        mpz_mod_ui(tmp, tmp, ideal->ideal->r);
        ASSERT(mpz_cmp_ui(tmp, ideal->ideal->r) <= 0);
        Tqr[j] = mpz_get_ui(tmp);
        i = j;
      } else {

#ifdef PRINT_TQR
        printf("0, ");
#endif // PRINT_TQR

        Tqr[j] = 0;
      }
    } else {

#ifdef PRINT_TQR
      mpz_t tmp2;
      mpz_init(tmp2);
      mpz_mul(tmp2, tmp, invert);
      mpz_mod_ui(tmp2, tmp2, ideal->ideal->r);
      gmp_printf("%Zd, ", tmp2);
      mpz_clear(tmp2);
#endif // PRINT_TQR

      mpz_mul_ui(tmp, tmp, Tqr[i]);
      mpz_mod_ui(tmp, tmp, ideal->ideal->r);
      Tqr[j] = mpz_get_ui(tmp);
    }
  }
  mpz_clear(tmp);

#ifdef PRINT_TQR
  printf("\b\b] mod %" PRIu64 ", h: ", ideal->ideal->r);
  mpz_poly_fprintf(stdout, ideal->ideal->h);
  mpz_clear(invert);
#endif // PRINT_TQR
}

/*
  WARNING: Need to assert if the function does what it is supposed to do.
  Compute Tqr for an ideal (r, h) with deg(h) > 1.

  Tqr: Tqr is a matrix in this case.
  matrix: MqLLL.
  t: dimension of the lattice.
  ideal: (r, h).
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

/*
  Sieve for a special-q of degree 1 and Tqr with a zero coefficient at the
  first place. This function sieve a c in the q lattice.

  array: in which we store the norms.
  c: element of the q lattice.
  ideal: an ideal with r < q.
  ci: the possible first coordinate of c to have c in the sieving region.
  H: the sieving interval.
  i: index of the first non-zero coefficient in Tqr.
  number_c_l: number of possible c with the same ci, ci+1, …, ct.
*/
void sieve_ci(array_ptr array, int64_vector_ptr c, ideal_1_srcptr ideal,
              int64_t ci, sieving_interval_srcptr H, unsigned int i,
              uint64_t number_c_l)
{
  uint64_t index = 0;

  if (ci <= (int64_t)H->h[i]) {
    int64_vector_setcoordinate(c, i, ci);
    array_int64_vector_index(&index, c, H, array->number_element);
    array->array[index] = array->array[index] - ideal->log;

#ifdef NUMBER_HIT
    if (number_c_l == 0) {
      number_of_hit = number_of_hit + 1;
    } else {
      number_of_hit = number_of_hit + number_c_l;
    }
#endif // NUMBER_HIT

#ifdef TRACE_POS
    if (index == TRACE_POS) {
      fprintf(file, "The ideal is: ");
      ideal_fprintf(file, ideal->ideal);
      fprintf(file, "log: %u\n", ideal->log);
      fprintf(file, "The new value of the norm is %u.\n", array->array[index]);
    }
#endif

    for (uint64_t k = 1; k < number_c_l; k++) {
      array->array[index + k] = array->array[index + k] - ideal->log;

#ifdef TRACE_POS
      if (index + k == TRACE_POS) {
        fprintf(file, "The ideal is: ");
        ideal_fprintf(file, ideal->ideal);
        fprintf(file, "log: %u\n", ideal->log);
        fprintf(file, "The new value of the norm is %u.\n",
                array->array[index + k]);
      }
#endif
    }

    int64_t tmp = ci;
    tmp = tmp + (int64_t)ideal->ideal->r;

    while(tmp <= (int64_t)H->h[i]) {
      index = index + ideal->ideal->r + number_c_l;

      array->array[index] = array->array[index] - ideal->log;

#ifdef NUMBER_HIT
    if (number_c_l == 0) {
      number_of_hit = number_of_hit + 1;
    } else {
      number_of_hit = number_of_hit + number_c_l;
    }
#endif // NUMBER_HIT

#ifdef TRACE_POS
      if (index == TRACE_POS) {
        fprintf(file, "The ideal is: ");
        ideal_fprintf(file, ideal->ideal);
        fprintf(file, "log: %u\n", ideal->log);
        fprintf(file, "The new value of the norm is %u.\n",
                array->array[index]);
      }
#endif

      for (uint64_t k = 1; k < number_c_l; k++) {
        array->array[index + k] = array->array[index + k] - ideal->log;

#ifdef TRACE_POS
        if (index + k == TRACE_POS) {
          fprintf(file, "The ideal is: ");
        ideal_fprintf(file, ideal->ideal);
        fprintf(file, "log: %u\n", ideal->log);
          fprintf(file, "The new value of the norm is %u.\n",
                  array->array[index + k]);
        }
#endif

      }
      tmp = tmp + (int64_t)ideal->ideal->r;
    }
  }
}

void sieve_1(array_ptr array, int64_vector_ptr c, uint64_t * Tqr,
                 ideal_1_srcptr ideal, sieving_interval_srcptr H,
                 unsigned int i, uint64_t number_c_l, int64_t * ci,
                 unsigned int pos)
{
  ASSERT(pos >= i);

  for (unsigned int j = i + 1; j < pos; j++) {
    * ci = * ci - ((int64_t)Tqr[j] * 2 * (int64_t)H->h[j]);
    if (* ci >= (int64_t)ideal->ideal->r || * ci < 0) {
      * ci = * ci % (int64_t)ideal->ideal->r;
    }
  }
  * ci = * ci + Tqr[pos];
  if (* ci >= (int64_t)ideal->ideal->r) {
    * ci = * ci - (int64_t)ideal->ideal->r;
  }
  if (* ci < 0) {
    * ci = * ci + (int64_t)ideal->ideal->r;
  }
  ASSERT(* ci >= 0 && * ci < (int64_t)ideal->ideal->r);

  int64_t lb = (i < H->t - 1) ? -(int64_t)H->h[i] : 0;
  int64_t k = (lb - * ci) / (-(int64_t) ideal->ideal->r);
  if ((lb - * ci) > 0) {
    k--;
  }
  int64_t ci_tmp = -k * (int64_t)ideal->ideal->r + * ci;
  ASSERT(ci_tmp >= lb);

  sieve_ci(array, c, ideal, ci_tmp, H, i, number_c_l);
}

#ifdef SIEVE_U
void sieve_u(array_ptr array, mpz_t ** Tqr, ideal_u_srcptr ideal,
             mpz_vector_srcptr c, sieving_interval_srcptr H)
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
    array_mpz_vector_index(&index, c, H, array->number_element);
    array->array[index] = array->array[index] -
      ideal->log;

#ifdef TRACE_POS
    if (index == TRACE_POS) {
      fprintf(file, "The ideal is: ");
      ideal_u_fprintf(file, ideal, H->t);
      fprintf(file, "The new value of the norm is %u.\n", array->array[index]);
    }
#endif
  }

  for (int row = 0; row < ideal->ideal->h->deg;
       row++) {
    mpz_clear(tmp[row]);
  }
  free(tmp);
}
#endif // SIEVE_U

void special_q_sieve(array_ptr array, mat_Z_srcptr matrix,
                     factor_base_srcptr fb, sieving_interval_srcptr H)

{
#ifdef TIMER_SIEVE
  double time_tqr = 0;
  double time_sieve = 0;
  double sec = 0;
#endif // TIMER_SIEVE

  for (uint64_t i = 0; i < fb->number_element_1; i++) {

#ifdef TIMER_SIEVE
    if (fb->factor_base_1[i]->ideal->r > TIMER_SIEVE) {
      sec = seconds();
    }
#endif // TIMER_SIEVE

    uint64_t * Tqr = (uint64_t *) malloc(sizeof(uint64_t) * (H->t));
    compute_Tqr_1(Tqr, matrix, H->t, fb->factor_base_1[i]);

#ifdef TIMER_SIEVE
    if (fb->factor_base_1[i]->ideal->r > TIMER_SIEVE) {
      time_tqr = time_tqr + seconds() - sec;
    }
#endif // TIMER_SIEVE

    int64_vector_t c;
    int64_vector_init(c, H->t);
    for (unsigned int j = 0; j < H->t - 1; j++) {
      int64_vector_setcoordinate(c, j, -(int64_t)H->h[j]);
    }
    int64_vector_setcoordinate(c, H->t - 1, 0);

    if (Tqr[0] != 0) {
      unsigned int pos = 0;

      int64_t c0 = 0;
      if (H->t == 2) {
        int64_vector_setcoordinate(c, 1, (int64_t)(- 1));
      } else {
        int64_vector_setcoordinate(c, 1, c->c[1] - 1);
      }
      for (unsigned int j = 1; j < H->t; j++) {
        c0 = c0 + (int64_t)Tqr[j] * c->c[j];
        if (c0 >= (int64_t)fb->factor_base_1[i]->ideal->r || c0 < 0) {
          c0 = c0 % (int64_t)fb->factor_base_1[i]->ideal->r;
        }
      }
      if (c0 < 0) {
        c0 = c0 + (int64_t)fb->factor_base_1[i]->ideal->r;
      }
      ASSERT(c0 >= 0);
      uint64_t number_c = array->number_element / (2 * H->h[0] + 1);

#ifdef TIMER_SIEVE
      if (fb->factor_base_1[i]->ideal->r > TIMER_SIEVE) {
        sec = seconds();
      }
#endif // TIMER_SIEVE

      pos = int64_vector_add_one_i(c, 1, H);
      sieve_1(array, c, Tqr, fb->factor_base_1[i], H, 0, 0, &c0, pos);

      for (uint64_t j = 1; j < number_c;  j++) {
        pos = int64_vector_add_one_i(c, 1, H);
        sieve_1(array, c, Tqr, fb->factor_base_1[i], H, 0, 0, &c0, pos);
      }

#ifdef TIMER_SIEVE
      if (fb->factor_base_1[i]->ideal->r > TIMER_SIEVE) {
        time_sieve = time_sieve + seconds() - sec;
      }
#endif // TIMER_SIEVE

#ifdef NUMBER_HIT
      printf("Number of hits: %" PRIu64 " for r: %" PRIu64 "\n", number_of_hit,
             fb->factor_base_1[i]->ideal->r);
      printf("Estimated number of hits: %u.\n",
             (unsigned int) ceil((double) array->number_element /
                  (double) fb->factor_base_1[i]->ideal->r));
      number_of_hit = 0;
#endif // NUMBER_HIT

    } else {
      unsigned int index = 1;
      while (Tqr[index] == 0 && index < H->t) {
        index++;
      }

      if (index == H->t - 1) {
        uint64_t number_c_u = array->number_element;
        uint64_t number_c_l = 1;
        for (uint64_t j = 0; j < index; j++) {
          number_c_u = number_c_u / (2 * H->h[j] + 1);
          number_c_l = number_c_l * (2 * H->h[j] + 1);
        }

	sieve_ci(array, c, fb->factor_base_1[i], 0, H, index, number_c_l);
      } else {
        int64_t ci = 0;
        if (index + 1 < H->t - 1) {
          int64_vector_setcoordinate(c, index + 1, c->c[index + 1] - 1);
        } else if (index + 1 == H->t - 1) {
          int64_vector_setcoordinate(c, index + 1, -1);
        }

        if (index < H->t - 1) {
          for (unsigned int j = index + 1; j < H->t; j++) {
            ci = ci + (int64_t)Tqr[j] * c->c[j];
            if (ci >= (int64_t)fb->factor_base_1[i]->ideal->r || ci < 0) {
              ci = ci % (int64_t)fb->factor_base_1[i]->ideal->r;
            }
          }
          if (ci < 0) {
            ci = ci + (int64_t)fb->factor_base_1[i]->ideal->r;
          }
          ASSERT(ci >= 0);
        }

        uint64_t number_c_u = array->number_element;
        uint64_t number_c_l = 1;
        for (uint64_t j = 0; j < index; j++) {
          number_c_u = number_c_u / (2 * H->h[j] + 1);
          number_c_l = number_c_l * (2 * H->h[j] + 1);
        }
        number_c_u = number_c_u / (2 * H->h[index] + 1);

	unsigned int pos = 0;
        pos = int64_vector_add_one_i(c, index + 1, H);
        sieve_1(array, c, Tqr, fb->factor_base_1[i], H, index,
                number_c_l, &ci, pos);
        for (uint64_t j = 1; j < number_c_u; j++) {
          pos = int64_vector_add_one_i(c, index + 1, H);
          sieve_1(array, c, Tqr, fb->factor_base_1[i], H, index,
                  number_c_l, &ci, pos);
        }
      }
    }

/* #ifdef PROJECTIVE_ROOT */
/*     mpz_t mod; */
/*     mpz_init(mod); */
/*     mpz_mod_ui(mod, mpz_poly_lc_const(f), fb->factor_base_1[i]->ideal->r); */
/*     if (mpz_cmp_ui(mod, 0) == 0) { */
/*       mpz_t g0; */
/*       mpz_init(g0); */
/*       mpz_set(g0, fb->factor_base_1[i]->ideal->g->coeff[0]); */
/*       /\* WARNING: if g0 is equal to 0, the projective root maps to infinity, *\/ */
/*       /\* therefore the polynomial do not exist. *\/ */
/*       if (mpz_cmp_ui(g0, 0) != 0) { */
/*         mpz_t q; */
/*         mpz_init(q); */
/*         mpz_set_ui(r, special_q->ideal->r); */
/*         mpz_invert(g0, g0, q); */

/*         //Build Tqr */
/*         mpz_t * Tqr = (mpz_t *) malloc(sizeof(mpz_t) * (H->t)); */
/*         mpz_set(Tqr[0], q); */
/*         mpz_set(Tqr[1], g0); */
/*         mpz_mul_si(g0, -1); */
/*         for (unsigned int j = 2; j < H->t; j++) { */
/*           mpz_mul(Tqr[j], g0, Tqr[j - 1]); */
/*           mpz_mod_ui(Tqr[j], Tqr[j], ideal->ideal->r); */
/*         } */
/*         mpz_clear(q); */

/*         //TODO: to be continued. */

/*       } */

/*       mpz_clear(g0); */
/*     } */
/*     mpz_clear(mod); */
/* #endif */

    int64_vector_clear(c);
    free(Tqr);
  }

#ifdef TIMER_SIEVE
  printf("Bound: %d\n", TIMER_SIEVE);
  printf("Build Tqr: %f\n", time_tqr);
  printf("Perform sieve: %f\n", time_sieve);
#endif // TIMER_SIEVE

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
}

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

void printf_relation(factor_ptr factor0, factor_ptr factor1, mpz_poly_srcptr a,
                     int t)
{
  printf("# ");
  mpz_poly_fprintf(stdout, a);
  for (int i = 0; i < a->deg; i++) {
    gmp_printf("%Zd,", a->coeff[i]);
  }
  if (t - 1 == a->deg) {
    gmp_printf("%Zd:", a->coeff[a->deg]);
  } else {
    gmp_printf("%Zd,", a->coeff[a->deg]);
    for (int i = a->deg + 1; i < t - 1; i++) {
      printf("0,");
    }
    printf("0:");
  }
  for (unsigned int i = 0; i < factor0->number - 1; i++) {
    gmp_printf("%Zd,", factor0->factorization[i]);
  }
  gmp_printf("%Zd:", factor0->factorization[factor0->number - 1]);
  for (unsigned int i = 0; i < factor1->number - 1; i++) {
    gmp_printf("%Zd,", factor1->factorization[i]);
  }
  gmp_printf("%Zd\n", factor1->factorization[factor1->number - 1]);
}

void good_polynomial(mpz_poly_srcptr a, mpz_poly_t * f, mpz_t * lpb, int t)
{
  mpz_t res0;
  mpz_t res1;
  mpz_init(res0);
  mpz_init(res1);

  mpz_poly_resultant(res0, f[0], a);
  mpz_poly_resultant(res1, f[1], a);
  mpz_abs(res0, res0);
  mpz_abs(res1, res1);

  factor_t factor0;
  factor_t factor1;

  gmp_factorize(factor0, res0);
  gmp_factorize(factor1, res1);

  if (factor_is_smooth(factor0, lpb[0]) && factor_is_smooth(factor1, lpb[1])) {
    printf_relation(factor0, factor1, a, t);
  }
  factor_clear(factor0);
  factor_clear(factor1);
  mpz_clear(res0);
  mpz_clear(res1);
}

void find_relation(uint64_array_t * indexes, uint64_t number_element, mpz_t *
                   lpb, mat_Z_srcptr matrix, mpz_poly_t * f,
                   sieving_interval_srcptr H)
{
  if (0 != MIN(indexes[0]->length, indexes[1]->length)) {
    //Verify that if there is an element in both indexes. If not, do other thing.
    uint64_t stop = MAX(indexes[0]->length, indexes[1]->length);
    //Master side.
    int mside = 0;
    //Other side.
    int oside = 1;
    if (stop == indexes[1]->length) {
      mside = 1;
      oside = 0;
    }

    uint64_t j = 0;

    for (uint64_t i = 0; i < stop; i++) {
      if(indexes[oside]->array[j] < indexes[mside]->array[i]) {
        if (j < indexes[oside]->length - 1) {
          j++;
        }
      } else if(indexes[oside]->array[j] == indexes[mside]->array[i]) {

        mpz_vector_t c;
        mpz_t gcd;
        mpz_init(gcd);
        mpz_vector_init(c, H->t);
        array_index_mpz_vector(c, indexes[oside]->array[j], H,
                               number_element);

        mpz_poly_t a;
        mpz_poly_init(a, 0);
        mat_Z_mul_mpz_vector_to_mpz_poly(a, matrix, c);
        mpz_poly_content(gcd, a);

        //a must be irreducible.
        if (mpz_cmp_ui(gcd, 1) == 0 && a->deg > 0 &&
            mpz_cmp_ui(mpz_poly_lc_const(a), 0) > 0) {
          good_polynomial(a, f, lpb, (int)H->t);
        }

        mpz_poly_clear(a);

        mpz_clear(gcd);
        mpz_vector_clear(c);

        if (j < indexes[oside]->length - 1) {
          j++;
        }
      }
    }

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
}

/*
  Initialise the parameters of the special-q sieve.
  f: the V functions to define the number fields.
  fbb: the V factor base bounds.
  t: dimension of the lattice.
  H: sieving interval.
  q_min: lower bound of the special-q range.
  q_max: upper bound of the special-q range.
  thresh: the V threshold.
  lpb: the V large prime bounds.
  array: array in which the norms are stored.
  matrix: the Mq matrix (set with zero coefficients).
  q_side: side of the special-q.
  V: number of number fields.
 */

void initialise_parameters(int argc, char * argv[], mpz_poly_t ** f,
                           uint64_t ** fbb, factor_base_t ** fb,
                           sieving_interval_ptr H,
                           uint64_t * q_min, uint64_t * q_max,
                           unsigned char ** thresh, mpz_t ** lpb,
                           array_ptr array, mat_Z_ptr matrix,
                           unsigned int * q_side , unsigned int * V)
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
  ASSERT(* V >= 2 && * V < 10);

  unsigned int t;
  int * r;
  param_list_parse_int_list_size(pl, "H", &r, &t, ".,");
  ASSERT(t >= 2);
  sieving_interval_init(H, t);
  for (unsigned int i = 0; i < t; i++) {
    sieving_interval_set_hi(H, i, (unsigned int) r[i]);
  }
  free(r);

  * fbb = malloc(sizeof(uint64_t * ) * (* V));
  * thresh = malloc(sizeof(unsigned char *) * (* V));
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
    param_list_parse_uint(pl, str, (unsigned int *) &thresh[0][i]);
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

  uint64_t number_element = 0;
  sieving_interval_number_element(&number_element, H);
  ASSERT(number_element >= 6);

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
  sieving_interval_t H;
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
  printf("H =\n");
  sieving_interval_fprintf(stdout, H);
  printf("V = %u\n", V);
  for (unsigned int i = 0; i < V; i++) {
    printf("fbb%u = %" PRIu64 "\n", i, fbb[i]);
  }
  for (unsigned int i = 0; i < V; i++) {
    printf("thresh%u = %u\n", i, (unsigned int)thresh[i]);
  }
  for (unsigned int i = 0; i < V; i++) {
    gmp_printf("lpb%u = %Zd\n", i, lpb[i]);
  }
  for (unsigned int i = 0; i < V; i++) {
    printf("f%u = ", i);
    mpz_poly_fprintf(stdout, f[i]);
  }
  printf("q_min = %" PRIu64 "\n", q_min);
  printf("q_max = %" PRIu64 "\n", q_min);
  printf("q_side = %u\n", q_side);
#endif // PRINT_PARAMETERS

  uint64_array_t * indexes = malloc(sizeof(uint64_array_t) * V);

  double ** pre_compute = malloc(sizeof(double * ) * V);
  for (unsigned int i = 0; i < V; i++) {
    pre_compute[i] = malloc((H->t) * sizeof(double));
    pre_computation(pre_compute[i], f[i], H->t);
  }

  double ** time = malloc(sizeof(double * ) * V);
  for (unsigned int i = 0; i < V; i++) {
    time[i] = malloc(sizeof(double) * 3);
  }
  double sec_tot;
  double sec_cofact;

#ifdef TRACE_POS
  file = fopen("TRACE_POS.txt", "w+");
  fprintf(file, "TRACE_POS: %d\n", TRACE_POS);
#endif

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
#endif

        sec = seconds();
        sec_tot = sec;

        //LLL part
        build_Mq_ideal_1(matrix, special_q);

#ifdef TRACE_POS
        fprintf(file, "Mq:\n");
        mat_Z_fprintf(file, matrix);
#endif

        mat_Z_LLL_transpose(matrix);

#ifdef TRACE_POS
        fprintf(file, "MqLLL:\n");
        mat_Z_fprintf(file, matrix);
#endif

#ifdef MEAN_NORM_BOUND
        double * norms_bound = malloc(sizeof(double) * V);
#endif // MEAN_NORM_BOUND

#ifdef MEAN_NORM
        double * norms = malloc(sizeof(double) * V);
#endif // MEAN_NORM

#ifdef MEAN_NORM_BOUND_SIEVE
        double * norms_bound_sieve = malloc(sizeof(double) * V);
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
          special_q_sieve(array, matrix, fb[j], H);
          time[j][1] = seconds() - sec;
          sec = seconds();
          /* find_index(indexes[j], array, thresh[j]); */

#ifdef MEAN_NORM_BOUND_SIEVE
        norms_bound_sieve[j] = norm_bound_sieve;
#endif // MEAN_NORM_BOUND_SIEVE

          time[j][2] = seconds() - sec;
          sec = seconds();

#ifdef TRACE_POS
            fprintf(file, "********************\n");
#endif
        }

        sec = seconds();
        //TODO: do the MNFS way fo find_relation.
        /* find_relation(indexes, array->number_element, lpb, matrix, f, H); */
        sec_cofact = seconds() - sec;

        for (unsigned j = 0; j < V; j++) {
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
        }
        printf("# Time to factorize: %fs.\n", sec_cofact);
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
#endif

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
#endif

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
  sieving_interval_clear(H);
  free(fbb);
  free(thresh);

  return 0;
}
