#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include "utils_norm.h"
#include "mat_int64.h"

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
 * q: value of the special-q (ie, q^(deg(g))), 1.0 if no special-q.
 */
static void mean_norm(mpz_poly_srcptr f, int64_poly_srcptr a, double q)
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
static void trace_pos_init(uint64_t i, int64_vector_srcptr vector, int64_poly_srcptr a,
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

#if 0
/*
 * Contains all the mode we can activate during the initialisation of the norms. *
 * special_q: 1 if there is a special_q, 0 otherwise.
 * f: the polynomial that defines the number field.
 * a: the polynomial.
 * q: the q of the special-q.
 * log_bound_resulant: log2 of an approximation of the resultant between a
 *  and f.
 * i: position in the array we want to follow.
 * vector: the vector corresponding to the ith position.
 * array: the array in which we store the norm.
 * deg_g: the degree of the g of the special-q.
 */
static void mode_init_norm(MAYBE_UNUSED int special_q,
    MAYBE_UNUSED mpz_poly_srcptr f, MAYBE_UNUSED int64_poly_srcptr a,
    MAYBE_UNUSED ideal_spq_srcptr spq, MAYBE_UNUSED double log_bound_resultant,
    MAYBE_UNUSED uint64_t i, MAYBE_UNUSED int64_vector_srcptr vector,
    MAYBE_UNUSED array_srcptr array)
{
#ifdef MEAN_NORM
  if (special_q == 1) {
    mean_norm(f, a, pow(((double) ideal_spq_get_q(spq), (double)
            ideal_spq_get_deg_h(spq)));
  } else {
    mean_norm(f, a, 1.0);
  }
#endif // MEAN_NORM

#ifdef MEAN_NORM_BOUND
  if (special_q == 1) {
    norm_bound = norm_bound + pow(2, log_bound_resultant) /
    pow(((double) ideal_spq_get_q(spq), (double) ideal_spq_get_deg_h(spq)));
  } else {
    norm_bound = norm_bound + pow(2, log_bound_resultant);
  }
#endif // MEAN_NORM_BOUND

#ifdef TRACE_POS
  trace_pos_init(i, vector, a, f, bound_resultant, special_q, array);
#endif // TRACE_POS
}
#endif // 0

#ifdef ASSERT_NORM
static unsigned char log_norm(const int * current_indexes, mat_int64_srcptr M,
    double_poly_srcptr f, ideal_spq_srcptr spq, int special_q,
    MAYBE_UNUSED unsigned int size)
{
  ASSERT(size == M->NumRows);
  ASSERT(size == M->NumCols);

  mpz_poly_t poly;
  mpz_poly_init(poly, M->NumRows);

  for (unsigned int i = 0; i < M->NumRows; i++) {
    int64_t tmp = 0;
    for (unsigned int k = 0; k < M->NumCols; k++) {
      tmp = tmp + M->coeff[i + 1][k + 1] * (int64_t) current_indexes[k];
    }
    mpz_poly_setcoeff_int64(poly, i, tmp);
  }

  mpz_poly_t f_Z;
  mpz_poly_init(f_Z, f->deg);
  mpz_poly_set_double_poly(f_Z, f);
  mpz_t norm;
  mpz_init(norm);
  norm_poly(norm, poly, f_Z);
  mpz_poly_clear(poly);
  mpz_poly_clear(f_Z);

  unsigned char res = (unsigned char) round(log2(mpz_get_d(norm)));

  if (special_q) {
    ASSERT(special_q == 1);

    res = res - ideal_spq_get_log(spq);
  }

  mpz_clear(norm);
  return res;
}

void assert_norm(array_srcptr array, sieving_bound_srcptr H, mpz_poly_srcptr f,
    mat_Z_srcptr matrix)
{
  mpz_vector_t c;
  mpz_vector_init(c, H->t);
  mpz_poly_t a;
  mpz_poly_init(a, 0);

  mpz_t norm;
  mpz_init(norm);

  for (uint64_t i = 0; i < array->number_element; i++) {
    array_index_mpz_vector(c, i, H, array->number_element);
    mat_Z_mul_mpz_vector_to_mpz_poly(a, matrix, c);
    norm_poly(norm, f, a);

    double log_norm = (double)mpz_sizeinbase(norm, 2);

    if (fabs((double)array->array[i] - log_norm) > 2.0) {
      fprintf(stderr, "# Error in norm at index %" PRIu64 "\n", i);
      fprintf(stderr, "# Vector c: ");
      mpz_vector_fprintf(stderr, c);
      fprintf(stderr, "# Polynomial a: ");
      mpz_poly_fprintf(stderr, a);
      fprintf(stderr, "# Value in the array: %u. Log of the norm: %f\n",
          array->array[i], log_norm);
    }
  }
  mpz_vector_clear(c);
  mpz_poly_clear(a);
  mpz_clear(norm);
}
#endif // ASSERT_NORM

#ifdef OLD_NORM
void pre_computation(double * pre_compute, mpz_poly_srcptr f, unsigned int d)
{
  pre_compute[0] = 0;
  mpz_t infinity_norm_f_tmp;
  mpz_init(infinity_norm_f_tmp);
  //Infinite norm of f.
  mpz_poly_infinity_norm(infinity_norm_f_tmp, f);
  double infinity_norm_f = mpz_get_d(infinity_norm_f_tmp);
  mpz_clear(infinity_norm_f_tmp);
  //Do the computation for each i.
  for (unsigned int i = 1; i < d; i++) {
    pre_compute[i] = (double)i * log2(infinity_norm_f) + ((double)f->deg / 2) *
      log2((double)(i + 1)) + ((double)i / 2) * log2((double)(f->deg + 1));
  }
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
void init_each_cell(array_ptr array, uint64_t i, int64_poly_ptr a,
    double * pre_compute, sieving_bound_srcptr H, mat_int64_srcptr matrix,
    unsigned int beg, mpz_poly_srcptr f, ideal_spq_srcptr spq, int special_q,
    MAYBE_UNUSED int64_vector_srcptr vector)
{
  double log_bound_resultant = 0;
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
  ASSERT(int64_poly_equal(a, a_tmp));
  int64_poly_clear(a_tmp);
#endif // NDEBUG

  //Compute an approximation of the norm.
  if (a->deg > 0) {
    uint64_t tmp = int64_poly_infinity_norm(a);
    log_bound_resultant =
      (double)pre_compute[a->deg] + log2((double)tmp) * (double)f->deg;
    
  } else {
    log_bound_resultant = 0;
  }
  /*
   * Store in array the value of the norm, without the contribution of the
   *  special-q if it is set in this number field.
   */
  if (special_q) {
    ASSERT(special_q == 1);

    array->array[i] = (unsigned char)log_bound_resultant -
      ideal_spq_get_log(spq);
    mode_init_norm(special_q, f, a, spq, log_bound_resultant, i,
        vector, array);
  } else {
    array->array[i] = (unsigned char)log_bound_resultant;
    mode_init_norm(special_q, f, a, spq, log_bound_resultant, i, vector,
        array);
  }
}

#else // OLD_NORM

static void next_bottom_left_cube_in_hypercube(int * bottom_left_cube,
    const unsigned int * length, const unsigned int * use_length,
    const int * bottom_left_hypercube, const unsigned int * length_hypercube,
    unsigned int t)
{
#ifndef NDEBUG
  int tmp = 0;
  int tmp2 = 0;
  for (unsigned int j = 0; j < t; j++) {
    tmp = tmp + bottom_left_cube[j];
    tmp2 = tmp2 + bottom_left_hypercube[j] - (int)length[j] +
      (int)length_hypercube[j];
  }
  ASSERT(tmp <= tmp2);
#endif

  unsigned int k = 0;
  while(k < t) {
    if (bottom_left_cube[k] + (int)length[k] == bottom_left_hypercube[k] +
        (int)length_hypercube[k]) {
      bottom_left_cube[k] = bottom_left_hypercube[k];
      k++;
    } else {
      break;
    }
  }
  if (k < t && use_length[k]) {
    bottom_left_cube[k] = bottom_left_cube[k] + (int)length[k];
  }

#ifndef NDEBUG
  for (unsigned int j = 0; j < t - 1; j++) {
    tmp = bottom_left_cube[j];
    ASSERT(tmp >= bottom_left_hypercube[j]);
    ASSERT(tmp + (int)length[j] - 1 < bottom_left_hypercube[j] +
        (int)length_hypercube[j]);
  }
  tmp = bottom_left_cube[t - 1];
  ASSERT(tmp >= 0);
  ASSERT(tmp + (int)length[t - 1] - 1 < bottom_left_hypercube[t - 1] +
      (int)length_hypercube[t - 1]);
#endif
}

void add_boolean(int * tab, unsigned int size)
{
#ifndef NDEBUG
  for (unsigned int j = 0; j < size; j++) {
    ASSERT(tab[j] < 2);
    ASSERT(tab[j] > -2);
  }
#endif

  unsigned int k = 0;
  while(k < size) {
    if (tab[k] == 1) {
      tab[k] = 0;
      k++;
    } else {
      break;
    }
  }
  if (k < size) {
    tab[k]++;
  }

#ifndef NDEBUG
  for (unsigned int j = 0; j < size; j++) {
    ASSERT(tab[j] < 2);
    ASSERT(tab[j] > -1);
  }
#endif
}

void add_arrays(int * current_indexes, const int * current_values,
    const int * true_false, const unsigned int * increment, unsigned int size)
{
  for (unsigned int i = 0; i < size; i++) {
    current_indexes[i] = current_values[i] + (unsigned int) true_false[i] *
      (increment[i] - 1);
  }
}

unsigned char log_norm_double(const int * current_indexes, mat_int64_srcptr M,
    double_poly_srcptr f, ideal_spq_srcptr spq, int special_q,
    MAYBE_UNUSED unsigned int size)
{
  ASSERT(size == M->NumRows);
  ASSERT(size == M->NumCols);

  double_poly_t poly;
  double_poly_init(poly, M->NumRows);

  int deg = -1;

  for (unsigned int i = 0; i < M->NumRows; i++) {
    int64_t tmp = 0;
    for (unsigned int k = 0; k < M->NumCols; k++) {
      tmp = tmp + M->coeff[i + 1][k + 1] * (int64_t) current_indexes[k];
    }
    poly->coeff[i] = (double) tmp;
    if (tmp != 0) {
      deg = (int) i;
    }
  }
  poly->deg = deg;

  double resultant = fabs(double_poly_resultant(poly, f));
  double_poly_clear(poly);

  ASSERT(resultant >= 0.0);

  unsigned char res = (unsigned char) round(log2(resultant));

  if (special_q) {
    ASSERT(special_q == 1);

    res = res - ideal_spq_get_log(spq);
  }

  return res;
}

void update_mini_maxi(unsigned char * mini, unsigned char * maxi,
    unsigned char val)
{
  if (val < * mini) {
    * mini = val;
  }
  if (val > * maxi) {
    * maxi = val;
  }
}

int generate_random_number(int a, int b)
{
  int random = rand() % (b - a) + a;

  ASSERT(random >= a);
  ASSERT(random < b);

  return random;
}

void add_one_values(int * current_indexes, const int * bottom_left_cube,
    const unsigned int * increment, unsigned int size)
{
#ifndef NDEBUG
  int tmp = 0;
  int tmp2 = 0;
  for (unsigned int j = 0; j < size; j++) {
    tmp = tmp + current_indexes[j];
    tmp2 = tmp2 + bottom_left_cube[j] + (int)increment[j] - 1;
  }
  ASSERT(tmp <= tmp2);
#endif

  unsigned int k = 0;
  while(k < size) {
    if (current_indexes[k] == bottom_left_cube[k] + (int)increment[k] - 1) {
      current_indexes[k] = bottom_left_cube[k];
      k++;
    } else {
      break;
    }
  }
  if (k < size) {
    current_indexes[k]++;
  }

#ifndef NDEBUG
  for (unsigned int j = 0; j < size - 1; j++) {
    tmp = current_indexes[j];
    ASSERT(tmp >= bottom_left_cube[j]);
    ASSERT(tmp < bottom_left_cube[j] + (int) increment[j]);
  }
  tmp = current_indexes[size - 1];
  ASSERT(tmp >= bottom_left_cube[size - 1]);
  ASSERT(tmp < bottom_left_cube[size - 1] + (int) increment[size - 1]);
#endif
}

static void init_cells(array_ptr array, const int * bottom_left_cube,
    const unsigned int * length, sieving_bound_srcptr H, double_poly_srcptr f,
    mat_int64_srcptr Mq, ideal_spq_srcptr spq, int special_q)
{
#ifndef NDEBUG
  for (unsigned int i = 0; i < H->t; i++) {
    ASSERT(length[i] >= 2);
  }
#endif // NDEBUG


  int * true_false = (int *) malloc(sizeof(int) * H->t);
  memset(true_false, 0, sizeof(int) * H->t);
  true_false[0] = -1;

  int * current_indexes = (int *) malloc(sizeof(int) * H->t);
  memset(current_indexes, 0, sizeof(int) * H->t);

  //Compute the number of vertices.
  uint64_t stop = 1 << H->t;

  unsigned char tmp = 0;
  unsigned char maxi = 0;
  unsigned char mini = UCHAR_MAX;

  //Compute each norm of the vertices of the cube.
  for (uint64_t i = 0; i < stop; i++) {
    add_boolean(true_false, H->t);
    add_arrays(current_indexes, bottom_left_cube, true_false, length, H->t);

    //Get the possible value of current_indexes.
    tmp = array_get_at(array, current_indexes, H);
    if (tmp == 255) {
      tmp = log_norm_double(current_indexes, Mq, f, spq, special_q, H->t);

      array_set_at(array, current_indexes, tmp, H);

    } 
    update_mini_maxi(&mini, &maxi, tmp);
  }
  free(true_false);

  unsigned int * new_length = (unsigned int * ) malloc(sizeof(unsigned int) *
      H->t);
  memset(new_length, 0, sizeof(unsigned int) * H->t);
  unsigned int * use_new_length = (unsigned int * ) malloc(
      sizeof(unsigned int) * H->t);
  memset(use_new_length, 0, sizeof(unsigned int) * H->t);
  unsigned int all_length_2 = 0;

  for (unsigned int i = 0; i < H->t; i++) {
    //Generate a random element inside the cube.
    current_indexes[i] = generate_random_number(bottom_left_cube[i],
        bottom_left_cube[i] + length[i]);
    if (length[i] > 2) {
      new_length[i] = length[i] / 2;
      use_new_length[i] = 1;
    } else {
      ASSERT(length[i] <= 2);

      new_length[i] = 2;
      all_length_2++;
      use_new_length[i] = 0;
    }
  }


  //If all_length_2 == H->t, all the cells are initialized.
  if (all_length_2 != H->t) {

    tmp = array_get_at(array, current_indexes, H);
    if (tmp == 255) {
      tmp = log_norm_double(current_indexes, Mq, f, spq, special_q, H->t);

      array_set_at(array, current_indexes, tmp, H);

    } 
    update_mini_maxi(&mini, &maxi, tmp);

    //TODO: arbitrary value. If maxi - mini <= 2, all the cells are intialzed to
    //maxi.
    if (maxi - mini > 2) {
      unsigned int i = 0;
      while (new_length[i] == 2) {
        i++;
      }
      uint64_t stop = 1 << (H->t - all_length_2);

      for (i = 0; i < H->t; i++) {
        current_indexes[i] = bottom_left_cube[i];
      }
      current_indexes[0] = current_indexes[0] - new_length[0];
      for (i = 0; i < stop; i++) {
        //TODO: not H but other think (ie length)
        next_bottom_left_cube_in_hypercube(current_indexes, new_length,
            use_new_length, bottom_left_cube, length, H->t);
        init_cells(array, current_indexes, new_length, H, f, Mq, spq,
            special_q);
      }
    } else {
      stop = 1;
      for (unsigned int i = 0; i < H->t; i++) {
        current_indexes[i] = bottom_left_cube[i];
        stop = stop * (uint64_t)length[i];
      }
      current_indexes[0]--;

      //TODO: use memset for the first coordinate.
      for (uint64_t i = 0; i < stop; i++) {
        add_one_values(current_indexes, bottom_left_cube, length, H->t);
        tmp = array_get_at(array, current_indexes, H);
        if (tmp == 255) {
          array_set_at(array, current_indexes, maxi, H);
        } 
      }
    }
  }
  free(current_indexes);
  free(new_length);
  free(use_new_length);
}
#endif // OLD_NORM

#ifdef OLD_NORM
void init_norm(array_ptr array, double * pre_compute,
    sieving_bound_srcptr H, mat_Z_srcptr matrix, mpz_poly_srcptr f,
    ideal_spq_srcptr spq, int special_q)
#else // OLD_NORM
void init_norm(array_ptr array, sieving_bound_srcptr H, mat_Z_srcptr matrix,
    mpz_poly_srcptr f, ideal_spq_srcptr spq, int special_q)
#endif // OLD_NORM
{
  ASSERT(special_q == 0 || special_q == 1);

  //TODO: error if H = [1, 1, …]
#ifndef OLD_NORM
  srand(time(NULL));

  /*Coordinate of the hypercube*/
  int * bottom_left_hypercube = (int *) malloc(sizeof(int) * H->t);
  unsigned int * length_hypercube = (unsigned int *)
    malloc(sizeof(unsigned int) * H->t);

  /*First division of the hypercube*/
  //Length of the first division of the sieving region in hypercubes.
  unsigned int * length = (unsigned int * ) malloc(sizeof(unsigned int)
      * H->t);
  //Say if the length must be used, ie if the corresponding value in length is
  //not equal to 2 for the first time.
  unsigned int * use_length = (unsigned int * ) malloc(sizeof(unsigned int)
      * H->t);
  //Coordinate of the an extremal point. In 2 dimension, this point will be the
  //point at the bottom left of a square.
  int * bottom_left_cube = (int * ) malloc(sizeof(int) * H->t);

  //TODO: why uint64_t?
  uint64_t stop = (uint64_t) 1 << H->t;

  //TODO: eventually, think about H->h[i] == 2.
  //TODO: arbitrary value. The sieving bound gives 2^H->t smaller hypercubes.
  for (uint64_t i = 0; i < H->t - 1; i++) {
    length[i] = H->h[i];
    use_length[i] = 1;
    bottom_left_cube[i] = -(int) H->h[i];
    bottom_left_hypercube[i] = -(int) H->h[i];
    length_hypercube[i] = 2 * H->h[i];
  }
  length[H->t - 1] = H->h[H->t - 1] / 2;
  use_length[H->t - 1] = 1;
  bottom_left_cube[H->t - 1] = 0;
  bottom_left_hypercube[H->t - 1] = 0;
  length_hypercube[H->t - 1] = H->h[H->t - 1];

  bottom_left_cube[0] = -(int) H->h[0] - (int)length[0];
  bottom_left_cube[H->t - 1] = 0;

  double_poly_t f_d;
  double_poly_init(f_d, f->deg);
  double_poly_set_const_mpz_poly(f_d, f);

  mat_int64_t Mq;
  mat_int64_init(Mq, matrix->NumRows, matrix->NumCols);
  mat_Z_to_mat_int64(Mq, matrix);  

  //Iterate on all the hypercube.
  for (unsigned int i = 0; i < stop; i++) {
    next_bottom_left_cube_in_hypercube(bottom_left_cube, length, use_length,
        bottom_left_hypercube, length_hypercube, H->t);
    init_cells(array, bottom_left_cube, length, H, f_d, Mq, spq, special_q);
  }

  mat_int64_clear(Mq);
  double_poly_clear(f_d);
  free(bottom_left_cube);
  free(length);
  free(use_length);
  free(bottom_left_hypercube);
  free(length_hypercube);

#else // OLD_NORM

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
  init_each_cell(array, 0, a, pre_compute, H, matrix_int, beg, f, spq,
      special_q, vector);

  for (uint64_t i = 1; i < array->number_element; i++) {
    beg = int64_vector_add_one_i(vector, 0, H);
    init_each_cell(array, i, a, pre_compute, H, matrix_int, beg, f, spq,
        special_q, vector);
  }

  int64_poly_clear(a);
  int64_vector_clear(vector);
  mat_int64_clear(matrix_int);
#endif // OLD_NORM
}
