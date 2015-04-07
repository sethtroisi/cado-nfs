#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include "utils_lattice.h"
#include "ideal.h"
#include "utils_mpz.h"
#include "sieving_bound.h"
#include "mat_double.h"
#include "double_vector.h"

#define DEFAULT_LENGTH_LIST_VECTOR 3
/*
 * List of int64_vector
 */

void list_int64_vector_init(list_int64_vector_ptr list)
{
  list->length = 0;
  list->v = (int64_vector_t * ) malloc(sizeof(int64_vector_t) *
      DEFAULT_LENGTH_LIST_VECTOR);
}

void list_int64_vector_add_int64_vector(list_int64_vector_ptr list, int64_vector_srcptr v)
{
  if ((list->length % DEFAULT_LENGTH_LIST_VECTOR) == 0 && list->length != 0) {
  list->v = realloc(list->v, sizeof(int64_vector_t) * (list->length +
        DEFAULT_LENGTH_LIST_VECTOR));
  }
  int64_vector_init(list->v[list->length], v->dim);
  int64_vector_set(list->v[list->length], v);
  list->length++;
}

void list_int64_vector_clear(list_int64_vector_ptr list)
{
  for (unsigned int i = 0; i < list->length; i++) {
    ASSERT(list->v[i]->dim != 0);
    int64_vector_clear(list->v[i]);
  }
  free(list->v);
  list->length = 0;
}

void list_int64_vector_fprintf(FILE * file, list_int64_vector_srcptr list)
{
  fprintf(file, "[\n");
  for (unsigned int i = 0; i < list->length - 1; i++) {
    int64_vector_fprintf(file, list->v[i]);
  }
  if (list->length != 0) {
    int64_vector_fprintf(file, list->v[list->length - 1]);
  }
  fprintf(file, "]\n");
}

void list_int64_vector_extract_mat_int64(list_int64_vector_ptr list,
    mat_int64_srcptr matrix)
{
  ASSERT(list->length == 0);

  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, matrix->NumRows);

  for (unsigned int i = 0; i < matrix->NumCols; i++) {
    mat_int64_extract_vector(v_tmp, matrix, i);
    list_int64_vector_add_int64_vector(list, v_tmp);
  }

  int64_vector_clear(v_tmp);
}

/*
 * From PNpoly
 */
static int int64_vector_in_list_int64_vector(int64_vector_srcptr vec, list_int64_vector_srcptr list)
{
  ASSERT(list->length > 2);
  ASSERT(vec->dim > 1);
#ifndef NDEBUG
  for (unsigned int i = 0; i < list->length; i++) {
    ASSERT(vec->dim == list->v[i]->dim);
    for (unsigned int j = 2; j < vec->dim; j++) {
      ASSERT(vec->c[j] == list->v[i]->c[j]);
    }
  }
#endif // NDEBUG

  for (unsigned int i = 0; i < list->length; i++) {
    if (int64_vector_equal(vec, list->v[i]) == 1) {
      return 1;
    }  
  }

  int c = 0;
  unsigned i = 0;
  unsigned j = 0;
  for (i = 0, j = list->length - 1; i < list->length; j =
      i++) {
    if ( ((list->v[i]->c[1] > vec->c[1]) != (list->v[j]->c[1] > vec->c[1])) &&
        (vec->c[0] < (list->v[j]->c[0] - list->v[i]->c[0]) *
         (vec->c[1] - list->v[i]->c[1]) / (list->v[j]->c[1] - list->v[i]->c[1]) +
         list->v[i]->c[0]) ) {
      c = !c;
    }
  }

  return c;
}

/*
 * List of double vector.
 */
void list_double_vector_init(list_double_vector_ptr list)
{
  list->length = 0;
  list->v = (double_vector_t * ) malloc(sizeof(double_vector_t) *
      DEFAULT_LENGTH_LIST_VECTOR);
}

void list_double_vector_add_double_vector(list_double_vector_ptr list, double_vector_srcptr v)
{
  if ((list->length % DEFAULT_LENGTH_LIST_VECTOR) == 0 && list->length != 0) {
    list->v = realloc(list->v, sizeof(double_vector_t) * (list->length +
          DEFAULT_LENGTH_LIST_VECTOR));
  }
  double_vector_init(list->v[list->length], v->dim);
  double_vector_set(list->v[list->length], v);
  list->length++;
}

void list_double_vector_clear(list_double_vector_ptr list)
{
  for (unsigned int i = 0; i < list->length; i++) {
    ASSERT(list->v[i]->dim != 0);
    double_vector_clear(list->v[i]);
  }
  free(list->v);
  list->length = 0;
}

void list_double_vector_fprintf(FILE * file, list_double_vector_srcptr list)
{
  fprintf(file, "[\n");
  for (unsigned int i = 0; i < list->length - 1; i++) {
    double_vector_fprintf(file, list->v[i]);
  }
  if (list->length != 0) {
    double_vector_fprintf(file, list->v[list->length - 1]);
  }
  fprintf(file, "]\n");
}

void list_double_vector_extract_mat_int64(list_double_vector_ptr list,
    mat_int64_srcptr matrix)
{
  ASSERT(list->length == 0);

  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, matrix->NumRows);
  double_vector_t vd_tmp;
  double_vector_init(vd_tmp, v_tmp->dim);

  for (unsigned int i = 0; i < matrix->NumCols; i++) {
    mat_int64_extract_vector(v_tmp, matrix, i);
    int64_vector_to_double_vector(vd_tmp, v_tmp);
    list_double_vector_add_double_vector(list, vd_tmp);
  }

  int64_vector_clear(v_tmp);
  double_vector_clear(vd_tmp);
}

/*
 * Compute angle between two vectors.
 *
 * v0: first vector.
 * v1: second vector.
 */
static double angle_2_coordinate(int64_vector_srcptr v0, int64_vector_srcptr v1)
{
  ASSERT(v0->dim == v1->dim);
  ASSERT(v0->dim >= 2);

#ifndef NDEBUG
  for (unsigned int i = 2; i < v0->dim; i++) {
    ASSERT(v0->c[i] == v1->c[i]);
  } 
#endif // NDEBUG

  double dot = 0;
  for (unsigned int i = 0; i < 2; i++) {
    dot += (double)(v0->c[i] * v1->c[i]);
  }
  dot = dot / (sqrt((double)(v0->c[0] * v0->c[0] + v0->c[1] * v0->c[1])));
  dot = dot / (sqrt((double)(v1->c[0] * v1->c[0] + v1->c[1] * v1->c[1])));
  return acos(dot);
}

/*
 * Compute the determinant of a 2 * 2 matrix.
 */
static int determinant2(mat_int64_srcptr matrix)
{
  return (int) matrix->coeff[1][1] * matrix->coeff[2][2] - matrix->coeff[2][1] *
    matrix->coeff[1][2];
}

int gauss_reduction(int64_vector_ptr v0, int64_vector_ptr v1,
    mat_int64_ptr U, int64_vector_srcptr v0_root, int64_vector_srcptr v1_root)
{
  //Verify that the vector are in the same plane.
  ASSERT(v0_root->dim == v1_root->dim);
  ASSERT(v0_root->dim >= 2);
#ifndef NDEBUG
  for (unsigned int i = 3; i < v1->dim; i++) {
    ASSERT(v0_root->c[i] == v1_root->c[i]);
  }

  if (U != NULL) {
    ASSERT(U->NumRows == U->NumCols);
    ASSERT(U->NumRows == 2);
  }

  int64_vector_t v0_tmp;
  int64_vector_init(v0_tmp, v0->dim);
  int64_vector_set(v0_tmp, v0);

  int64_vector_t v1_tmp;
  int64_vector_init(v1_tmp, v1->dim);
  int64_vector_set(v1_tmp, v1);
#endif // NDEBUG
  ASSERT(v0->dim == v1->dim);
  ASSERT(v0_root->dim >= v1->dim);

  int64_vector_set(v0, v0_root);
  int64_vector_set(v1, v1_root);

  //determinant of U, 0 if U is NULL.
  int det = 0;
  
  double B0 = int64_vector_norml2sqr(v0);
  double m = (double)int64_vector_dot_product(v0, v1) / B0;
  for (unsigned int i = 0; i < 2; i++) {
    v1->c[i] = v1->c[i] - (int64_t)round(m) * v0->c[i];
  }

  if (U != NULL) {
    U->coeff[1][1] = 1;
    U->coeff[1][2] = -(int64_t)round(m);
    U->coeff[2][1] = 0;
    U->coeff[2][2] = 1;
    det = 1;
  }

  double B1 = int64_vector_norml2sqr(v1);

  while (B1 < B0) {
    int64_vector_swap(v0, v1);
    B0 = B1;
    m = (double)int64_vector_dot_product(v0, v1) / B0;
    for (unsigned int i = 0; i < 2; i++) {
      v1->c[i] = v1->c[i] - (int64_t)round(m) * v0->c[i];
    }

    //Update U by U = U * transformation_matrix.
    if (U != NULL) {
      int64_t tmp = U->coeff[1][2];
      U->coeff[1][2] = U->coeff[1][1] - (int64_t)round(m) * U->coeff[1][2];
      U->coeff[1][1] = tmp;
      tmp = U->coeff[2][2];
      U->coeff[2][2] = U->coeff[2][1] - (int64_t)round(m) * U->coeff[2][2];
      U->coeff[2][1] = tmp;
      det = -det;
    }

    B1 = int64_vector_norml2sqr(v1);

    //To have an acute basis.
  }

  if (M_PI_2 < angle_2_coordinate(v0, v1))
  {
    for (unsigned int i = 0; i < 2; i++) {
      v1->c[i] = -v1->c[i];
    }

    if (U != NULL) {
      U->coeff[1][2] = -U->coeff[1][2];
      U->coeff[2][2] = -U->coeff[2][2];
      det = -det;
    }
  }

#ifndef NDEBUG
  if (U != NULL) {
    ASSERT(v0_tmp->c[0] * U->coeff[1][1] + v1_tmp->c[0] * U->coeff[2][1]
        == v0->c[0]);
    ASSERT(v0_tmp->c[1] * U->coeff[1][1] + v1_tmp->c[1] * U->coeff[2][1]
        == v0->c[1]);
    ASSERT(v0_tmp->c[0] * U->coeff[1][2] + v1_tmp->c[0] * U->coeff[2][2]
        == v1->c[0]);
    ASSERT(v0_tmp->c[1] * U->coeff[1][2] + v1_tmp->c[1] * U->coeff[2][2]
        == v1->c[1]);
    ASSERT(determinant2(U) == det);
    ASSERT(ABS(det) == 1);
  }

  int64_vector_clear(v0_tmp);
  int64_vector_clear(v1_tmp);
#endif // NDEBUG

  return det;
}

int gauss_reduction_zero(int64_vector_ptr v0, int64_vector_ptr v1,
    mat_int64_ptr U, int64_vector_srcptr v0_root, int64_vector_srcptr v1_root)
{
  int det = gauss_reduction(v0, v1, U, v0_root, v1_root);

  for (unsigned int i = 3; i < v1->dim; i++) {
    v0->c[i] = 0;
    v1->c[i] = 0;
  }

  return det; 
}

//TODO: why not merge with gauss reduction?
int reduce_qlattice(int64_vector_ptr v0, int64_vector_ptr v1,
    int64_vector_srcptr v0_root, int64_vector_srcptr v1_root, int64_t I)
{
  ASSERT(v0_root->c[1] == 0);
  ASSERT(v1_root->c[1] == 1);
  ASSERT(v0_root->dim == v1_root->dim);
  ASSERT(v0->dim == v1->dim);
  ASSERT(v0_root->dim == v1->dim);
  ASSERT(I > 0);
#ifndef NDEBUG
  for (unsigned int i = 3; i < v1->dim; i++) {
    ASSERT(v0_root->c[i] == v1_root->c[i]);
  }
#endif // NDEBUG

  int64_vector_set(v0, v0_root);
  int64_vector_set(v1, v1_root);

  int64_t a0 = -v0->c[0], b0 = v1->c[0], a1 = 0, b1 = 1, k;

  const int64_t hI = I;
  const int64_t mhI = -hI;

  while (b0 >= hI) {
    k = a0 / b0; a0 %= b0; a1 -= k * b1;

    if (a0 > mhI) {
      break;
    }
    k = b0 / a0; b0 %= a0; b1 -= k * a1;

    if (b0 < hI) {
      break;
    }
    k = a0 / b0; a0 %= b0; a1 -= k * b1;

    if (a0 > mhI) {
      break;
    }
    k = b0 / a0; b0 %= a0; b1 -= k * a1;

  }
  k = b0 - hI - a0;
  if (b0 > -a0) {
    if (!a0) {
      return 0;
    }
    k /= a0; b0 -= k * a0; b1 -= k * a1;

  } else {
    if (!b0) {
      return 0;
    }
    k /= b0; a0 += k * b0; a1 += k * b1;

  }

  ASSERT(a0 > mhI);
  ASSERT(0 >= a0);
  ASSERT(0 <= b0);
  ASSERT(b0 < hI);
  ASSERT(a1 > 0);
  ASSERT(b1 > 0);

  v0->c[0] = a0;
  v0->c[1] = a1;
  v1->c[0] = b0;
  v1->c[1] = b1;

  return 1;
}

int reduce_qlattice_zero(int64_vector_ptr v0, int64_vector_ptr v1,
    int64_vector_srcptr v0_root, int64_vector_srcptr v1_root, int64_t I)
{
  int res = reduce_qlattice(v0, v1, v0_root, v1_root, I);

  for (unsigned int i = 3; i < v1->dim; i++) {
    v0->c[i] = 0;
    v1->c[i] = 0;
  }

  return res;
}

/*
 * For instance, it is SV3.
 */
void SV4(list_int64_vector_ptr SV, int64_vector_srcptr v0_root,
    int64_vector_srcptr v1_root, int64_vector_srcptr v2)
{
  ASSERT(v0_root->dim == v1_root->dim);
#ifndef NDEBUG
  for (unsigned int j = 3; j < v0_root->dim; j++) {
    ASSERT(v0_root->c[j] == v1_root->c[j]);
    ASSERT(v0_root->c[j] == 0);
  }
#endif // NDEBUG
  ASSERT(v2->c[1] == 0);
  ASSERT(v2->c[2] == 1);

  int64_vector_t u;
  int64_vector_init(u, v0_root->dim);

  int64_vector_t v0;
  int64_vector_t v1;
  int64_vector_init(v0, v0_root->dim);
  int64_vector_init(v1, v1_root->dim);
  int64_vector_set(v0, v0_root);
  int64_vector_set(v1, v1_root);

  mat_int64_t U;
  mat_int64_init(U, 2, 2);

  int det = gauss_reduction(v0, v1, U, v0, v1);
  ASSERT(det != 0);

  double target0 = (double)v2->c[0] / (double)v0_root->c[0];
  //U^-1 * [target0, 0] = (v2_new_base_x, v2_new_base_y).
  double v2_new_base_x = (double)det * (double)U->coeff[2][2] * target0;
  double v2_new_base_y = -(double)det * (double)U->coeff[2][1] * target0;
  int64_t a = (int64_t)round(v2_new_base_x);
  int64_t b = (int64_t)round(v2_new_base_y);

  for (unsigned int i = 0; i < 2; i++) {
    u->c[i] = a * v0->c[i] + b * v1->c[i];
  }
  u->c[2] = 0;
  list_int64_vector_add_int64_vector(SV, u);

  //Build a triangle around the projection of v2 in the plane z = 0.
  if (v2_new_base_x - (double)a < 0) {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] - v0->c[i];
    }
    u->c[2] = 0;
    list_int64_vector_add_int64_vector(SV, u);
  } else {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] + v0->c[i];
    }
    u->c[2] = 0;
    list_int64_vector_add_int64_vector(SV, u);
  }

  if (v2_new_base_y - (double)b < 0) {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] - v1->c[i];
    }
    u->c[2] = 0;
    list_int64_vector_add_int64_vector(SV, u);
  } else {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] + v1->c[i];
    }
    u->c[2] = 0;
    list_int64_vector_add_int64_vector(SV, u);
  }

#ifndef NDEBUG
  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, v2->dim);
  int64_vector_set(v_tmp, v2);
  for (unsigned int i = 2; i < v_tmp->dim; i++) {
    v_tmp->c[i] = 0;
  }
  ASSERT(int64_vector_in_list_int64_vector(v_tmp, SV) == 1);
  int64_vector_clear(v_tmp);
#endif // NDEBUG

  for (unsigned int i = 0; i < 3; i++) {
    int64_vector_sub(SV->v[i], v2, SV->v[i]);
  }

  int64_vector_clear(v0);
  int64_vector_clear(v1);
  int64_vector_clear(u);
  mat_int64_clear(U);
}

void SV4_Mqr(list_int64_vector_ptr SV, mat_int64_srcptr Mqr)
{
  ASSERT(Mqr->NumRows == Mqr->NumCols);
  ASSERT(Mqr->NumRows == 3);

  int64_vector_t * v = malloc(sizeof(int64_vector_t) * Mqr->NumRows);
  for (unsigned int i = 0; i < Mqr->NumCols; i++) {
    int64_vector_init(v[i], Mqr->NumRows);
    mat_int64_extract_vector(v[i], Mqr, i);
  }
  SV4(SV, v[0], v[1], v[2]);

  for (unsigned int i = 0; i < Mqr->NumCols; i++) {
    int64_vector_clear(v[i]);
  }
  free(v);
}

/*
 * Return the gap between the x coordinate and the closer border of the sieving
 * region defined by the sieving bound H.
 *
 * v: current vector.
 * H: sieving bound.
 */
static int64_t difference_bound_x(int64_vector_srcptr v,
    sieving_bound_srcptr H)
{
  return MIN(ABS(-(int64_t)H->h[0] - v->c[0]), ABS((int64_t)H->h[0] - 1 -
        v->c[0]));
}

/*
 * Return the position of the vector with a x coordinate minimal, according to
 * the classification defined by the stamp array and the value of the
 * classification val_stamp.
 *
 * SV: list of vector.
 * H: sieving bound.
 * stamp: array with length equal to SV, gives the classification of the
 *  vectors.
 * val_stamp: value we want for the stamp of the vector.
 */
unsigned int find_min_x(list_int64_vector_srcptr SV, sieving_bound_srcptr H,
    unsigned char * stamp, unsigned char val_stamp)
{
  int64_t x = -1; 
  unsigned int pos = 0;

  for (unsigned int i = 1; i < SV->length; i++) {
    if (stamp[i] == val_stamp) {
      int64_t tmp = difference_bound_x(SV->v[i], H);
      if (x == -1) {
        x = tmp;
      }
      if (x >= tmp) {
        x = tmp;
        pos = i;
      }
    }
  }
  return pos;
}

/*
 * Add an FK vector (e0 or e1) if v is outside of the sieving region defined by
 * H to have the x coordinate of v in [-H0, H0[.
 *
 * v: current vector.
 * e0: a vector given by the Franke-Kleinjung algorithm.
 * e1: a vector given by the Franke-Kleinjung algorithm.
 * H: sieving bound.
 */
void add_FK_vector(int64_vector_ptr v, int64_vector_srcptr e0,
    int64_vector_srcptr e1, sieving_bound_srcptr H)
{
  ASSERT(v->c[0] < -(int64_t)H->h[0] || v->c[0] >= (int64_t)H->h[0]);

  int64_t dist = difference_bound_x(v, H);
  if (-e0->c[0] > e1->c[0]) {
    unsigned int nb = (unsigned int)ceil((double)dist / -(double)e0->c[0]);
    if (v->c[0] > 0) {
      for (unsigned int i = 0; i < nb; i++) {
        int64_vector_add(v, v, e0);
      }
    } else {
      for (unsigned int i = 0; i < nb; i++) {
        int64_vector_sub(v, v, e0);
      }
    }
  } else {
    unsigned int nb = (unsigned int)ceil((double)dist / (double)e1->c[0]);
    if (v->c[0] > 0) {
      for (unsigned int i = 0; i < nb; i++) {
        int64_vector_sub(v, v, e1);
      }
    } else {
      for (unsigned int i = 0; i < nb; i++) {
        int64_vector_add(v, v, e1);
      }
    }
  }
  
  ASSERT(v->c[0] < (int64_t)H->h[0]);
  ASSERT(v->c[0] >= -(int64_t)H->h[0]);
}

/*
 * Return 0 if v0 is added, 1 if v1 is added, 2 if v0 + v1 is added.
 */
unsigned int enum_pos_with_FK(int64_vector_ptr v, int64_vector_srcptr v_old,
    int64_vector_srcptr v0, int64_vector_srcptr v1, int64_t A, int64_t I)
{
  //v0 = (alpha, beta) and v1 = (gamma, delta)
  ASSERT(I > 0);
  ASSERT(v0->c[0] > -I);
  ASSERT(v0->c[0] <= 0);
  ASSERT(v0->c[1] > 0);
  ASSERT(v1->c[0] < I);
  ASSERT(v1->c[0] >= 0);
  ASSERT(v1->c[1] > 0);
  ASSERT(v1->c[0] - v0->c[0] >= I);
#ifndef NDEBUG
  for (unsigned int j = 3; j < v0->dim; j++) {
    ASSERT(v0->c[j] == 0);
    ASSERT(v1->c[j] == 0);
  }
#endif // NDEBUG
  ASSERT(v_old->c[0] < A + I);
  ASSERT(v_old->c[0] >= A);

  if (v_old->c[0] >= A - v0->c[0]) {
    int64_vector_add(v, v_old, v0);
    return 0;
  } else if (v_old->c[0] < A + I - v1->c[0]) {
    int64_vector_add(v, v_old, v1);
    return 1;
  } else {
    int64_vector_add(v, v_old, v0);
    int64_vector_add(v, v, v1);
    return 2;
  }
}

unsigned int enum_neg_with_FK(int64_vector_ptr v, int64_vector_srcptr v_old,
    int64_vector_srcptr v0, int64_vector_srcptr v1, int64_t A, int64_t I)
{
  //v0 = (alpha, beta) and v1 = (gamma, delta)
  ASSERT(I > 0);
  ASSERT(v0->c[0] > -I);
  ASSERT(v0->c[0] <= 0);
  ASSERT(v0->c[1] > 0);
  ASSERT(v1->c[0] < I);
  ASSERT(v1->c[0] >= 0);
  ASSERT(v1->c[1] > 0);
  ASSERT(v1->c[0] - v0->c[0] >= I);
#ifndef NDEBUG
  for (unsigned int j = 3; j < v0->dim; j++) {
    ASSERT(v0->c[j] == 0);
    ASSERT(v1->c[j] == 0);
  }
#endif // NDEBUG
  ASSERT(v_old->c[0] > -A - I);
  ASSERT(v_old->c[0] <= -A);

  if (v_old->c[0] <= -A + v0->c[0]) {
    int64_vector_sub(v, v_old, v0);
    return 0;
  } else if (v_old->c[0] > -A - I + v1->c[0]) {
    int64_vector_sub(v, v_old, v1);
    return 1;
  } else {
    int64_vector_sub(v, v_old, v0);
    int64_vector_sub(v, v, v1);
    return 2;
  }
}

void coordinate_FK_vector(uint64_t * coord_v0, uint64_t * coord_v1,
    int64_vector_srcptr v0, int64_vector_srcptr v1, sieving_bound_srcptr H,
    uint64_t number_element)
{
  ASSERT(2*H->h[0] > 0);
  ASSERT(v0->c[0] > -2*(int64_t)H->h[0]);
  ASSERT(v0->c[0] <= 0);
  ASSERT(v0->c[1] > 0);
  ASSERT(v1->c[0] < 2*(int64_t)H->h[0]);
  ASSERT(v1->c[0] >= 0);
  ASSERT(v1->c[1] > 0);
  ASSERT(v1->c[0] - v0->c[0] >= 2*(int64_t)H->h[0]);
#ifndef NDEBUG
  for (unsigned int j = 3; j < v0->dim; j++) {
    ASSERT(v0->c[j] == 0);
    ASSERT(v1->c[j] == 0);
  }
#endif // NDEBUG

  int64_vector_t e0;
  int64_vector_init(e0, v0->dim);
  int64_vector_set(e0, v0);
  int64_vector_t e1;
  int64_vector_init(e1, v1->dim);
  int64_vector_set(e1, v1);
 
  e0->c[0] = e0->c[0] + (int64_t)(H->h[0] - 1);
  e1->c[0] = e1->c[0] - (int64_t)(H->h[0]);
  e0->c[1] = e0->c[1] - (int64_t)(H->h[1]);
  e1->c[1] = e1->c[1] - (int64_t)(H->h[1]);

  if (int64_vector_in_sieving_region(e0, H)) {
    * coord_v0 = array_int64_vector_index(e0, H, number_element) - (2 *
        (int64_t) H->h[0] - 1);
  } else {
    * coord_v0 = 0; 
  }
  if (int64_vector_in_sieving_region(e1, H)) {
    * coord_v1 = array_int64_vector_index(e1, H, number_element);
  } else {
    * coord_v1 = 0;
  }

  int64_vector_clear(e0);
  int64_vector_clear(e1);
}

void plane_sieve_next_plane(int64_vector_ptr vs, list_int64_vector_srcptr SV,
    int64_vector_srcptr e0, int64_vector_srcptr e1, sieving_bound_srcptr H)
{
  list_int64_vector_t list_int64_vector_tmp;
  list_int64_vector_init(list_int64_vector_tmp);
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
    list_int64_vector_add_int64_vector(list_int64_vector_tmp, v_tmp);

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
    for (unsigned int i = 1; i < list_int64_vector_tmp->length; i++) {
      if (assert[i] == 0) {
        int64_t tmp= ABS(vs->c[0] - list_int64_vector_tmp->v[i]->c[0]);
        if (min == -1) {
          min = tmp;
        }
        if (min >= tmp) {
          min = tmp;
          pos = i;
        }
      }
    }
    int64_vector_set(vs, list_int64_vector_tmp->v[pos]);
    ASSERT(vs->c[0] < (int64_t)H->h[0]);
    ASSERT(vs->c[0] >= -(int64_t)H->h[0]);
  } else if (found == 1) {
    pos = find_min_x(list_int64_vector_tmp, H, assert, 1);
    int64_vector_set(vs, list_int64_vector_tmp->v[pos]);
    ASSERT(vs->c[0] < (int64_t)H->h[0]);
    ASSERT(vs->c[0] >= -(int64_t)H->h[0]);
  } else {
    ASSERT(found == 2);
    pos = find_min_x(list_int64_vector_tmp, H, assert, 2);
    int64_vector_set(vs, list_int64_vector_tmp->v[pos]);
    add_FK_vector(vs, e0, e1, H);
    ASSERT(vs->c[0] < (int64_t)H->h[0]);
    ASSERT(vs->c[0] >= -(int64_t)H->h[0]);
  }

  free(assert);
  list_int64_vector_clear(list_int64_vector_tmp);
}
#ifdef MAIN
void double_vector_gram_schmidt(list_double_vector_ptr list_new,
    mat_double_ptr m, list_double_vector_srcptr list_old)
{
  for (unsigned int i = 0; i < list_old->length; i++) {
    list_double_vector_add_double_vector(list_new, list_old->v[i]);
  }

  double_vector_t v_tmp;
  double_vector_init(v_tmp, list_new->v[0]->dim);
  for (unsigned int i = 0; i < list_new->length; i++) {
    m->coeff[i + 1][i + 1] = 1;
  }

  for (unsigned int i = 0; i < list_new->length; i++) {
    for (unsigned int j = 0; j < i; j++) {
      m->coeff[i + 1][j + 1] = double_vector_orthogonal_projection(v_tmp,
          list_new->v[j], list_new->v[i]); 
      double_vector_sub(list_new->v[i], list_new->v[i], v_tmp);
    }
  }

  double_vector_clear(v_tmp);
}

static double sum_mi(int64_vector_srcptr x, mat_double_srcptr m, unsigned int i)
{
  double sum = 0;
  for (unsigned int j = i + 1; j < x->dim; j++) {
    sum = sum + (double)x->c[j] * m->coeff[j + 1][i + 1];
  }
  return sum;
}

static double sum_li(double_vector_srcptr l, unsigned int i)
{
  double sum = 0;
  for (unsigned int j = i; j < l->dim; j++) {
    sum = sum + l->c[j];
  }
  return sum;
}

static void construct_v(int64_vector_ptr v, int64_vector_srcptr x,
    list_int64_vector_srcptr list)
{
  ASSERT(x->dim == list->length);
  ASSERT(v->dim == list->v[0]->dim);

  for (unsigned int i = 0; i < v->dim; i++) {
    for (unsigned int j = 0; j < x->dim; j++) {
      v->c[i] = v->c[i] + x->c[j] * list->v[j]->c[i];
    }
  }
}

/*
 * cf http://perso.ens-lyon.fr/guillaume.hanrot/Papers/iwcc.pdf
 */
void enum_lattice(mat_int64_srcptr Mqr, sieving_bound_srcptr H)
{
  ASSERT(Mqr->NumRows == Mqr->NumCols);
  ASSERT(Mqr->NumRows == 3);
  ASSERT(H->t == Mqr->NumRows);

  //Original basis of the lattice.
  list_int64_vector_t b_root;
  list_int64_vector_init(b_root);
  list_int64_vector_extract_mat_int64(b_root, Mqr);

  /* To compute and store the result of Gram-Schmidt. */
  list_double_vector_t list_e;
  list_double_vector_init(list_e);
  list_double_vector_extract_mat_int64(list_e, Mqr);

  //Matrix with the coefficient mu_{i, j}.
  mat_double_t M;
  mat_double_init(M, Mqr->NumRows, Mqr->NumCols);
  mat_double_set_zero(M);

  //Gram Schmidt orthogonalisation.
  list_double_vector_t list;
  list_double_vector_init(list);
  double_vector_gram_schmidt(list, M, list_e);
  
  /* Compute the square of the L2 norm for all the Gram-Schmidt vectors. */
  double_vector_t b;
  double_vector_init(b, list->length);
  for (unsigned int i = 0; i < b->dim; i++) {
    b->c[i] = double_vector_norml2sqr(list->v[i]);
  }

  /* Center of the cuboid. */
  double_vector_t t;
  double_vector_init(t, Mqr->NumRows);
  double_vector_set_zero(t);
  t->c[t->dim - 1] = round((double) H->h[t->dim - 1] / 2);

  //This A is equal to the A^2 in the paper we cite.
  double A = 0;
  for (unsigned int i = 0; i < t->dim - 1; i++) {
    A = A + (double)(H->h[0] * H->h[0]);
  }
  A = A + (t->c[t->dim - 1] * t->c[t->dim - 1]);

  //Verify if the matrix M has the good properties.
#ifdef NDEBUG
  for (unsigned int j = 0; j < list->v[0]->dim; j++) {
    if (j == 0) {
      ASSERT(0 != list->v[0]->c[j]);
    } else {
      ASSERT(0 == list->v[0]->c[j]);
    }
  }
  for(unsigned int i = 1; i < list->length; i++) {
    for (unsigned int j = 0; j < list->v[i]->dim; j++) {
      if (j == i) {
        ASSERT(1 == list->v[i]->c[j]);
      } else {
        ASSERT(0 == list->v[i]->c[j]);
      }
    }
  }
#endif

  //Coefficient of t in the Gram-Schmidt basis.
  double_vector_t ti;
  double_vector_init(ti, t->dim);
  double_vector_set_zero(ti);
  ti->c[ti->dim - 1] = t->c[ti->dim - 1];

  //The vector in the classical basis.
  int64_vector_t x;
  int64_vector_init(x, list->length);
  int64_vector_set_zero(x);
  x->c[x->dim - 1] = (int64_t) ceil(ti->c[x->dim - 1] -
      sqrt(A) / sqrt(b->c[x->dim - 1]));
  unsigned int i = list->length - 1;

  double_vector_t l;
  double_vector_init(l, x->dim);
  double_vector_set_zero(l);

  /*printf("Vector b: \n");*/
  /*list_double_vector_fprintf(stdout, list);*/
  /*printf("A: %f\n", A);*/
  /*printf("Vector t: ");*/
  /*double_vector_fprintf(stdout, t);*/
  /*printf("Mu: \n");*/
  /*mat_double_fprintf(stdout, M);*/
  /*printf("Norm^2 b: ");*/
  /*double_vector_fprintf(stdout, b);*/
  /*printf("Coeff ti: ");*/
  /*double_vector_fprintf(stdout, ti);*/
  /*printf("l: ");*/
  /*double_vector_fprintf(stdout, l);*/
  /*printf("x: ");*/
  /*int64_vector_fprintf(stdout, x);*/

  //TODO: is it true else if? is it true sqrt(A) in the last condition?
  while (i < list->length) {
    double tmp = (double)x->c[i] - ti->c[i] + sum_mi(x, M, i);
    l->c[i] = (tmp * tmp) * b->c[i];
    
    if (i == 0 && sum_li(l, 0) <= A) {
      int64_vector_t v_h;
      int64_vector_init(v_h, x->dim);
      int64_vector_set_zero(v_h);
      //x * orig_base.
      construct_v(v_h, x, b_root);
      if (int64_vector_in_sieving_region(v_h, H)) {
        int64_vector_fprintf(stdout, v_h);
      }
      int64_vector_clear(v_h);
      x->c[0] = x->c[0] + 1;
    } else if (i != 0 && sum_li(l, i) <= A) {
      i = i - 1;
      x->c[i] = (int64_t)ceil(ti->c[i] -
          sum_mi(x, M, i) - sqrt((A - sum_li(l, i + 1)) / b->c[i]));
    } else if (sum_li(l, i)> sqrt(A)) {
      i = i + 1;
      if (i < list->length) {
        x->c[i] = x->c[i] + 1;
      }
    }
  }

  double_vector_clear(l);
  double_vector_clear(b);
  int64_vector_clear(x);
  double_vector_clear(ti);
  mat_double_clear(M);
  double_vector_clear(t);
  list_double_vector_clear(list);
  list_double_vector_clear(list_e);
  list_int64_vector_clear(b_root);
}

// Tqr = [1, 144994, 40026]

int main()
{
  /*double_vector_t v0;*/
  /*double_vector_t v1;*/
  /*double_vector_t v2;*/
  /*mat_double_t M;*/
  /*double_vector_init(v0, 3);*/
  /*double_vector_init(v1, 3);*/
  /*double_vector_init(v2, 3);*/
  /*mat_double_init(M, 3, 3);*/
  /*mat_double_set_zero(M);*/

  /*v0->c[0] = 262109;*/
  /*v0->c[1] = 117115;*/
  /*v0->c[2] = 222083;*/
  /*v1->c[0] = 0;*/
  /*v1->c[1] = 1;*/
  /*v1->c[2] = 0;*/
  /*v2->c[0] = 0;*/
  /*v2->c[1] = 0;*/
  /*v2->c[2] = 1;*/

  /*list_double_vector_t lo;*/
  /*list_double_vector_t ln;*/
  /*list_double_vector_init(lo);*/
  /*list_double_vector_init(ln);*/
  /*list_double_vector_add_double_vector(lo, v0);*/
  /*list_double_vector_add_double_vector(lo, v1);*/
  /*list_double_vector_add_double_vector(lo, v2);*/

  /*double_vector_gram_schmidt(ln, M, lo);*/
  /*list_double_vector_fprintf(stdout, lo);*/
  /*mat_double_fprintf(stdout, M);*/
  /*list_double_vector_fprintf(stdout, ln);*/

  /*double_vector_clear(v0);*/
  /*double_vector_clear(v1);*/
  /*double_vector_clear(v2);*/
  /*mat_double_clear(M);*/
  /*list_double_vector_clear(lo);*/
  /*list_double_vector_clear(ln);*/
 
  sieving_bound_t H;
  sieving_bound_init(H, 3);
  H->h[0] = 128;
  H->h[1] = 128;
  H->h[2] = 128;
  /*H->h[0] = 1;*/
  /*H->h[1] = 1;*/
  /*H->h[2] = 1;*/
  mat_int64_t Mqr;
  mat_int64_init(Mqr, 3, 3);
  Mqr->coeff[1][1] = 262109;
  Mqr->coeff[1][2] = 117115;
  Mqr->coeff[1][3] = 222083;
  Mqr->coeff[2][1] = 0;
  Mqr->coeff[2][2] = 1;
  Mqr->coeff[2][3] = 0;
  Mqr->coeff[3][1] = 0;
  Mqr->coeff[3][2] = 0;
  Mqr->coeff[3][3] = 1;

  /*mat_int64_fprintf(stdout, Mqr);*/

  enum_lattice(Mqr, H);

  sieving_bound_clear(H);
  mat_int64_clear(Mqr);

  return 0;
}

#endif // MAIN
