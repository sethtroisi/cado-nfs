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
#include "list_int64_vector.h"
#include "list_double_vector.h"
#include "list_int64_vector_index.h"

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
MAYBE_UNUSED static int determinant2(mat_int64_srcptr matrix)
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

void skew_LLL(mat_int64_ptr MSLLL, mat_int64_srcptr M_root,
    int64_vector_srcptr skewness)
{
  ASSERT(skewness->dim == M_root->NumCols);
  ASSERT(M_root->NumRows == M_root->NumCols);
  ASSERT(MSLLL->NumRows == M_root->NumCols);
  ASSERT(MSLLL->NumRows == M_root->NumCols);

  mat_int64_t I_s;
  mat_int64_init(I_s, M_root->NumRows, M_root->NumCols);
  mat_int64_set_diag(I_s, skewness);

  mat_int64_t M;
  mat_int64_init(M, M_root->NumRows, M_root->NumCols);
  mat_int64_mul_mat_int64(M, I_s, M_root);
  
  mat_int64_t U;
  mat_int64_init(U, M_root->NumRows, M_root->NumCols);
  mat_int64_LLL_unimodular_transpose(U, M);
  mat_int64_mul_mat_int64(MSLLL, M_root, U);

  mat_int64_clear(U);
  mat_int64_clear(I_s);
  mat_int64_clear(M);
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
  ASSERT(int64_vector_in_polytop_list_int64_vector(v_tmp, SV) == 1);
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
 * Return the gap between the x coordinate and the closest border of the sieving
 * region defined by the sieving bound H.
 *
 * v: current vector.
 * H: sieving bound.
 */
/*static uint64_t difference_bound_x(int64_vector_srcptr v,*/
    /*sieving_bound_srcptr H)*/
/*{*/
  /*return (uint64_t) MIN(ABS(-(int64_t)H->h[0] - v->c[0]),*/
      /*ABS((int64_t)H->h[0] - 1 - v->c[0]));*/
/*}*/

static int64_t compute_k(int64_t x, int64_t xfk, int64_t H0)
{
  /*printf("%" PRId64 ", %" PRId64 ", %" PRId64 "\n", x, xfk, H0);*/
  /*printf("%f\n", (double)(-H0 - x) / (double)xfk);*/
  /*printf("%f\n", (double)(H0 - 1 - x) / (double)xfk);*/

  /*ASSERT(ceil((double)(-H0 - x) / (double)xfk) ==*/
      /*floor((double)(H0 - 1 - x) / (double)xfk - 1));*/

  return (int64_t)ceil((double)(-H0 - x) / (double)xfk);
}

void add_FK_vector(int64_vector_ptr v, list_int64_vector_srcptr list,
    int64_vector_srcptr e0, int64_vector_srcptr e1, sieving_bound_srcptr H)
{
  ASSERT(v->dim == 3);

  int64_vector_t v_use;
  int64_vector_init(v_use, e0->dim);

  //Select the Franke-Kleinjung vector which have the larger x coordinate.
  if (-e0->c[0] > e1->c[0]) {
    ASSERT(v->dim == e0->dim);
    //Set v_use = -e0;
    for (unsigned int i = 0; i < e0->dim; i++) {
      v_use->c[i] = -e0->c[i];
    }
  } else {
    int64_vector_set(v_use, e1);
  }

  ASSERT(v_use->c[0] > 0);

  //Compute the distance.
  unsigned int pos = 0;
  int64_t k = compute_k(list->v[0]->c[0], v_use->c[0], (int64_t)H->h[0]);

  ASSERT(list->v[0]->c[0] + k * v_use->c[0] < (int64_t)H->h[0]);
  ASSERT(list->v[0]->c[0] + k * v_use->c[0] >= -(int64_t)H->h[0]);

  uint64_t min_y = (uint64_t)ABS(list->v[0]->c[1] + k * v_use->c[1]);

  for (unsigned int i = 1; i < list->length; i++) {
    int64_t k_tmp = compute_k(list->v[i]->c[0], v_use->c[0], (int64_t)H->h[0]);
   
    ASSERT(list->v[i]->c[0] + k_tmp * v_use->c[0] < (int64_t)H->h[0]);
    ASSERT(list->v[i]->c[0] + k_tmp * v_use->c[0] >= -(int64_t)H->h[0]);

    uint64_t tmp = (uint64_t)ABS(list->v[i]->c[1] + k_tmp * v_use->c[1]);
    if (tmp < min_y) {
      min_y = tmp;
      pos = i;
      k = k_tmp;
    }
  }

  int64_vector_set(v, list->v[pos]);

  for (unsigned int i = 0; i < 2; i++) {
    v->c[i] = v->c[i] + k * v_use->c[i];
  }

  ASSERT(v->c[0] < (int64_t)H->h[0]);
  ASSERT(v->c[0] >= -(int64_t)H->h[0]);

  int64_vector_clear(v_use);
}

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

//TODO: is it a good idea to retun 0 if fail?
//TODO: index is not signed, why?
uint64_t index_vector(int64_vector_srcptr v, sieving_bound_srcptr H,
    uint64_t number_element)
{
  uint64_t index = 0;

  int64_vector_t v_start;
  int64_vector_init(v_start, v->dim);
  int64_vector_set_zero(v_start);
  int64_vector_t v_end;
  int64_vector_init(v_end, v->dim);
  int64_vector_set(v_end, v);
  for (unsigned int i = 0; i < v->dim - 1; i++) {
    if (v->c[i] >= 0) {
      v_start->c[i] = -(int64_t)H->h[i];
    } else {
      v_start->c[i] = (int64_t)(H->h[i] - 1);
    }
    v_end->c[i] = v_start->c[i] + v->c[i];
  }
  if (int64_vector_in_sieving_region(v_end, H)) {
    uint64_t index_end = array_int64_vector_index(v_end, H, number_element);
    uint64_t index_start =
      array_int64_vector_index(v_start, H, number_element);
    if (index_end > index_start) {
      index = index_end - index_start;
    } else {
      index = index_start - index_end;
    }
  } else {
    index = 0;
  }

  ASSERT(index < number_element);
  int64_vector_clear(v_start);
  int64_vector_clear(v_end);

  return index;
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

  ASSERT(* coord_v0 == index_vector(v0, H, number_element));
  ASSERT(* coord_v1 == index_vector(v1, H, number_element));

  int64_vector_clear(e0);
  int64_vector_clear(e1);
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
unsigned int find_min_y(list_int64_vector_srcptr SV, unsigned char * stamp,
    unsigned char val_stamp)
{
  int64_t y = -1; 
  unsigned int pos = 0;

  for (unsigned int i = 0; i < SV->length; i++) {
    if (stamp[i] == val_stamp) {
      //Find the closest y to 0.
      //TODO: warning, with this, we forget the problem in H1 - 1.
      int64_t tmp = ABS(SV->v[i]->c[1]);
      if (y == -1) {
        y = tmp;
      }
      if (y >= tmp) {
        y = tmp;
        pos = i;
      }
    }
  }
  return pos;
}

void plane_sieve_next_plane(int64_vector_ptr vs, list_int64_vector_srcptr SV,
    int64_vector_srcptr e0, int64_vector_srcptr e1, sieving_bound_srcptr H)
{
  // Contain all the possible vectors to go from z=d to z=d+1.
  list_int64_vector_t list;
  list_int64_vector_init(list, 3);
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
    list_int64_vector_add_int64_vector(list, v_tmp);

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

  if (found < 2) {
    unsigned pos = find_min_y(list, assert, found);
    int64_vector_set(vs, list->v[pos]);
  } else {
    ASSERT(found == 2);
    add_FK_vector(vs, list, e0, e1, H);
  }
  
  ASSERT(vs->c[0] < (int64_t)H->h[0]);
  ASSERT(vs->c[0] >= -(int64_t)H->h[0]);
  free(assert);
  
  list_int64_vector_clear(list);
}

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

static int space_sieve_good_vector(int64_vector_srcptr v,
    sieving_bound_srcptr H)
{
  ASSERT(v->dim == H->t);

  for (unsigned int i = 0; i < v->dim - 1; i++) {
    if ((-2 * (int64_t)H->h[i]) > v->c[i] ||
        (2 * (int64_t)H->h[i]) < v->c[i]) {
      return 0;
    }
  }
  if (0 > v->c[v->dim - 1] || (int64_t)H->h[v->dim - 1] <= v->c[v->dim - 1]) {
    return 0;
  }
  return 1; 
}

static int int64_vector_in_list_zero(int64_vector_srcptr v_tmp,
    list_int64_vector_index_srcptr list_zero)
{
  for (unsigned int i = 0; i < list_zero->length; i++) {
    for (unsigned int j = 0; j < v_tmp->dim - 1; j++) {
      if (ABS(v_tmp->c[j]) == ABS(list_zero->v[i]->vec->c[j])) {
        return 1;
      }
    }
  }
  return 0;
}

//TODO: change that.
static unsigned int space_sieve_linear_combination(
    list_int64_vector_index_ptr list, list_int64_vector_index_ptr list_zero,
    mat_int64_srcptr MSLLL, sieving_bound_srcptr H, uint64_t number_element)
{
  unsigned int vector_1 = 0;

  list_int64_vector_t list_tmp;
  list_int64_vector_init(list_tmp, 3);
  list_int64_vector_extract_mat_int64(list_tmp, MSLLL);

  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, MSLLL->NumRows);
  for (int64_t l = -1; l < 2; l++) {
    for (int64_t m = -1; m < 2; m++) {
      for (int64_t n = -1; n < 2; n++) {
        //v_tmp = l * v[0] + n * v[1] + m * v[2]
        int64_vector_mul(v_tmp, list_tmp->v[0], (int64_t)l);
        int64_vector_addmul(v_tmp, v_tmp, list_tmp->v[1], (int64_t)m);
        int64_vector_addmul(v_tmp, v_tmp, list_tmp->v[2], (int64_t)n);
        if (space_sieve_good_vector(v_tmp, H)) {
          if (v_tmp->c[2] == 0) {
            if (v_tmp->c[1] != 0 || v_tmp->c[0] != 0) {
              int64_vector_reduce(v_tmp, v_tmp);
              if (!int64_vector_in_list_zero(v_tmp, list_zero)) {
                list_int64_vector_index_add_int64_vector_index(list_zero,
                    v_tmp, index_vector(v_tmp, H, number_element));
              }
            }
          } else {
            list_int64_vector_index_add_int64_vector_index(list, v_tmp, 0);
            if (v_tmp->c[2] == 1) {
              vector_1 = 1;
            }
          }
        }
      }
    }
  }

  for (unsigned int i = 0; i < list_zero->length; i++) {
    ASSERT(int64_vector_gcd(list_zero->v[i]->vec) == 1);
  }

  int64_vector_clear(v_tmp);
  list_int64_vector_clear(list_tmp);

  return vector_1;
}
int int64_vector_in_sieving_region_dim(int64_vector_srcptr v,
    sieving_bound_srcptr H) {
  ASSERT(v->dim == H->t);

  for (unsigned int i = 0; i < v->dim - 1; i++) {
    if (-(int64_t)H->h[i] > v->c[i] || ((int64_t)H->h[i]) <= v->c[i]) {
      return 0;
    }
  }
  return 1; 
}

void plane_sieve_1_enum_plane_incomplete(int64_vector_ptr v_in,
    int64_vector_srcptr e0, int64_vector_srcptr e1, sieving_bound_srcptr H)
{
  int64_vector_t v;
  int64_vector_init(v, v_in->dim);
  int64_vector_set(v, v_in);
  //Perform Franke-Kleinjung enumeration.
  //x increases.
  while (v->c[1] < (int64_t)H->h[1]) {
    if (v->c[1] >= -(int64_t)H->h[1]) {
      int64_vector_set(v_in, v);
      int64_vector_clear(v);
      return;
    }
    enum_pos_with_FK(v, v, e0, e1, -(int64_t)H->h[0], 2 * (int64_t)
        H->h[0]);
  }

  //x decreases.
  int64_vector_set(v, v_in);
  enum_neg_with_FK(v, v, e0, e1, -(int64_t)H->h[0] + 1, 2 * (int64_t)H->h[0]);
  while (v->c[1] >= -(int64_t)H->h[1]) {
    if (v->c[1] < (int64_t)H->h[1]) {
      int64_vector_set(v_in, v);
      int64_vector_clear(v);
      return;
    }
    enum_neg_with_FK(v, v, e0, e1, -(int64_t)H->h[0] + 1, 2 *
        (int64_t) H->h[0]);
  }
  int64_vector_clear(v);
}

//SPACE_SIEVE_CONTRIBUTION is broken.
void plane_sieve_1_incomplete(int64_vector_ptr s_out, int64_vector_srcptr s,
    MAYBE_UNUSED mat_int64_srcptr Mqr, sieving_bound_srcptr H,
    list_int64_vector_srcptr list_FK, list_int64_vector_srcptr list_SV)
{
  ASSERT(Mqr->NumRows == Mqr->NumCols);
  ASSERT(Mqr->NumRows == 3);

  int64_vector_set(s_out, s);

  plane_sieve_next_plane(s_out, list_SV, list_FK->v[0], list_FK->v[1], H);
  //Enumerate the element of the sieving region.
  for (unsigned int d = (unsigned int) s->c[2] + 1; d < H->h[2]; d++) {
    plane_sieve_1_enum_plane_incomplete
      (s_out, list_FK->v[0], list_FK->v[1], H);
    if (int64_vector_in_sieving_region(s_out, H)) {
      return;
    }

    //Jump in the next plane.
    plane_sieve_next_plane(s_out, list_SV, list_FK->v[0], list_FK->v[1], H);
  }
}

void space_sieve_1_3D(array_ptr array, ideal_1_srcptr r, mat_int64_srcptr Mqr,
    sieving_bound_srcptr H, uint64_t threshold_hit)
{
  int64_vector_t skewness;
  int64_vector_init(skewness, 3);
  skewness->c[0] = 1;
  skewness->c[1] = 1;
  //TODO: bound?
  skewness->c[2] = (int64_t)(4 * H->h[0] * H->h[0] * H->h[0]) / (int64_t)r->ideal->r;

  mat_int64_t MSLLL;
  mat_int64_init(MSLLL, Mqr->NumRows, Mqr->NumCols);
  skew_LLL(MSLLL, Mqr, skewness);
  int64_vector_clear(skewness);

  uint64_t nb_hit = 0;
  uint64_t expected_hit = (uint64_t)(4 * H->h[0] * H->h[1] * H->h[2]) /
    r->ideal->r;

  for (unsigned int col = 1; col <= MSLLL->NumCols; col++) {
    if (MSLLL->coeff[3][col] < 0) {
      for (unsigned int row = 1; row <= MSLLL->NumRows; row++) {
        MSLLL->coeff[row][col] = - MSLLL->coeff[row][col];
      }
    }
  }
  
  //TODO: we must generate all the zero vectors!
  list_int64_vector_index_t list_vec;
  list_int64_vector_index_init(list_vec, 3);
  list_int64_vector_index_t list_vec_zero;
  list_int64_vector_index_init(list_vec_zero, 3);
  unsigned int vector_1 = space_sieve_linear_combination(
      list_vec, list_vec_zero, MSLLL, H, array->number_element);
  mat_int64_clear(MSLLL);

  list_int64_vector_index_sort_last(list_vec);
  int plane_sieve = 0;

  list_int64_vector_t list_FK;
  list_int64_vector_init(list_FK, 3);
  list_int64_vector_t list_SV;
  list_int64_vector_init(list_SV, 3);

#ifdef SPACE_SIEVE_OUT
  list_int64_vector_t list_out;
  list_int64_vector_init(list_out, 3);
#else
  list_int64_vector_t list_s;
  list_int64_vector_init(list_s, 3);
#endif // SPACE_SIEVE_OUT

  int64_vector_t s;
  int64_vector_init(s, Mqr->NumRows);
  int64_vector_set_zero(s);
  nb_hit++;

#ifdef SPACE_SIEVE_CONTRIBUTION
  uint64_t index_s = array_int64_vector_index(s, H, array->number_element);
  array->array[index_s] = array->array[index_s] - r->log;
#endif // SPACE_SIEVE_CONTRIBUTION

#ifdef SPACE_SIEVE_OUT
  unsigned int first_s = list_int64_vector_add_int64_vector(list_out, s);
  ASSERT(first_s == 0);
#else
  list_int64_vector_add_int64_vector(list_s, s);
#endif // SPACE_SIEVE_OUT

  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, s->dim);
  while (s->c[2] < (int64_t)H->h[2]) {

#ifdef SPACE_SIEVE_OUT
    ASSERT(int64_vector_equal(s, list_out->v[first_s]) == 0);
#else
    ASSERT(int64_vector_equal(s, list_s->v[0]) == 0);
#endif // SPACE_SIEVE_OUT

    for (unsigned int i = 0; i < list_vec_zero->length; i++) {
      int64_vector_add(v_tmp, s, list_vec_zero->v[i]->vec);

#ifdef SPACE_SIEVE_CONTRIBUTION
      uint64_t index_tmp = index_s;
#endif // SPACE_SIEVE_CONTRIBUTION

      while (int64_vector_in_sieving_region(v_tmp, H)) {
        nb_hit++;
#ifdef SPACE_SIEVE_OUT
        list_int64_vector_add_int64_vector(list_out, v_tmp);
#else
        list_int64_vector_add_int64_vector(list_s, v_tmp);
#endif // SPACE_SIEVE_OUT
        
#ifdef SPACE_SIEVE_CONTRIBUTION
        if (list_vec_zero->v[i]->vec->c[1] < 0) {
          index_tmp = index_tmp - list_vec_zero->v[i]->index;
        } else {
          index_tmp = index_tmp + list_vec_zero->v[i]->index;
        }
        array->array[index_tmp] = array->array[index_tmp] - r->log;
#endif // SPACE_SIEVE_CONTRIBUTION

        int64_vector_add(v_tmp, v_tmp, list_vec_zero->v[i]->vec);
      }
      int64_vector_sub(v_tmp, s, list_vec_zero->v[i]->vec);
 
#ifdef SPACE_SIEVE_CONTRIBUTION
      index_tmp = index_s;
#endif // SPACE_SIEVE_CONTRIBUTION

      while (int64_vector_in_sieving_region(v_tmp, H)) {
        nb_hit++;
#ifdef SPACE_SIEVE_OUT
        list_int64_vector_add_int64_vector(list_out, v_tmp);
#else
        list_int64_vector_add_int64_vector(list_s, v_tmp);
#endif // SPACE_SIEVE_OUT

#ifdef SPACE_SIEVE_CONTRIBUTION
        if (list_vec_zero->v[i]->vec->c[1] < 0) {
          index_tmp = index_tmp + list_vec_zero->v[i]->index;
        } else {
          index_tmp = index_tmp - list_vec_zero->v[i]->index;
        }
        array->array[index_tmp] = array->array[index_tmp] - r->log;
#endif // SPACE_SIEVE_CONTRIBUTION

        int64_vector_sub(v_tmp, v_tmp, list_vec_zero->v[i]->vec);
      }
    }
 
#ifdef SPACE_SIEVE_CONTRIBUTION
    unsigned int index_vec = 0;
    unsigned int s_change = 0;
#endif // SPACE_SIEVE_CONTRIBUTION

    unsigned int hit = 0;
    int64_vector_t s_tmp;
    int64_vector_init(s_tmp, s->dim);
    for (unsigned int i = 0; i < s->dim; i++) {
      s_tmp->c[i] = (int64_t)H->h[i];
    }
#ifdef SPACE_SIEVE_OUT
    for (unsigned int i = first_s; i < list_out->length; i++) {
      ASSERT(list_out->v[i]->c[2] == s->c[2]);
      for (unsigned int j = 0; j < list_vec->length; j++) {
        int64_vector_add(v_tmp, list_out->v[i], list_vec->v[j]->vec);
#else
    for (unsigned int i = 0; i < list_s->length; i++) {
      /*printf("list%u: ", i);*/
      /*int64_vector_fprintf(stdout, list_s->v[i]);*/
      /*printf("s: ");*/
      /*int64_vector_fprintf(stdout, s);*/
      ASSERT(list_s->v[i]->c[2] == s->c[2]);
      for (unsigned int j = 0; j < list_vec->length; j++) {
        int64_vector_add(v_tmp, list_s->v[i], list_vec->v[j]->vec);
#endif // SPACE_SIEVE_OUT
        if (int64_vector_in_sieving_region_dim(v_tmp, H)) {
          if (v_tmp->c[2] < s_tmp->c[2]) {
            int64_vector_set(s_tmp, v_tmp);

#ifdef SPACE_SIEVE_CONTRIBUTION
            index_vec = j;
#ifdef SPACE_SIEVE_OUT
            if (int64_vector_equal(s, list_out->v[i])) {
#else
            if (int64_vector_equal(s, list_s->v[i])) {
#endif // SPACE_SIEVE_OUT
              s_change = 1;
            }
#endif // SPACE_SIEVE_CONTRIBUTION

          }
          hit = 1;
        }
      }
    }
    if (hit) {
      int64_vector_set(s, s_tmp);

      if (s->c[2] < (int64_t)H->h[2]) {
        nb_hit++;
#ifdef SPACE_SIEVE_OUT
        first_s = list_int64_vector_add_int64_vector(list_out, s);
#else
        list_int64_vector_delete_elements(list_s);
        list_int64_vector_add_int64_vector(list_s, s);
#endif // SPACE_SIEVE_OUT

#ifdef SPACE_SIEVE_CONTRIBUTION
        if (!s_change) {
          if (list_vec->v[index_vec]->index == 0) {
            list_vec->v[index_vec]->index =
              index_vector(list_vec->v[index_vec]->vec, H,
                  array->number_element);
          }
          index_s = index_s + list_vec->v[index_vec]->index;
        } else {
          ASSERT(s_change == 1);

          index_s = array_int64_vector_index(s, H, array->number_element);
        }
        array->array[index_s] = array->array[index_s] - r->log;
#endif // SPACE_SIEVE_CONTRIBUTION
      }
    }
    if (!hit && expected_hit - nb_hit >= threshold_hit) {
      ASSERT(hit == 0);

      /*fprintf(stderr, "# Need to plane sieve.\n");*/
      if (!plane_sieve) {
        ASSERT(plane_sieve == 0);

        int64_vector_t * vec = malloc(sizeof(int64_vector_t) * Mqr->NumRows);
        for (unsigned int i = 0; i < Mqr->NumCols; i++) {
          int64_vector_init(vec[i], Mqr->NumRows);
          mat_int64_extract_vector(vec[i], Mqr, i);
        }
        if (vector_1 && list_vec->length != 0) {
          unsigned int cpt = 0;
          while (cpt < list_vec->length && list_vec->v[cpt]->vec->c[2] < 2) {
            if (list_vec->v[cpt]->vec->c[2] == 1) {
              list_int64_vector_add_int64_vector(list_SV,
                  list_vec->v[cpt]->vec);
            }
            cpt++;
          }
        } else {
          ASSERT(vector_1 == 0);
          SV4(list_SV, vec[0], vec[1], vec[2]);
        }
        //TODO: Warning, reduce q lattice can ouput nothing.
        int boolean = reduce_qlattice(vec[0], vec[1], vec[0], vec[1],
            (uint64_t)(2 * H->h[0]));
        if (boolean == 0) {
          fprintf(stderr,
              "# Plane sieve (called by space sieve does not support this type of Mqr.\n");
          mat_int64_fprintf_comment(stderr, Mqr);

          for (unsigned int i = 0; i < Mqr->NumCols; i++) {
            int64_vector_clear(vec[i]);
          }
          free(vec);
          int64_vector_clear(v_tmp);
          int64_vector_clear(s);
#ifdef SPACE_SIEVE_OUT
          list_int64_vector_clear(list_out);
#else
          list_int64_vector_clear(list_s);
#endif // SPACE_SIEVE_OUT
          list_int64_vector_index_clear(list_vec);
          list_int64_vector_index_clear(list_vec_zero);
          list_int64_vector_clear(list_FK);
          list_int64_vector_clear(list_SV);

          return;
        }
        list_int64_vector_add_int64_vector(list_FK, vec[0]);
        list_int64_vector_add_int64_vector(list_FK, vec[1]);
        if (space_sieve_good_vector(vec[0], H)) {
          if (!int64_vector_in_list_zero(vec[0], list_vec_zero)) {
            list_int64_vector_index_add_int64_vector_index(list_vec_zero,
                vec[0], index_vector(vec[0], H, array->number_element));
          }
        }
        if (space_sieve_good_vector(vec[1], H)) {
          if (!int64_vector_in_list_zero(vec[1], list_vec_zero)) {
            list_int64_vector_index_add_int64_vector_index(list_vec_zero,
                vec[1], index_vector(vec[1], H, array->number_element));
          }
        }
        for (unsigned int i = 0; i < Mqr->NumCols; i++) {
          int64_vector_clear(vec[i]);
        }
        free(vec);
        plane_sieve = 1;
      }
      int64_vector_t s_out;
      int64_vector_init(s_out, s->dim);
      int64_vector_t v_new;
      int64_vector_init(v_new, s->dim);
      plane_sieve_1_incomplete(s_out, s, Mqr, H, list_FK, list_SV);
      if (int64_vector_in_sieving_region(s_out, H)) {
        nb_hit++;
#ifdef SPACE_SIEVE_OUT
        first_s = list_int64_vector_add_int64_vector(list_out, s_out);
#else
        list_int64_vector_delete_elements(list_s);
        list_int64_vector_add_int64_vector(list_s, s_out);
#endif // SPACE_SIEVE_OUT
        int64_vector_sub(v_new, s_out, s);
        uint64_t index_new = index_vector(v_new, H, array->number_element);
        list_int64_vector_index_add_int64_vector_index(list_vec, v_new,
            index_new);

#ifdef SPACE_SIEVE_CONTRIBUTION
        index_s = index_s + index_new;
        array->array[index_s] = array->array[index_s] - r->log;
#endif // SPACE_SIEVE_CONTRIBUTION

        list_int64_vector_index_sort_last(list_vec);
      }
      int64_vector_set(s, s_out);
      int64_vector_clear(s_out);
      int64_vector_clear(v_new);
    }
    int64_vector_clear(s_tmp);
  }

  int64_vector_clear(v_tmp);
  int64_vector_clear(s);
#ifdef SPACE_SIEVE_OUT
  list_int64_vector_clear(list_out);
#else
  list_int64_vector_clear(list_s);
#endif // SPACE_SIEVE_OUT
  list_int64_vector_index_clear(list_vec);
  list_int64_vector_index_clear(list_vec_zero);
  list_int64_vector_clear(list_FK);
  list_int64_vector_clear(list_SV);
}
