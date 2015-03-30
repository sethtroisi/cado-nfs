#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include "utils_lattice.h"
#include "ideal.h"
#include "utils_mpz.h"
#include "sieving_bound.h"

/*
 * To deal with the list of short vectors.
 */
void SV_init(SV_ptr SV)
{
  SV->length = 0;
  SV->v = (int64_vector_t * ) malloc(sizeof(int64_vector_t) * 4);
}

void SV_add_int64_vector(SV_ptr SV, int64_vector_srcptr v)
{
  if ((SV->length % 4) == 0) {
    SV->v = realloc(SV->v, sizeof(int64_vector_t) * (SV->length + 4));
  }
  int64_vector_init(SV->v[SV->length], v->dim);
  int64_vector_set(SV->v[SV->length], v);
  SV->length++;
}

void SV_clear(SV_ptr SV)
{
  for (unsigned int i = 0; i < SV->length; i++) {
    ASSERT(SV->v[i]->dim != 0);
    int64_vector_clear(SV->v[i]);
  }
  free(SV->v);
  SV->length = 0;
}

void SV_fprintf(FILE * file, SV_srcptr SV)
{
  fprintf(file, "[\n");
  for (unsigned int i = 0; i < SV->length - 1; i++) {
    int64_vector_fprintf(file, SV->v[i]);
  }
  if (SV->length != 0) {
    int64_vector_fprintf(file, SV->v[SV->length - 1]);
  }
  fprintf(file, "]\n");
}

/*
 * From PNpoly
 */
int int64_vector_in_SV(int64_vector_srcptr vec, SV_srcptr SV)
{
  ASSERT(SV->length > 2);
  ASSERT(vec->dim > 1);
#ifndef NDEBUG
  for (unsigned int i = 0; i < SV->length; i++) {
    ASSERT(vec->dim == SV->v[i]->dim);
    for (unsigned int j = 2; j < vec->dim; j++) {
      ASSERT(vec->c[j] == SV->v[i]->c[j]);
    }
  }
#endif // NDEBUG

  for (unsigned int i = 0; i < SV->length; i++) {
    if (int64_vector_equal(vec, SV->v[i]) == 1) {
      return 1;
    }  
  }

  int c = 0;
  unsigned i = 0;
  unsigned j = 0;
  for (i = 0, j = SV->length - 1; i < SV->length; j =
      i++) {
    if ( ((SV->v[i]->c[1] > vec->c[1]) != (SV->v[j]->c[1] > vec->c[1])) &&
        (vec->c[0] < (SV->v[j]->c[0] - SV->v[i]->c[0]) *
         (vec->c[1] - SV->v[i]->c[1]) / (SV->v[j]->c[1] - SV->v[i]->c[1]) +
         SV->v[i]->c[0]) ) {
      c = !c;
    }
  }

  return c;
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
void SV4(SV_ptr SV, int64_vector_srcptr v0_root,
    int64_vector_srcptr v1_root, int64_vector_srcptr v2)
{
  ASSERT(v0_root->dim == v1_root->dim);
#ifndef NDEBUG
  for (unsigned int j = 3; j < v0_root->dim; j++) {
    ASSERT(v0_root->c[j] == v1_root->c[j]);
  }
#endif // NDEBUG

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
  SV_add_int64_vector(SV, u);

  //Build a triangle around the projection of v2 in the plane z = 0.
  if (v2_new_base_x - (double)a < 0) {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] - v0->c[i];
    }
    u->c[2] = 0;
    SV_add_int64_vector(SV, u);
  } else {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] + v0->c[i];
    }
    u->c[2] = 0;
    SV_add_int64_vector(SV, u);
  }

  if (v2_new_base_y - (double)b < 0) {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] - v1->c[i];
    }
    u->c[2] = 0;
    SV_add_int64_vector(SV, u);
  } else {
    for (unsigned int i = 0; i < 2; i++) {
      u->c[i] = SV->v[0]->c[i] + v1->c[i];
    }
    u->c[2] = 0;
    SV_add_int64_vector(SV, u);
  }

#ifndef NDEBUG
  int64_vector_t v_tmp;
  int64_vector_init(v_tmp, v2->dim);
  int64_vector_set(v_tmp, v2);
  for (unsigned int i = 2; i < v_tmp->dim; i++) {
    v_tmp->c[i] = 0;
  }
  ASSERT(int64_vector_in_SV(v_tmp, SV) == 1);
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

void SV4_Mqr(SV_ptr SV, mat_int64_srcptr Mqr)
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
static unsigned int find_min_x(SV_srcptr SV, sieving_bound_srcptr H,
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
static void add_FK_vector(int64_vector_ptr v, int64_vector_srcptr e0,
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

static void coordinate_FK_vector(uint64_t * coord_v0, uint64_t * coord_v1,
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

void plane_sieve_array(array_ptr array, ideal_1_srcptr r,
    mat_int64_srcptr Mqr, sieving_bound_srcptr H)
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
  SV_t SV;
  SV_init(SV);
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
    SV_clear(SV);

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

#ifndef NDEBUG
          uint64_t index_tmp = array_int64_vector_index(v, H,
              array->number_element);
          ASSERT(index_tmp == index_v);
#endif // NDEBUG

        }
        array->array[index_v] = array->array[index_v] - r->log;
        /*int64_vector_fprintf(stdout, v);*/
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

#ifndef NDEBUG
          uint64_t index_tmp = array_int64_vector_index(v, H,
              array->number_element);
          ASSERT(index_tmp == index_v);
#endif // NDEBUG

        }
        array->array[index_v] = array->array[index_v] - r->log;
        /*int64_vector_fprintf(stdout, v);*/
      }
      if (ABS(v->c[0]) < ABS(vs_tmp->c[0])) {
        int64_vector_set(vs_tmp, v);
      }
      FK_value = enum_neg_with_FK(v, v, e0, e1, -(int64_t)H->h[0] + 1, 2 *
          (int64_t) H->h[0]);
    }
    int64_vector_set(vs, vs_tmp);

    /*printf("# Starting point: ");*/
    /*int64_vector_fprintf(stdout, vs);*/

    //Jump in the next plane.
    SV_t SV_tmp;
    SV_init(SV_tmp);
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
      SV_add_int64_vector(SV_tmp, v_tmp);

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
      for (unsigned int i = 1; i < SV_tmp->length; i++) {
        if (assert[i] == 0) {
          int64_t tmp= ABS(vs->c[0] - SV_tmp->v[i]->c[0]);
          if (min == -1) {
            min = tmp;
          }
          if (min >= tmp) {
            min = tmp;
            pos = i;
          }
        }
      }
      int64_vector_set(vs, SV_tmp->v[pos]);
      ASSERT(vs->c[0] < (int64_t)H->h[0]);
      ASSERT(vs->c[0] >= -(int64_t)H->h[0]);
    } else if (found == 1) {
      pos = find_min_x(SV_tmp, H, assert, 1);
      int64_vector_set(vs, SV_tmp->v[pos]);
      ASSERT(vs->c[0] < (int64_t)H->h[0]);
      ASSERT(vs->c[0] >= -(int64_t)H->h[0]);
    } else {
      ASSERT(found == 2);
      pos = find_min_x(SV_tmp, H, assert, 2);
      int64_vector_set(vs, SV_tmp->v[pos]);
      add_FK_vector(vs, e0, e1, H);
      ASSERT(vs->c[0] < (int64_t)H->h[0]);
      ASSERT(vs->c[0] >= -(int64_t)H->h[0]);
    }

    free(assert);
    SV_clear(SV_tmp);
  }

  SV_clear(SV);
  int64_vector_clear(v);
  int64_vector_clear(vs);
  int64_vector_clear(vs_tmp);
  int64_vector_clear(e0);
  int64_vector_clear(e1);
}
