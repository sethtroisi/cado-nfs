#ifndef UTILS_SIEVE_H
#define UTILS_SIEVE_H 

#include <stdint.h>
#include "cado.h"
#include "utils.h"
#include "int64_vector.h"
#include "mat_int64.h"
#include "array.h"
#include "ideal.h"

typedef struct
{
  int64_vector_t * v;
  unsigned int length;
}  s_list_vector_t;

typedef s_list_vector_t list_vector_t[1];
typedef s_list_vector_t * list_vector_ptr;
typedef const s_list_vector_t * list_vector_srcptr;

void list_vector_init(list_vector_ptr SV);

void list_vector_add_int64_vector(list_vector_ptr SV, int64_vector_srcptr v);

void list_vector_clear(list_vector_ptr SV);

void list_vector_fprintf(FILE * file, list_vector_srcptr SV);

/*
 * Compute the acute Gauss reduction of two vector (v0_root and v1_root) if they
 *  are in the same plane. Return the determinant of U if U is not set as NULL,
 *  0 otherwise.
 *
 * v0: the first output vector.
 * v1: the second output vector.
 * U: unimodular matrix to go from the two first coordinate of v0_root and
 *  v1_root to the two first coordinate of v0 and v1. (v0_root, v1_root) * U =
 *  (v0, v1).
 * v0_root: the first input vector.
 * v1_root: the second input vector.
 */
int gauss_reduction(int64_vector_ptr v0, int64_vector_ptr v1, mat_int64_ptr U,
    int64_vector_srcptr v0_root, int64_vector_srcptr v1_root);

/*
 * Do the same thing as gauss_reduction, but set all of the other two first
 * coordinates of v0 and v1 to zero.
 */
int gauss_reduction_zero(int64_vector_ptr v0, int64_vector_ptr v1,
    mat_int64_ptr U, int64_vector_srcptr v0_root, int64_vector_srcptr v1_root);

/*
 * Use the Franke-Kleinjung to modify the two first coordinates of v0_root and
 *  v1_root, if they respect the conditions. See sieve/las-plattice.h
 *  (03/31/2015). Return 1 if the reduction is done, 0 otherwise.
 *
 * v0: the output first vector.
 * v1: the output second vector.
 * v0_root: the initial first vector.
 * v1_root: the initial second vector
 * I: minimal gap between the two x coordinate of v0 and v1. Usually, I = 2*H0.
 */
int reduce_qlattice(int64_vector_ptr v0, int64_vector_ptr v1,
    int64_vector_srcptr v0_root, int64_vector_srcptr v1_root, int64_t I);

/*
 * Same thing as reduce_qlattice above, but set the other two first coordinates
 *  to zero for v0 and v1.
 */
int reduce_qlattice_zero(int64_vector_ptr v0, int64_vector_ptr v1,
                         int64_vector_srcptr v0_root,
                         int64_vector_srcptr v1_root, int64_t I);

/*
 * The SV4 (for instance, is just SV3 but I keep the name for instance) to find
 *  some short vector in the lattice defined by (v0_root, v1_root, v2)  with a
 *  z coordinate equal to 1. Very dependent on the form of the vector: v0 = (a,
 *  b, 0, …, 0), v1 = (c, d, 0, …, 0) and v2 = (e, 0, 1, 0, …, 0).
 *
 * SV: the list of short vector close to v2.
 * v0_root, v1_root, v2: three vectors that gives a basis of a lattice.
 */
void SV4(list_vector_ptr SV, int64_vector_srcptr v0_root,
         int64_vector_srcptr v1_root, int64_vector_srcptr v2);

/*
 * Same thing as SV4, but the basis vector is in a Mqr matrix (vector in
 * columns).
 */
void SV4_Mqr(list_vector_ptr SV, mat_int64_srcptr Mqr);

/*
 * Enumerate with the Franke-Kleinjung algorithm for x increasing. v->c[0]
 *  is in [A, A + I[.
 *
 * v: v_old + l v0 + m v1, (l, m) in [0, 1]^2. v->c[0] > v_old->c[0].
 * v_old: the starting point vector.
 * v0: a vector given by reduce_qlattice, v0->c[0] < 0.
 * v1: a vector given by reduce_qlattice, v1->c[0] > 0.
 * A: how the interval is centered.
 * I: the width of the interval.
 */
unsigned int enum_pos_with_FK(int64_vector_ptr v, int64_vector_srcptr v_old,
    int64_vector_srcptr v0, int64_vector_srcptr v1, int64_t A, int64_t I);

/*
 * Enumerate with the same technique as in Franke-Kleinjung algorithm but for x
 *  decreasing. v->c[0] is in ]-A - I, -A].
 *
 * v: v_old + l v0 + m v1, (l, m) in [0, 1]^2. v->c[0] < v_old->c[0].
 * v_old: the starting point vector.
 * v0: a vector given by reduce_qlattice, v0->c[0] < 0.
 * v1: a vector given by reduce_qlattice, v1->c[0] > 0.
 * A: how the interval is centered.
 * I: the width of the interval.
 */
unsigned int enum_res_with_FK(int64_vector_ptr v, int64_vector_srcptr v_old,
    int64_vector_srcptr v0, int64_vector_srcptr v1, int64_t A, int64_t I);

/*
 * Plane sieve.
 *
 * array: the array in which we store the norms.
 * r: the ideal we consider.
 * Mqr: the Mqr matrix.
 * H: the sieving bound that defines the sieving region.
 */
void plane_sieve_array(array_ptr array, ideal_1_srcptr r, mat_int64_srcptr Mqr,
    sieving_bound_srcptr H);

#endif // UTILS_SIEVE_H
