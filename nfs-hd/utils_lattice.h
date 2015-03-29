#ifndef UTILS_SIEVE_H
#define UTILS_SIEVE_H 

#include <stdint.h>
#include "cado.h"
#include "utils.h"
#include "int64_vector.h"
#include "mat_int64.h"

typedef struct
{
  int64_vector_t * v;
  unsigned int length;
}  s_SV_t;

typedef s_SV_t SV_t[1];
typedef s_SV_t * SV_ptr;
typedef const s_SV_t * SV_srcptr;

void SV_init(SV_ptr SV);

void SV_add_int64_vector(SV_ptr SV, int64_vector_srcptr v);

void SV_clear(SV_ptr SV);

void SV_fprintf(FILE * file, SV_srcptr SV);

/*
 * Return the determinant of U if U is not set as NULL, 0 otherwise.
 */
int gauss_reduction(int64_vector_ptr v0, int64_vector_ptr v1, mat_int64_ptr U,
    int64_vector_srcptr v0_root, int64_vector_srcptr v1_root);

int gauss_reduction_zero(int64_vector_ptr v0, int64_vector_ptr v1,
    mat_int64_ptr U, int64_vector_srcptr v0_root, int64_vector_srcptr v1_root);

/*
  Use the Franke-Kleinjung to modify the two first coordinates of v0_root and v1_root, if they respect the conditions.
  v0: the output first vector
  v1: the output second vector
  v0_root: the initial first vector
  v1_root: the initial second vector H0: H0 such as defined for the sieving region. In the
   Franke-Kleinjung article, I = 2*H0.
*/
int reduce_qlattice(int64_vector_ptr v0, int64_vector_ptr v1,
    int64_vector_srcptr v0_root, int64_vector_srcptr v1_root, int64_t I);

int reduce_qlattice_zero(int64_vector_ptr v0, int64_vector_ptr v1,
                         int64_vector_srcptr v0_root,
                         int64_vector_srcptr v1_root, int64_t I);

void enum_pos_with_FK(int64_vector_ptr v, int64_vector_srcptr v_old,
                      int64_vector_srcptr v0, int64_vector_srcptr v1, int64_t A,
                      int64_t I);


void enum_neg_with_FK(int64_vector_ptr v, int64_vector_srcptr v_old,
                      int64_vector_srcptr v0, int64_vector_srcptr v1, int64_t A,
                      int64_t I);

void SV4(SV_ptr SV, int64_vector_srcptr v0_root,
         int64_vector_srcptr v1_root, int64_vector_srcptr v2);

void SV4_Mqr(SV_ptr SV, mat_int64_srcptr Mqr);

void plane_sieve(mat_int64_srcptr Mqr, sieving_bound_srcptr H);

#endif // UTILS_SIEVE_H
