#ifndef INT64_VECTOR_H
#define INT64_VECTOR_H

#include <stdio.h>
#include <math.h>
#include "cado.h"
#include "utils.h"
#include "sieving_interval.h"

typedef struct
{
  unsigned int dim;
  int64_t * c;
} int64_vector_struct_t;

typedef int64_vector_struct_t int64_vector_t[1];
typedef int64_vector_struct_t * int64_vector_ptr;
typedef const int64_vector_struct_t * int64_vector_srcptr;

void int64_vector_init(int64_vector_ptr v, unsigned int d);

void int64_vector_clear(int64_vector_ptr v);

void int64_vector_swap(int64_vector_ptr v0, int64_vector_ptr v1);

void int64_vector_set(int64_vector_ptr v, int64_vector_srcptr s);

static inline int64_t int64_vector_getcoeff(int64_vector_srcptr v,
                                            unsigned int i)
{
  ASSERT(i < v->dim);
  return v->c[i];
}

static inline void int64_vector_setcoeff(int64_vector_ptr v,
                                         unsigned int i, int64_t d)
{
  ASSERT(i < v->dim);
  v->c[i] = d;
}

void int64_vector_setcoordinate(int64_vector_ptr v, unsigned int i,
                                int64_t z);

int int64_vector_equal(int64_vector_srcptr a, int64_vector_srcptr b);

void int64_vector_fprintf(FILE * file, int64_vector_srcptr);

void int64_vector_add_one(int64_vector_ptr v, sieving_interval_srcptr H);

unsigned int int64_vector_add_one_i(int64_vector_ptr v, unsigned int i,
                                    sieving_interval_srcptr H);

void int64_vector_to_mpz_vector(mpz_vector_ptr a, int64_vector_srcptr b);

double int64_vector_norml2sqr(int64_vector_srcptr a);

double int64_vector_norml2(int64_vector_srcptr a);

int64_t int64_vector_dot_product(int64_vector_srcptr v0,
                                 int64_vector_srcptr v1);

#endif
