#ifndef INT64_VECTOR_H
#define INT64_VECTOR_H

#include <stdio.h>
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

void int64_vector_swap(int64_vector_ptr v1, int64_vector_ptr v2);

void int64_vector_set(int64_vector_ptr v, int64_vector_srcptr s);

void int64_vector_setcoordinate(int64_vector_ptr v, unsigned int i,
                                int64_t z);

void int64_vector_printf(int64_vector_srcptr v);

void int64_vector_fprintf(FILE * file, int64_vector_srcptr);

void int64_vector_add_one(int64_vector_ptr v, sieving_interval_srcptr H);

unsigned int int64_vector_add_one_i(int64_vector_ptr v, unsigned int i,
                                    sieving_interval_srcptr H);

#endif
