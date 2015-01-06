#ifndef UINT64_ARRAY_H
#define UINT64_ARRAY_H

#include <stdio.h>
#include <stdint.h>

typedef struct {
  uint64_t length;
  uint64_t * array;
} s_uint64_array_t;

typedef s_uint64_array_t uint64_array_t[1];
typedef s_uint64_array_t * uint64_array_ptr;
typedef const s_uint64_array_t * uint64_array_srcptr;

void uint64_array_init(uint64_array_ptr array, uint64_t length);

void uint64_array_set_coeff(uint64_array_ptr array, uint64_t index,
                            uint64_t coeff);

void uint64_array_realloc(uint64_array_ptr array, uint64_t number);

void uint64_array_clear(uint64_array_ptr array);

void uint64_array_fprintf(FILE * file, uint64_array_srcptr array);

#endif /* UINT64_ARRAY_H */
