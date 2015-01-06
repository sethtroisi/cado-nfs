#include "uint64_array.h"
#include <stdio.h>
#include <stdlib.h>
#include "macros.h"
#include <inttypes.h>

void uint64_array_init(uint64_array_ptr array, uint64_t length)
{
  array->length = length;
  array->array = malloc(sizeof(uint64_t) * length);
#ifndef NDEBUG
  for (uint64_t i = 0; i < length; i++) {
    array->array[i] = 0;
  }
#endif
}

void uint64_array_set_coeff(uint64_array_ptr array, uint64_t index,
                            uint64_t coeff)
{
  ASSERT(index < array->length);
  array->array[index] = coeff;
}

void uint64_array_realloc(uint64_array_ptr array, uint64_t number)
{
  array->array = realloc(array->array, sizeof(uint64_t) * number);
  array->length = number;
}

void uint64_array_clear(uint64_array_ptr array)
{
  free(array->array);
  array->length = 0;
}

void uint64_array_printf(uint64_array_srcptr array)
{
  printf("[");
  for (uint64_t i = 0; i < array->length - 1; i++) {
    printf("%" PRIu64 ", ", array->array[i]);
  }
  printf("%" PRIu64 "]\n", array->array[array->length - 1]);
}

void uint64_array_fprintf(FILE * file, uint64_array_srcptr array)
{
  fprintf(file, "[");
  for (uint64_t i = 0; i < array->length - 1; i++) {
    fprintf(file, "%" PRIu64 ", ", array->array[i]);
  }
  fprintf(file, "%" PRIu64 "]\n", array->array[array->length - 1]);
}
