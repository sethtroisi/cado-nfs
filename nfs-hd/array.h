#ifndef ARRAY_H
#define ARRAY_H

#include <stdio.h>
#include "cado.h"
#include "utils.h"
#include "int64_vector.h"

typedef struct {
  uint64_t number_element;
  unsigned char * array;
} s_array_t;

typedef s_array_t array_t[1];
typedef s_array_t * array_ptr;
typedef const s_array_t * array_srcptr;

/*
 * Initialise an array.
 *
 * array: the array.
 * number_element: number of elements in the array.
 */
void array_init(array_ptr array, uint64_t number_element);

/*
 * Delete an array.
 *
 * array: the array.
 */
void array_clear(array_ptr array);

/*
 * Get the ith element of array.
 *
 * array: the array.
 * i: index of the element.
 */
static inline unsigned char array_get(array_srcptr array, uint64_t i)
{
  ASSERT(i < array->number_element);
  return array->array[i];
}

/*
 * Set the ith element of array with val.
 *
 * array: the array.
 * i: index of the element.
 * val: value of the element.
 */
static inline void array_set(array_srcptr array, uint64_t i, unsigned char val)
{
  ASSERT(i < array->number_element);
  array->array[i] = val;
}

/*
 * Write an array in a file.
 *
 * filew: path to the file in which the array will be printed.
 * array: the array we want to print.
 */
void array_fprintf(FILE * filew, array_srcptr array);

/*
 * Give the vector (mpz_vector) associated with the index in the array.
 *
 * v: a vector.
 * index: index.
 * H: the sieving bound.
 * number_element: number of element in the sieving region.
 */
void array_index_mpz_vector(mpz_vector_ptr v, uint64_t index,
    sieving_bound_srcptr H, uint64_t number_element);

/*
 * Give the index associated with a vector (mpz_vector).
 * 
 * index: the index of the polynomial.
 * v: the vector.
 * H: the sieving bound.
 * number_element: number of element in the sieving region, only used in the
 *  NDEBUG mode.
 */
void array_mpz_vector_index(uint64_t * index, mpz_vector_srcptr v,
                            sieving_bound_srcptr H, uint64_t number_element);

/*
 * Give the index associated with a vector (int64_vector).
 *
 * index: the index of the polynomial.
 * v: the vector.
 * H: the sieving bound.
 * number_element: number of element in the sieving region, only used in the
 *  NDEBUG mode.
 */
void array_int64_vector_index(uint64_t * index, int64_vector_srcptr v,
    sieving_bound_srcptr H, uint64_t number_element);

#endif /* ARRAY_H */
