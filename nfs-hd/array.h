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
  Initialise an array.

  array: the array.
  number_element: number of elements in the array.
*/
void array_init(array_ptr array, uint64_t number_element);

/*
  Delete an array.

  array: the array.
*/
void array_clear(array_ptr array);

/*
  Print an array.

  array: the array we want to print.
*/
void array_printf(array_srcptr array);

/*
  Write an array in a file.

  filew: path to the file in which the array will be printed.
  array: the array we want to print.
*/
void array_fprintf(FILE * filew, array_srcptr array);

/*
  Give the vector (mpz_vector) associated with the index in the array.

  v: a vector.
  index: index.
  H: the sieving interval.
  number_element: number of element in the sieving region.
*/
void array_index_mpz_vector(mpz_vector_ptr v, uint64_t index,
                            sieving_interval_srcptr H, uint64_t number_element);

/*
  WARNING: unused function.
  Give the polynomial (mpz_poly) associated with the index in the array.

  f: a polynomial.
  index: index.
  H: the sieving interval.
  array: the array, to know the number of element in the sieving region.
*/
void array_index_mpz_poly(mpz_poly_ptr poly, uint64_t index,
                          sieving_interval_srcptr H, uint64_t number_element);

/*
  Give the index associated with a vector (mpz_vector).

  index: the index of the polynomial.
  v: the vector.
  H: the sieving interval.
  number_element: number of element in the sieving region, only used in the
   NDEBUG mode.
*/
void array_mpz_vector_index(uint64_t * index, mpz_vector_srcptr v,
                            sieving_interval_srcptr H, uint64_t number_element);

/*
  Give the index associated with a vector (int64_vector).

  index: the index of the polynomial.
  v: the vector.
  H: the sieving interval.
  number_element: number of element in the sieving region, only used in the
   NDEBUG mode.
*/
void array_int64_vector_index(uint64_t * index, int64_vector_srcptr v,
                              sieving_interval_srcptr H,
                              uint64_t number_element);

/*
  WARNING: unused function.
  Give the index associated with a polynomial (mpz_poly).

  index: the index of the polynomial.
  poly: the polynomial.
  H: the sieving interval.
*/
void array_mpz_poly_index(uint64_t * index, mpz_poly_srcptr poly,
                          sieving_interval_srcptr H, uint64_t number_element);

#endif /* ARRAY_H */
