#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "sieving_interval.h"
#include "array.h"
#include "macros.h"
#include "int64_vector.h"
#include <inttypes.h>

void array_init(array_ptr array, uint64_t number_element)
{
  ASSERT(number_element != 0);

  array->number_element = number_element;
  array->array = malloc(sizeof(unsigned char) * number_element);
}

void array_clear(array_ptr array)
{
  free(array->array);
  array->number_element = 0;
}

void array_printf(array_srcptr array)
{
  for (uint64_t i = 0; i < array->number_element; i++) {
    printf("%u\n", array->array[i]);
  }
}

void array_fprintf(FILE * filew, array_srcptr array)
{
  for (uint64_t i = 0; i < array->number_element; i++) {
    fprintf(filew, "%u\n", array->array[i]);
  }
}

void array_index_mpz_vector(mpz_vector_ptr v, uint64_t index,
                            sieving_interval_srcptr H, uint64_t number_element)
{
  ASSERT(v->dim == H->t);

  unsigned int k = v->dim - 1;
  uint64_t prod = number_element / ((uint64_t)H->h[k] + 1);
  uint64_t res = index / prod;
  mpz_vector_setcoordinate_si(v, k, res);
  index = index - res * prod;
  k = k - 1;
  prod = prod / (2 * H->h[k] + 1);
  for ( ; k > 0; k--) {
    res = index / prod;
    mpz_vector_setcoordinate_si(v, k, res - H->h[k]);
    index = index - res * prod;
    prod = prod / (2 * H->h[k - 1] + 1);
  }
  res = index / prod;
  mpz_vector_setcoordinate_si(v, k, res - H->h[k]);

#ifndef NDEBUG
  mpz_t tmp;
  mpz_init(tmp);
  for (unsigned int i = 0; i < v->dim - 1; i++) {
    mpz_set(tmp, v->c[i]);
    mpz_abs(tmp, tmp);
    ASSERT(mpz_cmp_ui(tmp, H->h[i]) <= 0);
  }
  mpz_set(tmp, v->c[v->dim - 1]);
  ASSERT(mpz_cmp_ui(tmp, 0) >= 0);
  ASSERT(mpz_cmp_ui(tmp, H->h[H->t - 1]) <= 0);
  mpz_clear(tmp);
#endif
}

void array_index_mpz_poly(mpz_poly_ptr poly, uint64_t index,
                          sieving_interval_srcptr H, uint64_t number_element)
{
  ASSERT(poly->deg < (int)H->t);
  unsigned int k = H->t - 1;
  uint64_t prod = number_element / ((H->h[k] + 1));
  uint64_t res = index / prod;
  mpz_poly_setcoeff_si(poly, k, res);
  index = index - res * prod;
  k = k - 1;
  prod = prod / (2 * H->h[k] + 1);
  for ( ; k > 0; k--) {
    res = index / prod;
    mpz_poly_setcoeff_si(poly, k, res - H->h[k]);
    index = index - res * prod;
    prod = prod / (2 * H->h[k - 1] + 1);
  }
  res = index / prod;
  mpz_poly_setcoeff_si(poly, k, res - H->h[k]);

#ifndef NDEBUG
  mpz_t tmp;
  mpz_init(tmp);
  for (int i = 0; i < poly->deg; i++) {
    mpz_set(tmp, poly->coeff[i]);
    mpz_abs(tmp, tmp);
    ASSERT(mpz_cmp_ui(tmp, H->h[i]) <= 0);
  }
  if (poly->deg == (int)(H->t - 1)) {
    mpz_set(tmp, poly->coeff[poly->deg]);
    ASSERT(mpz_cmp_ui(tmp, 0) >= 0);
    ASSERT(mpz_cmp_ui(tmp, H->h[poly->deg]) <= 0);
  } else if (poly->deg == -1) {

  } else {
    mpz_set(tmp, poly->coeff[poly->deg]);
    mpz_abs(tmp, tmp);
    ASSERT(mpz_cmp_ui(tmp, H->h[poly->deg]) <= 0);
  }
  mpz_clear(tmp);
#endif
}

void array_mpz_vector_index(uint64_t * index, mpz_vector_srcptr v,
                            sieving_interval_srcptr H, uint64_t number_element)
{
#ifndef NDEBUG
  ASSERT(v->dim == H->t);
  mpz_t tmp;
  mpz_init(tmp);
  for (unsigned int i = 0; i < v->dim - 1; i++) {
    mpz_set(tmp, v->c[i]);
    mpz_abs(tmp, tmp);
    ASSERT(mpz_cmp_ui(tmp, H->h[i]) <= 0);
  }
  mpz_set(tmp, v->c[v->dim - 1]);
  ASSERT(mpz_cmp_ui(tmp, 0) >= 0);
  ASSERT(mpz_cmp_ui(tmp, H->h[H->t - 1]) <= 0);
  mpz_clear(tmp);
#else
  mpz_t tmp;
#endif

  mpz_init(tmp);
  uint64_t prod = 1;
  unsigned int k = 0;
  mpz_set(tmp, v->c[k]);
  * index = (uint64_t)(mpz_get_si(tmp) + H->h[k]);

  for (k = 1; k < v->dim - 1; k++) {
    prod = prod * (2 * H->h[k - 1] + 1);
    mpz_set(tmp, v->c[k]);
    * index = * index + (uint64_t)(mpz_get_si(tmp) + H->h[k]) * prod;
  }
  prod = prod * (2 * H->h[k - 1] + 1);
  mpz_set(tmp, v->c[k]);
  * index = * index + (uint64_t)(mpz_get_si(tmp)) * prod;
  mpz_clear(tmp);

  ASSERT(* index < number_element);
}

void array_int64_vector_index(uint64_t * index, int64_vector_srcptr v,
                              sieving_interval_srcptr H,
                              uint64_t number_element)
{
#ifndef NDEBUG
  ASSERT(v->dim == H->t);
  int64_t tmp = 0;
  for (unsigned int i = 0; i < v->dim - 1; i++) {
    tmp = v->c[i];
    if (tmp < 0) {
      tmp = -tmp;
    }
    ASSERT(tmp <= (int64_t)H->h[i]);
  }
  tmp = v->c[v->dim - 1];
  ASSERT(tmp >= 0);
  ASSERT(tmp <= (int64_t)H->h[H->t - 1]);
#endif

  uint64_t prod = 1;
  unsigned int k = 0;
  * index = (uint64_t)v->c[k] + H->h[k];

  for (k = 1; k < v->dim - 1; k++) {
    prod = prod * (2 * H->h[k - 1] + 1);
    * index = * index + ((uint64_t)v->c[k] + H->h[k]) * prod;
  }
  prod = prod * (2 * H->h[k - 1] + 1);
  * index = * index + (uint64_t)(v->c[k]) * prod;

  ASSERT(* index < number_element);
}

void array_mpz_poly_index(uint64_t * index, mpz_poly_srcptr poly,
                          sieving_interval_srcptr H, uint64_t number_element)
{
#ifndef NDEBUG
  ASSERT(poly->deg < (int)H->t);
  mpz_t tmp;
  mpz_init(tmp);
  for (int i = 0; i < poly->deg; i++) {
    mpz_set(tmp, poly->coeff[i]);
    mpz_abs(tmp, tmp);
    ASSERT(mpz_cmp_ui(tmp, H->h[i]) <= 0);
  }
  if (poly->deg == (int)H->t - 1) {
    mpz_set(tmp, poly->coeff[poly->deg]);
    ASSERT(mpz_cmp_ui(tmp, 0) >= 0);
    ASSERT(mpz_cmp_ui(tmp, H->h[poly->deg]) <= 0);
  } else if (poly->deg == -1) {

  } else {
    mpz_set(tmp, poly->coeff[poly->deg]);
    mpz_abs(tmp, tmp);
    ASSERT(mpz_cmp_ui(tmp, H->h[poly->deg]) <= 0);
  }
  mpz_clear(tmp);
#else
  mpz_t tmp;
#endif

  if (poly->deg == (int)H->t - 1) {
    mpz_init(tmp);
    uint64_t prod = 1;
    int k = 0;
    mpz_set(tmp, poly->coeff[k]);
    * index = (uint64_t)(mpz_get_si(tmp) + H->h[k]);

    for (k = 1; k < (int)poly->deg; k++) {
      prod = prod * (2 * H->h[k - 1] + 1);
      mpz_set(tmp, poly->coeff[k]);
      * index = * index + (uint64_t)(mpz_get_si(tmp) + H->h[k]) * prod;
    }
    prod = prod * (2 * H->h[k - 1] + 1);
    mpz_set(tmp, poly->coeff[k]);
    * index = * index + (uint64_t)(mpz_get_si(tmp)) * prod;
    mpz_clear(tmp);
  } else if (0 <= poly->deg && poly->deg < (int)H->t - 1) {
    mpz_init(tmp);
    uint64_t prod = 1;
    int k = 0;
    mpz_set(tmp, poly->coeff[k]);
    * index = (uint64_t)(mpz_get_si(tmp) + H->h[k]);

    for (k = 1; k <= poly->deg; k++) {
      prod = prod * (2 * H->h[k - 1] + 1);
      mpz_set(tmp, poly->coeff[k]);
      * index = * index + (uint64_t)(mpz_get_si(tmp) + H->h[k]) * prod;
    }
    for (k = poly->deg + 1; k < (int)H->t - 1; k++) {
      prod = prod * (2 * H->h[k - 1] + 1);
      * index = * index + (uint64_t)(H->h[k]) * prod;
    }

    mpz_clear(tmp);
  } else if (poly->deg == -1) {
    uint64_t prod = 1;
    unsigned int k = 0;
    * index = (uint64_t)(H->h[k]);

    for (k = 1; k < H->t - 1; k++) {
      prod = prod * (2 * H->h[k - 1] + 1);
      * index = * index + (uint64_t)(H->h[k]) * prod;
    }
  }

  ASSERT(* index < number_element);
}
