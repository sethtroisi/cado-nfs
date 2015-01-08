#ifndef FACTOR_BASE_H
#define FACTOR_BASE_H

#include <stdio.h>
#include <stdint.h>
#include "ideal.h"

/*
  Representation of a factor base.
*/
typedef struct {
  ideal_1_t * factor_base_1;
  uint64_t number_element_1;
  ideal_u_t * factor_base_u;
  uint64_t number_element_u;
  ideal_pr_t * factor_base_pr;
  uint64_t number_element_pr;
} s_factor_base_t;

typedef s_factor_base_t factor_base_t[1];
typedef s_factor_base_t * factor_base_ptr;
typedef const s_factor_base_t * factor_base_srcptr;

/*
  Initialise a factor base.

  factor_base: the factor base.
  number_element: number of element (generally, an upper bound).
  t: dimension of the lattice.
*/
void factor_base_init(factor_base_ptr factor_base, uint64_t number_element_1,
                      uint64_t number_element_u, uint64_t number_element_pr);

/*
  Realloc the factor base. The new number of element must be less than the old
   number of element stored in the factor base.

  factor_base: the factor base.
  new_number_element: the number of element in the factor base.
*/
void factor_base_realloc(factor_base_ptr factor_base, uint64_t number_element_1,
                         uint64_t number_element_u, uint64_t number_element_pr);

/*
  Set an ideal in part at an index. Do not forget to define LINESIEVE if you
   want to set an ideal mod r.

  factor_base: the factor base.
  index: index in the factor base.
  r: the r of the ideal (r, h).
  h: the h of the ideal (r, h).
*/
void factor_base_set_ideal_1_part(factor_base_ptr factor_base,
                                  unsigned int index, uint64_t r,
                                  mpz_poly_srcptr h, unsigned int t);

/*
  Set an ideal in part at an index. Do not forget to define LINESIEVE if you
   want to set an ideal mod r.

  factor_base: the factor base.
  index: index in the factor base.
  r: the r of the ideal (r, h).
  h: the h of the ideal (r, h).
*/
void factor_base_set_ideal_u_part(factor_base_ptr factor_base,
                                  unsigned int index, uint64_t r,
                                  mpz_poly_srcptr h, unsigned int t);

/*
  Set an ideal in part at an index. Do not forget to define LINESIEVE if you
   want to set an ideal mod r.

  factor_base: the factor base.
  index: index in the factor base.
  r: the r of the ideal (r, h).
*/
void factor_base_set_ideal_pr(factor_base_ptr factor_base,
                              unsigned int index, uint64_t r,
                              unsigned int t);

/*
  Delete the factor base bound.

  factor_base: the factor base.
*/
void factor_base_clear(factor_base_ptr factor_base, unsigned int t);

/*
  Write the factor base bound.

  file: the file in which we want to write.
  factor_base: the factor base we want to write.
 */
void factor_base_fprintf(FILE * file, factor_base_srcptr factor_base, unsigned
                         int t);

#endif  /* FACTOR_BASE_H */
