#ifndef __IJVEC_H__
#define __IJVEC_H__

#include <stdint.h>

#include "types.h"


/* Vector in the reduced q-lattice as a pair of polynomials (i,j).
 * Note that, by convention, j should be kept monic.
 *****************************************************************************/

static inline
void ijvec_set_i(ijvec_ptr v, ij_srcptr i)
{
  ijvec_set_ij(v, i);
}

static inline
void ijvec_set_j(ijvec_ptr v, ij_srcptr j, unsigned I)
{
  ijvec_set_ij(v, j);
  ijvec_mul_ti(v, v, I);
}

static inline
void ijvec_set_i_j(ijvec_ptr v, ij_srcptr i, ij_srcptr j, unsigned I)
{
  ijvec_t ii, jj;
  ijvec_set_i(ii, i);
  ijvec_set_j(jj, j, I);
  ijvec_add_disjoint(v, ii, jj);
}

static inline
void ijvec_get_i(ij_ptr i, ijvec_srcptr v, unsigned I)
{
  ijvec_t ii;
  ijvec_mod_ti(ii, v, I);
  ij_set_ijvec(i, ii);
}

static inline
void ijvec_get_j(ij_ptr j, ijvec_srcptr v, unsigned I)
{
  ijvec_t jj;
  ijvec_div_ti(jj, v, I);
  ij_set_ijvec(j, jj);
}

static inline
void ijvec_get_i_j(ij_ptr i, ij_ptr j, ijvec_srcptr v, unsigned I)
{
  ijvec_get_i(i, v, I);
  ijvec_get_j(j, v, I);
}

// Corresponding position as an unsigned int.
typedef uint64_t ijpos_t;

// Return a strict higher bound on the position, given the degree bounds
// I and J.
static inline
ijpos_t ijvec_get_max_pos(unsigned I, unsigned J)
{
  ijvec_t t;
  ijvec_set_ti(t, I+J);
  return ijvec_get_ui(t, I, J);
}

// Return the position corresponding to an (i,j)-vector.
static inline
ijpos_t ijvec_get_pos(ijvec_srcptr v, unsigned I, unsigned J)
{
  return ijvec_get_ui(v, I, J);
}

// Return the starting position of the line corresponding to j.
// Once i is known, just add the offset ijvec_get_offset(i, I) to compute the
// full position of the vector (i,j).
static inline
ijpos_t ijvec_get_start_pos(ij_srcptr j, unsigned I, unsigned J)
{
  ijvec_t jj;
  ijvec_set_j(jj, j, I);
  return ijvec_get_ui(jj, I, J);
}

// Return the position offset corresponding to i.
static inline
ijpos_t ijvec_get_offset(ij_srcptr i, unsigned I)
{
  return ij_get_ui(i, I, 0);
}

// Convert a position to an (i,j)-vector.
// Return 1 if successful.
static inline
int ijvec_set_pos(ijvec_ptr v, ijpos_t pos, unsigned I, unsigned J)
{
  return ijvec_set_ui(v, pos, I, J);
}



/* Basis of the p-lattice seen as a GF(p)-vector space of (i,j)-vectors.
 *****************************************************************************/

// Compute the (i,j)-basis of a given p-lattice.
void ijbasis_compute_large(ijvec_t *basis,  unsigned *basis_dim,
                           unsigned I,      unsigned  J,
                           ijvec_t *euclid, unsigned *euclid_dim,
                           unsigned hatI,   unsigned  hatJ,
                           large_fbideal_srcptr gothp, fbprime_srcptr lambda);

void ijbasis_compute_small(ij_t *basis, ij_t *adjustment_basis,
        small_fbideal_srcptr gothp, fbprime_srcptr lambda,
        unsigned I, unsigned J);

#endif  /* __IJVEC_H__ */
