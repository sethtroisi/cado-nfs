#ifndef __IJVEC_H__
#define __IJVEC_H__

#include <stdint.h>

#include "types.h"


/* Vector in the reduced q-lattice as a pair of polynomials (i,j).
 * Note that, by convention, j should be kept monic.
 *****************************************************************************/

// Vector assignement, addition, mul_by_t^i, ...
static inline
void ijvec_add(ijvec_ptr r, ijvec_srcptr u, ijvec_srcptr v)
{
  ij_add(r->i, u->i, v->i);
  ij_add(r->j, u->j, v->j);
}

static inline
void ijvec_set_zero(ijvec_ptr r) {
  ij_set_zero(r->i);
  ij_set_zero(r->j);
}

static inline
void ijvec_set(ijvec_ptr r, ijvec_srcptr u) {
  ij_set(r->i, u->i);
  ij_set(r->j, u->j);
}

static inline
void ijvec_mul_ti(ijvec_ptr r, ijvec_srcptr u, unsigned i) {
  ij_mul_ti(r->i, u->i, i);
  ij_mul_ti(r->j, u->j, i);
}

// Corresponding position as an unsigned int.
typedef uint64_t ijpos_t;

// Return a strict higher bound on the position, given the degree bounds
// I and J.
static inline
ijpos_t ijvec_get_max_pos(unsigned I, unsigned J)
{
  fppol64_t t;
  fppol64_set_ti(t, I+J);
  return fppol64_get_ui(t, I, J);
}

// Return the position corresponding to an (i,j)-vector.
static inline
ijpos_t ijvec_get_pos(ijvec_srcptr v, unsigned I, unsigned J)
{
  fppol64_t ii, jj;
  fppol64_set_ij       (ii, v->i);
  fppol64_set_ij       (jj, v->j);
  fppol64_mul_ti       (jj, jj, I);
  fppol64_add_disjoint (jj, jj, ii);
  return fppol64_get_ui(jj, I, J);
}

// Return the starting position of the line corresponding to j.
// Once i is known, just add the offset ijvec_get_offset(i, I) to compute the
// full position of the vector (i,j).
static inline
ijpos_t ijvec_get_start_pos(ij_srcptr j, unsigned I, unsigned J)
{
  fppol64_t jj;
  fppol64_set_ij       (jj, j);
  fppol64_mul_ti       (jj, jj, I);
  return fppol64_get_ui(jj, I, J);
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
  fppol64_t ii, jj;
  if (!fppol64_set_ui(jj, pos, I, J))
    return 0;
  fppol64_mod_ti(ii, jj,  I);
  fppol64_div_ti(jj, jj,  I);
  ij_set_64(v->i, ii);
  ij_set_64(v->j, jj);
  return 1;
}



/* Basis of the p-lattice seen as a GF(p)-vector space of (i,j)-vectors.
 *****************************************************************************/

// Allocate memory for an (i,j)-basis of degrees bounded by I and J.
void ijbasis_init(ijbasis_ptr basis, unsigned I, unsigned J);

// Compute the (i,j)-basis of a given p-lattice.
void ijbasis_compute_large(ijbasis_ptr euclid, ijbasis_ptr basis,
        large_fbideal_srcptr gothp, fbprime_srcptr lambda);

void ijbasis_compute_small(ij_t *basis, ij_t *adjustment_basis,
        small_fbideal_srcptr gothp, fbprime_srcptr lambda,
        unsigned I, unsigned J);


// Clean up memory.
void ijbasis_clear(ijbasis_ptr basis);

#endif  /* __IJVEC_H__ */
