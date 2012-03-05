#ifndef __IJVEC_H__
#define __IJVEC_H__

#include "fppol.h"
#include "types.h"



/* Vector in the reduced q-lattice as a pair of polynomials (i,j).
 *
 * Note that, by convension, j should be kept monic.
 *****************************************************************************/

// Structure definition.
typedef struct {
  ij_t i;
  ij_t j;
} __ijvec_struct;

// Type and pointer shorthands.
typedef       __ijvec_struct  ijvec_t[1];
typedef       __ijvec_struct *ijvec_ptr;
typedef const __ijvec_struct *ijvec_srcptr;

// Vector addition.
static inline
void ijvec_add(ijvec_ptr r, ijvec_srcptr u, ijvec_srcptr v) {
  ij_add(r->i, u->i, v->i);
  ij_add(r->j, u->j, v->j);
}


// Corresponding position as a pair of unsigned ints.
typedef struct {
  unsigned i;
  unsigned j;
} ijpos2_t;

// Conversion to a pair of unsigned ints.
static inline
ijpos2_t ijvec_get_pos2(ijvec_srcptr u, int I, int J)
{ ijpos2_t pos = { .i = ij_get_ui      (i, I),
                   .j = ij_monic_get_ui(j, J) };
  return pos; }

// Conversion from a pair of unsigned ints.
// Return 1 if successful.
static inline
int ijvec_set_pos2(ijvec_ptr r, ijpos2_t pos, int I, int J)
{ return  ij_set_ui      (r->i, pos.i, I) &&
          ij_monic_set_ui(r->j, pos.j, J);


// Corresponding position as an unsigned int.
typedef unsigned ijpos_t;

// Conversion to an unsigned int.
static inline
ijpos_t ijvec_get_pos(ijvec_srcptr u, int I, int J)
{
  ij_t t;
  ij_mul_ti      (t, u->j, I);
  ij_add_disjoint(t, t, u->i);
  return ij_monic_get_ui(t, I+J);
}

// Conversion from an unsigned int.
// Return 1 if successful.
static inline
int ijvec_set_pos(ijvec_ptr r, ijpos_t pos, int I, int J)
{
  ij_t t;
  if (!ij_monic_set_ui(t, pos, I+J)) return 0;
  ij_mod_ti(u->i, t, I);
  ij_div_ti(u->j, t, I);
  return 1;
}

#endif  /* __IJVEC_H__ */
