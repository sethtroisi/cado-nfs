#ifndef __SUBLAT_H__
#define __SUBLAT_H__

#include "macros.h"
#include "fppol.h"
#include "types.h"
#include "ijvec.h"

typedef struct {
    int nb;             // the number of valid sublattices
    int deg;            // the degree of the modulus
    // The following fields are meaningless if nb = 1, deg = 0.
    int n;              // the index of the current sublattice (in .lat)
    fppol16_t modulus;  // the modulus used for sublatticing
    fppol16_t lat[9][2];        // a description of all sublattices.
} sublat_struct_t;

// Rem: we have hardcoded the constant 9. Changing it requires changing
// the initializers below, but that's the only thing to do.

typedef sublat_struct_t   sublat_t[1];
typedef sublat_struct_t * sublat_ptr;

#ifdef USE_F2
#define FPPOL16_ZERO {0}
#elif USE_F3
#define FPPOL16_ZERO {0,0}
#else
#error "Please define FPPOL16_ZERO for your field"
#endif

static sublat_t no_sublat = {{
    1,
    0,
    0,
    FPPOL16_ZERO,
    {
        {FPPOL16_ZERO, FPPOL16_ZERO},
        {FPPOL16_ZERO, FPPOL16_ZERO},
        {FPPOL16_ZERO, FPPOL16_ZERO},
        {FPPOL16_ZERO, FPPOL16_ZERO},
        {FPPOL16_ZERO, FPPOL16_ZERO},
        {FPPOL16_ZERO, FPPOL16_ZERO},
        {FPPOL16_ZERO, FPPOL16_ZERO},
        {FPPOL16_ZERO, FPPOL16_ZERO},
        {FPPOL16_ZERO, FPPOL16_ZERO}
    }
}};


// The following sublattice repartition is valid only in charac 2
#ifdef USE_F2
static sublat_t nine_sublat = {{
    9,
    2,
    0,
    { 6 },   // t^2+t = t*(t+1)
    {
        {{0},{1}},
        {{1},{0}}, {{1},{1}}, {{1},{2}}, {{1},{3}},
        {{2},{1}}, {{2},{3}},
        {{3},{1}}, {{3},{2}},
    }
}};
#endif

// One-liner that tells whether sublattices are active.
static inline int use_sublat(sublat_ptr sublat)
{
    return (sublat->nb > 1);
}


// Convert an (i,j)-vector inside a sublattice into the
// (hat i, hat j)-vector in the traditional q-lattice coordinates.
// The conversion to (a,b) is then the usual ij2ab() function.
//   hati = i*modulus + i0
//   hatj = j*modulus + j0
static inline
void ij_convert_sublat(ijvec_t W, ijvec_t V, sublat_ptr sublat)
{
    if (use_sublat(sublat)) {
#ifdef USE_F2
        // Implemented only for the case t*(t+1) in char 2.
        ASSERT(sublat == &nine_sublat[0]);
        ij_t tmp0, tmp1;
        ij_mul_ti(tmp0, V->i, 1);
        ij_mul_ti(tmp1, V->i, 2);
        ij_add(tmp0, tmp0, tmp1);
        ij_set_16(tmp1, sublat->lat[sublat->n][0]);
        ij_add(W->i, tmp0, tmp1);
        ij_mul_ti(tmp0, V->j, 1);
        ij_mul_ti(tmp1, V->j, 2);
        ij_add(tmp0, tmp0, tmp1);
        ij_set_16(tmp1, sublat->lat[sublat->n][1]);
        ij_add(W->j, tmp0, tmp1);
#else
        ASSERT(0);
#endif
    } else {
        ijvec_set(W, V);
    }
}




#endif   /* __SUBLAT_H__ */
