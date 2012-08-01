#ifndef __SUBLAT_H__
#define __SUBLAT_H__

#include "macros.h"
#include "fppol.h"
#include "types.h"
#include "ijvec.h"

#define MAX_SUBLAT 9

typedef struct {
    int nb;             // the number of valid sublattices
    int deg;            // the degree of the modulus
    // The following fields are meaningless if nb = 1, deg = 0.
    int n;              // the index of the current sublattice (in .lat)
    fppol16_t modulus;  // the modulus used for sublatticing
    fppol16_t lat[MAX_SUBLAT][2]; // a description of all sublattices.
} __sublat_struct;

typedef       __sublat_struct  sublat_t[1];
typedef       __sublat_struct *sublat_ptr;
typedef const __sublat_struct *sublat_srcptr;

static MAYBE_UNUSED sublat_t no_sublat = {{
    1,
    0,
    0,
    {0},
    {{{0},{0}}}
}};


// The following sublattice repartition is valid only in charac 2
#ifdef USE_F2
static MAYBE_UNUSED sublat_t nine_sublat = {{
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
static inline int use_sublat(sublat_srcptr sublat)
{
    return (sublat->nb > 1);
}


// Convert an (i,j)-vector inside a sublattice into the
// (hat i, hat j)-vector in the traditional q-lattice coordinates.
// The conversion to (a,b) is then the usual ij2ab() function.
//   hati = i*modulus + i0
//   hatj = j*modulus + j0

static inline
void ij_convert_sublat(ij_t hati, ij_t hatj, ij_t i, ij_t j,
        sublat_ptr sublat)
{
    if (use_sublat(sublat)) {
#ifdef USE_F2
        // Implemented only for the case t*(t+1) in char 2.
        ASSERT(sublat->modulus[0] == 6);
        ij_t tmp0, tmp1;
        ij_mul_ti(tmp0, i, 1);
        ij_mul_ti(tmp1, i, 2);
        ij_add(tmp0, tmp0, tmp1);
        ij_set_16(tmp1, sublat->lat[sublat->n][0]);
        ij_add(hati, tmp0, tmp1);
        ij_mul_ti(tmp0, j, 1);
        ij_mul_ti(tmp1, j, 2);
        ij_add(tmp0, tmp0, tmp1);
        ij_set_16(tmp1, sublat->lat[sublat->n][1]);
        ij_add(hatj, tmp0, tmp1);
#else
        ASSERT(0);
#endif
    } else {
        ij_set(hati, i);
        ij_set(hatj, j);
    }
}

#ifdef USE_F2
// Functions to compute x mod (t^2+t).
static inline
void sublat_mod(fppol8_t res, unsigned int x) {
    uint8_t xt = x & 1U;
    uint8_t xtt = __builtin_parity(x) ^ xt;
    res[0] = xt | (xtt<<1);
}

static inline
void sublat_mod_16(fppol8_t res, fppol16_t x) {
    sublat_mod(res, x[0]);
}

static inline
void sublat_mod_32(fppol8_t res, fppol32_t x) {
    sublat_mod(res, x[0]);
}

static inline
void sublat_mod_ij(fppol8_t res, ij_t x) {
    sublat_mod(res, x[0]);
}

// Functions to compute x div (t^2+t).
static inline
void sublat_div_16(fppol16_t res, fppol16_t x) {
    uint16_t xx = x[0];
    xx >>= 1; // divide by t
    // division by (t+1), using log-time algo.
    uint16_t xx2 = xx ^ (xx>>1);
    uint16_t xx4 = xx2 ^ (xx2>>2);
    uint16_t xx8 = xx4 ^ (xx4>>4);
    uint16_t xx16 = xx8 ^ (xx8>>8);
    res[0] = xx16 >> 1;
}

static inline
void sublat_div_32(fppol32_t res, fppol32_t x) {
    uint32_t xx = x[0];
    xx >>= 1; // divide by t
    // division by (t+1), using log-time algo.
    uint32_t xx2 = xx ^ (xx>>1);
    uint32_t xx4 = xx2 ^ (xx2>>2);
    uint32_t xx8 = xx4 ^ (xx4>>4);
    uint32_t xx16 = xx8 ^ (xx8>>8);
    uint32_t xx32 = xx16 ^ (xx16>>16);
    res[0] = xx32 >> 1;
}

static inline
void sublat_div_ij(ij_t res, ij_t x) {
#if __ij_SIZE == 32
    sublat_div_32(res, x);
#elif __ij_SIZE == 16
    sublat_div_16(res, x);
#else
    ASSERT_ALWAYS(0);
#endif
}

// Arithmetic modulo (t^2 + t), assuming reduced inputs.
static inline
void sublat_add(fppol8_t z, fppol8_t x, fppol8_t y) {
    z[0] = x[0] ^ y[0];
}

static inline
void sublat_sub(fppol8_t z, fppol8_t x, fppol8_t y) {
    sublat_add(z, x, y);
}

static inline
void sublat_mul(fppol8_t z, fppol8_t x, fppol8_t y) {
    uint8_t yy0, yy1;
    yy0 = y[0]; yy1 = yy0 << 1;
    uint8_t xx = x[0], zz;
    zz = (xx & 1u) ? yy0 : 0;
    yy1 = (xx & 2u) ? yy1 : 0;
    zz ^= yy1;  // zz contains now the unreduced result
    uint8_t rr = (zz & 4u); // coefficient in t^2
    z[0] = zz ^ rr ^ (rr>>1);
}

static inline
void sublat_addmul(fppol8_t Z, fppol8_t z, fppol8_t x, fppol8_t y) {
    uint8_t yy0, yy1;
    yy0 = y[0]; yy1 = yy0 << 1;
    uint8_t xx = x[0], zz;
    zz = (xx & 1u) ? yy0 : 0;
    yy1 = (xx & 2u) ? yy1 : 0;
    zz ^= yy1;  // zz contains now the unreduced result
    uint8_t rr = (zz & 4u); // coefficient in t^2
    Z[0] = z[0] ^ zz ^ rr ^ (rr>>1);
}

// Multiply an ij_t by an element of degree at most 1.
static inline
void ij_mul_sublat(ij_t z, ij_t x, fppol8_t y)
{
    z[0] = (y[0] & 1u) ? x[0] : 0;
    ij_t xx;
    xx[0] = x[0] << 1;
    z[0] ^= (y[0] & 2u) ? xx[0] : 0;
}

#endif


#endif   /* __SUBLAT_H__ */
