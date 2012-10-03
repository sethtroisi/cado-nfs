#ifndef PLINGEN_H_
#define PLINGEN_H_

/* This contains some definitions for lingen mod p.
 *
 * This code is to be understood as scotch tape around an old and quirky
 * implementation.
 */

#include <stdlib.h>

#include "abase.h"

/* We have flat lists of field elements all the way through */
typedef abelt * bw_mnpoly;
typedef abelt * bw_nbpoly;
typedef abelt * bw_nnpoly;
typedef abelt * bw_bbpoly;
typedef abelt * bw_mbpoly;
typedef abelt * bw_mmmat;
typedef abelt * bw_mnmat;
typedef abelt * bw_nbmat;
typedef abelt * bw_bbmat;
typedef abelt * bw_mbmat;

typedef abelt_ur * bw_mbmat_ur;

typedef struct { unsigned int m,n; abfield ab; } dims;

#ifdef __cplusplus
extern "C" {
#endif

static inline abelt * polymat_alloc(unsigned int m, unsigned int n, int len) {
    return (abelt*) malloc(m*n*len*sizeof(abelt));
}
static inline void polymat_zero(abelt * x, unsigned int m, unsigned int n, int len) {
    memset(x, 0, m*n*len*sizeof(abelt));
}
static inline void polymat_free(abelt * x /*, unsigned int m MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, int len MAYBE_UNUSED */) {
    free(x);
}
static inline abelt * polymat_part(abelt * x, unsigned int m, unsigned int n, unsigned int i, unsigned int j, unsigned int k) {
    /* Assume row-major in all circumstances. Old code used to support
     * various orderings, here we don't */
    return x+(k*m+i)*n+j;
}
static inline abdst_elt polymat_coeff(abelt * x, unsigned int m, unsigned int n, unsigned int i, unsigned int j, unsigned int k) {
    return *polymat_part(x,m,n,i,j,k);
}
static inline abelt_ur * polymat_alloc_ur(unsigned int m, unsigned int n, int len) {
    return (abelt_ur*) malloc(m*n*len*sizeof(abelt_ur));
}
static inline void polymat_zero_ur(abelt_ur * x, unsigned int m, unsigned int n, int len) {
    memset(x, 0, m*n*len*sizeof(abelt_ur));
}
static inline void polymat_free_ur(abelt_ur * x /*, unsigned int m MAYBE_UNUSED, unsigned int n MAYBE_UNUSED, int len MAYBE_UNUSED */) {
    free(x);
}
static inline abelt_ur * polymat_part_ur(abelt_ur * x, unsigned int m, unsigned int n, unsigned int i, unsigned int j, unsigned int k) {
    /* Assume row-major in all circumstances. Old code used to support
     * various orderings, here we don't */
    return x+(k*m+i)*n+j;
}
static inline abdst_elt polymat_coeff_ur(abelt_ur * x, unsigned int m, unsigned int n, unsigned int i, unsigned int j, unsigned int k) {
    return *polymat_part_ur(x,m,n,i,j,k);
}

static inline void bwmat_copy_coeffs(abdst_field ab MAYBE_UNUSED, abelt * x0, int stride0, abelt * x1, int stride1, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        abset(ab, *x0, *x1);
        x0+=stride0;
        x1+=stride1;
    }
}
static inline void bwmat_move_coeffs(abdst_field ab MAYBE_UNUSED, abelt * x0, int stride0, abelt * x1, int stride1, unsigned int n)
{
    ASSERT_ALWAYS(stride0 == stride1); /* Otherwise there's probably no point */
    if (x0 < x1) {
        for(unsigned int i = 0 ; i < n ; i++) {
            abset(ab, x0[i * stride0], x1[i * stride1]);
        }
    } else {
        for(unsigned int i = n ; i-- ; ) {
            abset(ab, x0[i * stride0], x1[i * stride1]);
        }
    }
}



#ifdef __cplusplus
}
#endif

#endif	/* PLINGEN_H_ */
