#ifndef POLYMAT_H_
#define POLYMAT_H_

#include "abase.h"

/* This is used only for plingen. */

struct polymat_s {
    unsigned int m;
    unsigned int n;
    size_t size;
    size_t alloc;
    abelt * x;
};
typedef struct polymat_s polymat[1];
typedef struct polymat_s * polymat_ptr;
typedef const struct polymat_s * polymat_srcptr;

struct polymat_ur_s {
    unsigned int m;
    unsigned int n;
    size_t size;
    size_t alloc;
    abelt_ur * x;
};
typedef struct polymat_ur_s polymat_ur[1];
typedef struct polymat_ur_s * polymat_ur_ptr;
typedef const struct polymat_ur_s * polymat_ur_srcptr;


#ifdef __cplusplus
extern "C" {
#endif
void polymat_init(polymat_ptr p, unsigned int m, unsigned int n, int len);
void polymat_realloc(polymat_ptr p, int newalloc);
void polymat_zero(polymat_ptr p);
void polymat_clear(polymat_ptr p);
static inline abelt * polymat_part(polymat_ptr p, unsigned int i, unsigned int j, unsigned int k);
static inline abdst_elt polymat_coeff(polymat_ptr p, unsigned int i, unsigned int j, unsigned int k);
void polymat_ur_init(polymat_ur_ptr p, unsigned int m, unsigned int n, int len);
void polymat_ur_realloc(polymat_ur_ptr p, int newalloc);
void polymat_ur_zero(polymat_ur_ptr p);
void polymat_ur_clear(polymat_ur_ptr p);
static inline abelt_ur * polymat_ur_part(polymat_ur_ptr p, unsigned int i, unsigned int j, unsigned int k);
static inline abdst_elt_ur polymat_ur_coeff(polymat_ur_ptr p, unsigned int i, unsigned int j, unsigned int k);
void bwmat_copy_coeffs(abdst_field ab MAYBE_UNUSED, abelt * x0, int stride0, abelt * x1, int stride1, unsigned int n);
void bwmat_move_coeffs(abdst_field ab MAYBE_UNUSED, abelt * x0, int stride0, abelt * x1, int stride1, unsigned int n);
void polymat_mul(abdst_field ab, polymat c, polymat a, polymat b);
void polymat_mp(abdst_field ab, polymat c, polymat a, polymat b);
void polymat_mp_raw(abdst_field ab,
        polymat c,
        unsigned int xc,
        polymat a,
        unsigned int xa, unsigned int na,
        polymat b,
        unsigned int xb, unsigned int nb);

#ifdef __cplusplus
}
#endif

/* {{{ access interface for polymat */
static inline abelt * polymat_part(polymat_ptr p, unsigned int i, unsigned int j, unsigned int k) {
    /* Assume row-major in all circumstances. Old code used to support
     * various orderings, here we don't */
    return p->x+(k*p->m+i)*p->n+j;
}
static inline abdst_elt polymat_coeff(polymat_ptr p, unsigned int i, unsigned int j, unsigned int k) {
    return *polymat_part(p,i,j,k);
}
/* }}} */
/* {{{ access interface for polymat_ur */
static inline abelt_ur * polymat_ur_part(polymat_ur_ptr p, unsigned int i, unsigned int j, unsigned int k) {
    /* Assume row-major in all circumstances. Old code used to support
     * various orderings, here we don't */
    return p->x+(k*p->m+i)*p->n+j;
}
static inline abdst_elt_ur polymat_ur_coeff(polymat_ur_ptr p, unsigned int i, unsigned int j, unsigned int k) {
    return *polymat_ur_part(p,i,j,k);
}
/* }}} */

#endif	/* POLYMAT_H_ */
