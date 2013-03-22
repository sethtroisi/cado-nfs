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

/* {{{ cut-off structures */
/* This structure is used to decide which algorithm to use for a given
 * input length. This essentially is a function from N*N to a finite set of
 * choices (so far, {0,1} only). The value returned for a balanced input
 * length x*x is:
 * if x >= cut: 1
 * if x < cut:
 *      if table == NULL:
 *              return 0
 *      find last (s,a) in table such that s <= x
 *      return a
 * For unbalanced input length x*y, we compare MIN(x,y) with subdivide.
 * At the moment there is no additional table designed to improve the
 * choice.
 *
 * Cutoffs are used only by polymat.c and bigpolymat.c -- however the
 * tuning program also uses it, so we expose it here too.
 */
/* Such a cutoff structure is used for finding which algorithm to used
 * for:
 *  - a balanced x*x product. -> basecase(0) or karatsuba(1)
 *  - an unbalanced x*y product. -> basecase(0) or split into balanced
 *  sizes(1)
 */

struct polymat_cutoff_info {
    unsigned int cut;
    unsigned int subdivide;
    unsigned int (*table)[2];
    unsigned int table_size;
};
#ifdef __cplusplus
extern "C" {
#endif
void polymat_cutoff_info_init(struct polymat_cutoff_info * c);
void polymat_cutoff_info_clear(struct polymat_cutoff_info * c);
void polymat_set_mul_kara_cutoff(const struct polymat_cutoff_info * new_cutoff, struct polymat_cutoff_info * old_cutoff);
void polymat_set_mp_kara_cutoff(const struct polymat_cutoff_info * new_cutoff, struct polymat_cutoff_info * old_cutoff);
#ifdef __cplusplus
}
#endif
/* }}} */

#ifdef __cplusplus
extern "C" {
#endif
void polymat_init(polymat_ptr p, unsigned int m, unsigned int n, int len);
int polymat_check_pre_init(polymat_srcptr p);
void polymat_realloc(polymat_ptr p, int newalloc);
void polymat_zero(polymat_ptr p);
void polymat_clear(polymat_ptr p);
void polymat_fill_random(abdst_field ab MAYBE_UNUSED, polymat_ptr a, unsigned int size, gmp_randstate_t rstate);
void polymat_swap(polymat_ptr a, polymat_ptr b);
static inline abelt * polymat_part(polymat_ptr p, unsigned int i, unsigned int j, unsigned int k);
static inline abdst_elt polymat_coeff(polymat_ptr p, unsigned int i, unsigned int j, unsigned int k);
void polymat_ur_init(polymat_ur_ptr p, unsigned int m, unsigned int n, int len);
void polymat_ur_realloc(polymat_ur_ptr p, int newalloc);
void polymat_ur_zero(polymat_ur_ptr p);
void polymat_ur_clear(polymat_ur_ptr p);
void polymat_ur_swap(polymat_ur_ptr a, polymat_ur_ptr b);
static inline abelt_ur * polymat_ur_part(polymat_ur_ptr p, unsigned int i, unsigned int j, unsigned int k);
static inline abdst_elt_ur polymat_ur_coeff(polymat_ur_ptr p, unsigned int i, unsigned int j, unsigned int k);


void bwmat_copy_coeffs(abdst_field ab MAYBE_UNUSED, abelt * x0, int stride0, abelt * x1, int stride1, unsigned int n);
void bwmat_move_coeffs(abdst_field ab MAYBE_UNUSED, abelt * x0, int stride0, abelt * x1, int stride1, unsigned int n);

void polymat_addmat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb);

void polymat_submat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb);

void polymat_reducemat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat_ur a, unsigned int ka);


void polymat_mul(abdst_field ab, polymat c, polymat a, polymat b);
void polymat_mp(abdst_field ab, polymat c, polymat a, polymat b);
void polymat_mp_raw(abdst_field ab,
        polymat c,
        unsigned int xc,
        polymat a,
        unsigned int xa, unsigned int na,
        polymat b,
        unsigned int xb, unsigned int nb,
        int transpose, int add);

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
