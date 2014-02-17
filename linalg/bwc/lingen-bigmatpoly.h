#ifndef BIGMATPOLY_H_
#define BIGMATPOLY_H_

#include "abase.h"
#include "lingen-matpoly.h"

#include "select_mpi.h"

/* This defines an MPI-shared polynomial matrix type */

struct bigmatpoly_s {
    /* XXX the first four fields must be compatible with bigmatpoly_ft_s */
    unsigned int m1;      /* number of block rows, index i */
    unsigned int n1;      /* number of block cols, index j */
    MPI_Comm row;        /* size == n1 */
    MPI_Comm col;        /* size == m1 */

    unsigned int m;     /* total number of rows */
    unsigned int n;     /* total number of cols */
    /* The following three are also in cells */
    unsigned int m0;      /* number of rows per block */
    unsigned int n0;      /* number of cols per block */
    size_t size;
    matpoly * cells;
};
typedef struct bigmatpoly_s bigmatpoly[1];
typedef struct bigmatpoly_s * bigmatpoly_ptr;
typedef const struct bigmatpoly_s * bigmatpoly_srcptr;


#ifdef __cplusplus
extern "C" {
#endif

void bigmatpoly_init(abdst_field ab, bigmatpoly_ptr p, bigmatpoly_srcptr model, unsigned int m, unsigned int n, int len);
void bigmatpoly_finish_init(abdst_field ab, bigmatpoly_ptr p, unsigned int m, unsigned int n, int len);
void bigmatpoly_init_model(bigmatpoly_ptr model, MPI_Comm comm, unsigned int m, unsigned int n);
int bigmatpoly_check_pre_init(bigmatpoly_srcptr p);
matpoly_ptr bigmatpoly_my_cell(bigmatpoly_ptr p);
// void bigmatpoly_realloc(bigmatpoly_ptr p, int newalloc);
void bigmatpoly_zero(abdst_field ab, bigmatpoly_ptr p);
void bigmatpoly_clear(abdst_field ab, bigmatpoly_ptr p);
void bigmatpoly_clear_model(bigmatpoly_ptr p);
void bigmatpoly_set_size(bigmatpoly_ptr p, size_t size);

void bigmatpoly_swap(bigmatpoly_ptr a, bigmatpoly_ptr b);
static inline matpoly * bigmatpoly_part(bigmatpoly_ptr p, unsigned int i, unsigned int j);
static inline matpoly_ptr bigmatpoly_cell(bigmatpoly_ptr p, unsigned int i, unsigned int j);

void bigmatpoly_truncate_loc(abdst_field ab, bigmatpoly_ptr dst, bigmatpoly_ptr src, unsigned int size);
void bigmatpoly_rshift(abdst_field ab, bigmatpoly_ptr dst, bigmatpoly_ptr src, unsigned int k);

void bigmatpoly_mul(abdst_field ab, bigmatpoly c, bigmatpoly a, bigmatpoly b);
void bigmatpoly_mp(abdst_field ab, bigmatpoly c, bigmatpoly a, bigmatpoly b);

void bigmatpoly_gather_mat(abdst_field ab, matpoly dst, bigmatpoly src);
void bigmatpoly_scatter_mat(abdst_field ab, bigmatpoly_ptr dst, matpoly_ptr src);

/* {{{ access interface for bigmatpoly */
static inline matpoly * bigmatpoly_part(bigmatpoly_ptr p, unsigned int i, unsigned int j) {
    return p->cells+i*p->n1+j;
}
static inline matpoly_ptr bigmatpoly_cell(bigmatpoly_ptr p, unsigned int i, unsigned int j) {
    return *bigmatpoly_part(p,i,j);
}
/* }}} */


#ifdef __cplusplus
}
#endif

#endif	/* BIGMATPOLY_H_ */
