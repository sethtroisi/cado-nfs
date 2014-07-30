#ifndef BIGPOLYMAT_H_
#define BIGPOLYMAT_H_

#include "mpfq_layer.h"
#include "lingen-polymat.h"

#include "select_mpi.h"

/* This defines an MPI-shared polynomial matrix type */

struct bigpolymat_s {
    unsigned int m;     /* total number of rows */
    unsigned int n;     /* total number of cols */
    unsigned int m1;      /* number of block rows, index i */
    unsigned int n1;      /* number of block cols, index j */
    /* The following three are also in cells */
    unsigned int m0;      /* number of rows per block */
    unsigned int n0;      /* number of cols per block */
    size_t size;
    polymat * cells;
    MPI_Comm row;       /* size == n */
    MPI_Comm col;       /* size == m */
};
typedef struct bigpolymat_s bigpolymat[1];
typedef struct bigpolymat_s * bigpolymat_ptr;
typedef const struct bigpolymat_s * bigpolymat_srcptr;


#ifdef __cplusplus
extern "C" {
#endif

void bigpolymat_bcast_polymat_cutoff(struct polymat_cutoff_info * slot, int rank, MPI_Comm comm);

void bigpolymat_init(abdst_field ab, bigpolymat_ptr p, bigpolymat_srcptr model, unsigned int m, unsigned int n, int len);
void bigpolymat_finish_init(abdst_field ab, bigpolymat_ptr p, unsigned int m, unsigned int n, int len);
void bigpolymat_init_model(bigpolymat_ptr model, MPI_Comm comm, unsigned int m, unsigned int n);
int bigpolymat_check_pre_init(bigpolymat_srcptr p);
polymat_ptr bigpolymat_my_cell(bigpolymat_ptr p);
// void bigpolymat_realloc(bigpolymat_ptr p, int newalloc);
void bigpolymat_zero(abdst_field ab, bigpolymat_ptr p);
void bigpolymat_clear(abdst_field ab, bigpolymat_ptr p);
void bigpolymat_clear_model(bigpolymat_ptr p);
void bigpolymat_set_size(bigpolymat_ptr p, size_t size);

void bigpolymat_swap(bigpolymat_ptr a, bigpolymat_ptr b);
static inline polymat * bigpolymat_part(bigpolymat_ptr p, unsigned int i, unsigned int j);
static inline polymat_ptr bigpolymat_cell(bigpolymat_ptr p, unsigned int i, unsigned int j);

void bigpolymat_truncate_loc(abdst_field ab, bigpolymat_ptr dst, bigpolymat_ptr src, unsigned int size);
void bigpolymat_rshift(abdst_field ab, bigpolymat_ptr dst, bigpolymat_ptr src, unsigned int k);

void bigpolymat_mul(abdst_field ab, bigpolymat c, bigpolymat a, bigpolymat b);
void bigpolymat_mp(abdst_field ab, bigpolymat c, bigpolymat a, bigpolymat b);

void bigpolymat_gather_mat(abdst_field ab, polymat dst, bigpolymat src);
void bigpolymat_scatter_mat(abdst_field ab, bigpolymat_ptr dst, polymat_ptr src);

/* {{{ access interface for bigpolymat */
static inline polymat * bigpolymat_part(bigpolymat_ptr p, unsigned int i, unsigned int j) {
    return p->cells+i*p->n1+j;
}
static inline polymat_ptr bigpolymat_cell(bigpolymat_ptr p, unsigned int i, unsigned int j) {
    return *bigpolymat_part(p,i,j);
}
/* }}} */


#ifdef __cplusplus
}
#endif

#endif	/* BIGPOLYMAT_H_ */
