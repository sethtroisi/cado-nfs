#include "cado.h"
#include <stdlib.h>
#include "macros.h"
#include "utils.h"
#include "lingen-matpoly.h"
#include "lingen-matpoly-ft.h"
#include "logline.h"

/* timings made on cochon, rev 6877b97 (buggy; fixed in 76dde6c) */
#define MP_FTI_DEPTH_ADJ_24_36_36 { { 1, 6 }, { 2, 4 }, { 3, 3 }, { 4, 2 }, { 10, 1 }, { 11, 2 }, { 13, 1 }, { 22, 2 }, { 28, 1 }, { 32, 2 }, { 33, 1 }, { 38, 0 }, { 39, 1 }, { 54, 0 }, { 55, 1 }, { 64, 0 }, { 65, 1 }, { 66, 0 }, { 103, 1 }, { 104, 0 }, { 107, 1 }, { 114, 0 }, { 115, 1 }, { 129, 0 }, }

#define MUL_FTI_DEPTH_ADJ_36_36_36 { { 1, 6 }, { 2, 3 }, { 3, 2 }, { 6, 1 }, { 7, 2 }, { 14, 1 }, { 23, 0 }, { 26, 1 }, { 44, 0 }, { 46, 1 }, { 54, 0 }, { 61, 1 }, { 62, 0 }, }


void matpoly_ft_init(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, unsigned int m, unsigned int n, struct fft_transform_info * fti)
{
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    t->m = m;
    t->n = n;
    t->data = malloc(m * n * fft_alloc_sizes[0]);
    memset(t->data, 0, m * n * fft_alloc_sizes[0]);
    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            fft_transform_prepare(tij, fti);
        }
    }
}

void matpoly_ft_clear(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, struct fft_transform_info * fti MAYBE_UNUSED)
{
    free(t->data);
    memset(t, 0, sizeof(*t));
}

void matpoly_ft_dft(abdst_field ab, matpoly_ft_ptr t, matpoly_ptr a, struct fft_transform_info * fti)
{
    /* Recall that this is for polynomials ! */
    ASSERT_ALWAYS(t->m == a->m);
    ASSERT_ALWAYS(t->n == a->n);

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    void * tt = malloc(fft_alloc_sizes[1]);
    memset(tt, 0, fft_alloc_sizes[1]);

    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            abvec aij = matpoly_part(ab, a, i, j, 0);
            /* ok, casting like this is a crude hack ! */
            fft_do_dft_fppol(tij, (mp_limb_t *) aij, a->size, tt, fti, ab->p);
            logline_printf(2, ".");
        }
        logline_printf(2, " | ");
    }
    free(tt);
}

void matpoly_ft_zero(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, struct fft_transform_info * fti)
{
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            fft_zero(tij, fti);
        }
    }
}

void matpoly_ft_export(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, struct fft_transform_info * fti)
{
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            fft_transform_export(tij, fti);
        }
    }
}

void matpoly_ft_import(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, struct fft_transform_info * fti)
{
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            fft_transform_import(tij, fti);
        }
    }
}

void matpoly_ft_add(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti)
{
    /* Recall that this is for polynomials ! */
    ASSERT_ALWAYS(t0->m == t1->m);
    ASSERT_ALWAYS(t0->n == t1->n);
    ASSERT_ALWAYS(t0->m == u->m);
    ASSERT_ALWAYS(t0->n == u->n);

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);

    for(unsigned int i = 0 ; i < t0->m ; i++) {
        for(unsigned int j = 0 ; j < t0->n ; j++) {
            size_t offset = (i*t0->n+j) * fft_alloc_sizes[0];
            void * t0ij = pointer_arith(t0->data, offset);
            void * t1ij = pointer_arith(t1->data, offset);
            void * uij  = pointer_arith(u->data, offset);
            fft_add(uij, t0ij, t1ij, fti);
        }
    }
}

void matpoly_ft_addmul(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti)
{
    /* Recall that this is for polynomials ! */
    ASSERT_ALWAYS(t0->n == t1->m);
    ASSERT_ALWAYS(t0->m == u->m);
    ASSERT_ALWAYS(t1->n == u->n);

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    size_t tsize = fft_alloc_sizes[0];
    void * qt = malloc(fft_alloc_sizes[1]);
    void * tt = malloc(fft_alloc_sizes[2]);
    memset(qt, 0, fft_alloc_sizes[1]);
    memset(tt, 0, fft_alloc_sizes[2]);

    for(unsigned int i = 0 ; i < t0->m ; i++) {
        for(unsigned int j = 0 ; j < t1->n ; j++) {
            void * uij  = pointer_arith(u->data, (i*u->n+j) * tsize);
            /* now only for mul
            memset(uij, 0, fft_alloc_sizes[0]);
            fft_transform_prepare(uij, fti);
            */
            for(unsigned int k = 0 ; k < t0->n ; k++) {
                void * t0ik = pointer_arith(t0->data, (i*t0->n+k) * tsize);
                void * t1kj = pointer_arith(t1->data, (k*t1->n+j) * tsize);
                fft_addmul(uij, t0ik, t1kj, tt, qt, fti);
            }
        }
    }
    free(qt);
    free(tt);
}

void matpoly_ft_mul(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti)
{
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    size_t tsize = fft_alloc_sizes[0];
    memset(u->data, 0, u->m * u->n * tsize);
    for(unsigned int i = 0 ; i < u->m ; i++) {
        for(unsigned int j = 0 ; j < u->n ; j++) {
            void * uij  = pointer_arith(u->data, (i*u->n+j) * tsize);
            fft_transform_prepare(uij, fti);
        }
    }
    matpoly_ft_addmul(ab, u, t0, t1, fti);
}

/*
void matpoly_ft_sub(abdst_field ab, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti);
*/

void matpoly_ft_ift(abdst_field ab, matpoly_ptr a, matpoly_ft_ptr t, struct fft_transform_info * fti)
{
    /* the caller must allocate first by himself ! */
    ASSERT_ALWAYS(!matpoly_check_pre_init(a));
    /* Recall that this is for polynomials ! */
    ASSERT_ALWAYS(t->m == a->m);
    ASSERT_ALWAYS(t->n == a->n);

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    void * tt = malloc(fft_alloc_sizes[1]);
    memset(tt, 0, fft_alloc_sizes[1]);

    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            abvec aij = matpoly_part(ab, a, i, j, 0);
            /* ok, casting like this is a crude hack ! */
            fft_do_ift_fppol((mp_limb_t *) aij, a->size, tij, tt, fti, ab->p);
            logline_printf(2, ".");
        }
        logline_printf(2, " | ");
    }
    free(tt);
}

void matpoly_ft_ift_mp(abdst_field ab, matpoly_ptr a, matpoly_ft_ptr t, unsigned int shift, struct fft_transform_info * fti)
{
    /* Recall that this is for polynomials ! */
    ASSERT_ALWAYS(t->m == a->m);
    ASSERT_ALWAYS(t->n == a->n);

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    void * tt = malloc(fft_alloc_sizes[1]);
    memset(tt, 0, fft_alloc_sizes[1]);

    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            abvec aij = matpoly_part(ab, a, i, j, 0);
            /* ok, casting like this is a crude hack ! */
            fft_do_ift_fppol_mp((mp_limb_t *) aij, a->size, tij, tt, fti, ab->p, shift);
            logline_printf(2, ".");
        }
        logline_printf(2, " | ");
    }
    free(tt);
}



void matpoly_mul_caching_adj(abdst_field ab, matpoly c, matpoly a, matpoly b, unsigned int adj)/*{{{*/
{
    matpoly_ft tc, ta, tb;
    mpz_t p;
    mpz_init(p);
    abfield_characteristic(ab, p);
    struct fft_transform_info fti[1];
    fft_get_transform_info_fppol(fti, p, a->size, b->size, a->n);
    if (adj != UINT_MAX) {
        fft_transform_info_adjust_depth(fti, adj);
    }
    matpoly_clear(ab, c);
    matpoly_init(ab, c, a->m, b->n, a->size + b->size - 1);
    matpoly_ft_init(ab, ta, a->m, a->n, fti);
    matpoly_ft_init(ab, tb, b->m, b->n, fti);
    matpoly_ft_init(ab, tc, a->m, b->n, fti);
    matpoly_ft_dft(ab, ta, a, fti);
    matpoly_ft_dft(ab, tb, b, fti);
    matpoly_ft_mul(ab, tc, ta, tb, fti);
    c->size = a->size + b->size - 1;
    ASSERT_ALWAYS(c->size <= c->alloc);
    matpoly_ft_ift(ab, c, tc, fti);
    matpoly_ft_clear(ab, ta, fti);
    matpoly_ft_clear(ab, tb, fti);
    matpoly_ft_clear(ab, tc,  fti);
    mpz_clear(p);
}/*}}}*/
void matpoly_mp_caching_adj(abdst_field ab, matpoly c, matpoly a, matpoly b, unsigned int adj)/*{{{*/
{
    matpoly_ft tc, ta, tb;
    mpz_t p;
    mpz_init(p);
    abfield_characteristic(ab, p);
    struct fft_transform_info fti[1];
    fft_get_transform_info_fppol_mp(fti, p, MIN(a->size, b->size), MAX(a->size, b->size), a->n);
    if (adj != UINT_MAX) {
        fft_transform_info_adjust_depth(fti, adj);
    }
    matpoly_clear(ab, c);
    matpoly_init(ab, c, a->m, b->n, MAX(a->size, b->size) - MIN(a->size, b->size) + 1);
    matpoly_ft_init(ab, ta, a->m, a->n, fti);
    matpoly_ft_init(ab, tb, b->m, b->n, fti);
    matpoly_ft_init(ab, tc, a->m, b->n, fti);
    matpoly_ft_dft(ab, ta, a, fti);
    matpoly_ft_dft(ab, tb, b, fti);
    matpoly_ft_mul(ab, tc, ta, tb, fti);
    c->size = MAX(a->size, b->size) - MIN(a->size, b->size) + 1;
    ASSERT_ALWAYS(c->size <= c->alloc);
    matpoly_ft_ift_mp(ab, c, tc, MIN(a->size, b->size) - 1, fti);
    matpoly_ft_clear(ab, ta, fti);
    matpoly_ft_clear(ab, tb, fti);
    matpoly_ft_clear(ab, tc,  fti);
    mpz_clear(p);
}/*}}}*/
