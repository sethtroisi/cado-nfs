#include "cado.h"
#include <stdlib.h>
#include "macros.h"
#include "utils.h"
#include "lingen-matpoly.h"
#include "lingen-matpoly-ft.h"

/* hackish.
 * when mpfq has a better field_characteristic function, we might
 * consider doing something better.
 */
#ifndef PTR
#define PTR(x) ((x)->_mp_d)
#endif
#ifndef SIZ
#define SIZ(x) ((x)->_mp_size)
#endif
#ifndef ALLOC
#define ALLOC(x) ((x)->_mp_alloc)
#endif


void matpoly_ft_init(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr t, unsigned int m, unsigned int n, struct fft_transform_info * fti)
{
    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    t->m = m;
    t->n = n;
    t->data = malloc(m * n * fft_alloc_sizes[0]);
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

    /* okay, this is downright disgusting */
    /* TODO: mpfq's interface is flawed. We should have access to an
     * mpz_t variable */ 
    mpz_t p;
    SIZ(p) = sizeof(abelt) / sizeof(mp_limb_t);
    ALLOC(p) = sizeof(abelt) / sizeof(mp_limb_t);
    PTR(p) = ab->p;

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    void * tt = malloc(fft_alloc_sizes[1]);

    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            abelt * aij = matpoly_part(a, i, j, 0);
            fft_do_dft_fppol(tij, (mp_limb_t *) aij, a->size, tt, fti, p);
        }
    }
    free(tt);
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

void matpoly_ft_mul(abdst_field ab MAYBE_UNUSED, matpoly_ft_ptr u, matpoly_ft_ptr t0, matpoly_ft_ptr t1, struct fft_transform_info * fti)
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

    for(unsigned int i = 0 ; i < t0->m ; i++) {
        for(unsigned int j = 0 ; j < t1->n ; j++) {
            void * uij  = pointer_arith(u->data, (i*u->n+j) * tsize);
            memset(uij, 0, fft_alloc_sizes[0]);
            fft_transform_prepare(uij, fti);
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

    /* okay, this is downright disgusting */
    /* TODO: mpfq's interface is flawed. We should have access to an
     * mpz_t variable */ 
    mpz_t p;
    SIZ(p) = sizeof(abelt) / sizeof(mp_limb_t);
    ALLOC(p) = sizeof(abelt) / sizeof(mp_limb_t);
    PTR(p) = ab->p;

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    void * tt = malloc(fft_alloc_sizes[1]);

    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            abelt * aij = matpoly_part(a, i, j, 0);
            fft_do_ift_fppol((mp_limb_t *) aij, a->size, tij, tt, fti, p);
        }
    }
    free(tt);
}

void matpoly_ft_ift_mp(abdst_field ab, matpoly_ptr a, matpoly_ft_ptr t, unsigned int shift, struct fft_transform_info * fti)
{
    /* Recall that this is for polynomials ! */
    ASSERT_ALWAYS(t->m == a->m);
    ASSERT_ALWAYS(t->n == a->n);

    /* okay, this is downright disgusting */
    /* TODO: mpfq's interface is flawed. We should have access to an
     * mpz_t variable */ 
    mpz_t p;
    SIZ(p) = sizeof(abelt) / sizeof(mp_limb_t);
    ALLOC(p) = sizeof(abelt) / sizeof(mp_limb_t);
    PTR(p) = ab->p;

    size_t fft_alloc_sizes[3];
    fft_get_transform_allocs(fft_alloc_sizes, fti);
    void * tt = malloc(fft_alloc_sizes[1]);

    for(unsigned int i = 0 ; i < t->m ; i++) {
        for(unsigned int j = 0 ; j < t->n ; j++) {
            size_t offset = (i*t->n+j) * fft_alloc_sizes[0];
            void * tij = pointer_arith(t->data, offset);
            abelt * aij = matpoly_part(a, i, j, 0);
            fft_do_ift_fppol_mp((mp_limb_t *) aij, a->size, tij, tt, fti, p, shift);
        }
    }
    free(tt);
}
