#include "cado.h"
#include "abase.h"
#include <stdlib.h>
#include <gmp.h>
#include "portability.h"
#include "macros.h"
#include "lingen-matpoly.h"

int matpoly_check_pre_init(matpoly_srcptr p)
{
    if (p->m && p->n && p->alloc)
        return 0;
    if (!p->m && !p->n && !p->alloc)
        return 1;
    abort();
    return 0;
}

/* with the exception of matpoly_realloc, all functions here are exactly
 * identical to those in lingen-polymat.c */
/* {{{ init/zero/clear interface for matpoly */
void matpoly_init(abdst_field ab, matpoly_ptr p, unsigned int m, unsigned int n, int len) {
    memset(p, 0, sizeof(matpoly));
    /* As a special case, we allow a pre-init state with m==n==len==0 */
    ASSERT(!m == !n);
    ASSERT(!m == !len);
    if (!m) return;
    p->m = m;
    p->n = n;
    p->alloc = len;
    p->size = 0;
    abvec_init(ab, &(p->x), m*n*p->alloc);
    abvec_set_zero(ab, p->x, m*n*p->alloc);
}
void matpoly_realloc(abdst_field ab, matpoly_ptr p, size_t newalloc) {
    matpoly_check_pre_init(p);
    /* zero out the newly added data */
    if (newalloc > p->alloc) {
        /* allocate new space, then deflate */
        abvec_reinit(ab, &(p->x), p->m*p->n*p->alloc, p->m*p->n*newalloc);
        abelt * rhead = p->x + p->m*p->n*p->alloc;
        abelt * whead = p->x + p->m*p->n*newalloc;
        for(unsigned int i = p->m ; i-- ; ) {
            for(unsigned int j = p->n ; j-- ; ) {
                whead -= newalloc;
                rhead -= p->alloc;
                abvec_set(ab, whead, rhead, p->alloc);
                abvec_set_zero(ab, whead+p->alloc, newalloc - p->alloc);
            }
        }
    } else {
        /* inflate, then free space */
        ASSERT_ALWAYS(p->size <= newalloc);
        abelt * rhead = p->x;
        abelt * whead = p->x;
        for(unsigned int i = 0 ; i < p->m ; i++) {
            for(unsigned int j = 0 ; j < p->n ; j++) {
                abvec_set(ab, whead, rhead, newalloc);
                rhead += p->alloc;
                whead += newalloc;
            }
        }
        abvec_reinit(ab, &(p->x), p->m*p->n*p->alloc, p->m*p->n*newalloc);
    }
    p->alloc = newalloc;
}
void matpoly_zero(abdst_field ab, matpoly_ptr p) {
    p->size = 0;
    abvec_set_zero(ab, p->x, p->m*p->n*p->alloc);
}
void matpoly_clear(abdst_field ab, matpoly_ptr p) {
    abvec_clear(ab, &(p->x), p->m*p->n*p->alloc);
    memset(p, 0, sizeof(matpoly));
}
void matpoly_swap(matpoly_ptr a, matpoly_ptr b)
{
    matpoly x;
    memcpy(x, a, sizeof(matpoly));
    memcpy(a, b, sizeof(matpoly));
    memcpy(b, x, sizeof(matpoly));
}
/* }}} */

void matpoly_fill_random(abdst_field ab MAYBE_UNUSED, matpoly_ptr a, unsigned int size, gmp_randstate_t rstate)
{
    ASSERT_ALWAYS(size <= a->alloc);
    a->size = size;
    abvec_random(ab, a->x, a->m*a->n*size, rstate);
}


/* polymat didn't have a proper add(), which is somewhat curious */
void matpoly_add(abdst_field ab,
        matpoly c,
        matpoly a,
        matpoly b)
{
    ASSERT_ALWAYS(c->size == a->size);
    ASSERT_ALWAYS(c->size == b->size);
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < a->n ; j++) {
            abvec_add(ab,
                    matpoly_part(c, i, j, 0),
                    matpoly_part(a, i, j, 0),
                    matpoly_part(b, i, j, 0), c->size);
        }
    }
}

void matpoly_multiply_column_by_x(abdst_field ab, matpoly_ptr pi, unsigned int j, unsigned int size)/*{{{*/
{
    ASSERT_ALWAYS((size + 1) <= pi->alloc);
    for(unsigned int i = 0 ; i < pi->m ; i++) {
        memmove(matpoly_part(pi, i, j, 1), matpoly_part(pi, i, j, 0), 
                size * sizeof(abelt));
        abset_ui(ab, matpoly_coeff(pi, i, j, 0), 0);
    }
}/*}}}*/

void matpoly_truncate(abdst_field ab, matpoly_ptr dst, matpoly_ptr src, unsigned int size)/*{{{*/
{
    ASSERT_ALWAYS(size <= src->alloc);
    if (dst == src) {
        /* Never used by the code, so far. We're leaving garbage coeffs
         * on top, could this be a problem ? */
        dst->size = size;
        return;
    }
    if (matpoly_check_pre_init(dst)) {
        matpoly_init(ab, dst, src->m, src->n, size);
    }
    ASSERT_ALWAYS(dst->m == src->m);
    ASSERT_ALWAYS(dst->n == src->n);
    ASSERT_ALWAYS(size <= dst->alloc);
    dst->size = size;
    /* XXX Much more cumbersome here than for polymat, of course */
    for(unsigned int i = 0 ; i < src->m ; i++) {
        for(unsigned int j = 0 ; j < src->n ; j++) {
            abvec_set(ab,
                    matpoly_part(dst, i, j, 0),
                    matpoly_part(src, i, j, 0),
                    size);
        }
    }
}/*}}}*/

/* XXX compared to polymat, our diffferent striding has a consequence,
 * clearly ! */
void matpoly_extract_column(abdst_field ab,/*{{{*/
        matpoly_ptr dst, unsigned int jdst, unsigned int kdst,
        matpoly_ptr src, unsigned int jsrc, unsigned int ksrc)
{
    ASSERT_ALWAYS(dst->m == src->m);
    for(unsigned int i = 0 ; i < src->m ; i++)
        abset(ab,
            matpoly_coeff(dst, i, jdst, kdst),
            matpoly_coeff(src, i, jsrc, ksrc));
}/*}}}*/

void matpoly_extract_row_fragment(abdst_field ab,/*{{{*/
        matpoly_ptr dst, unsigned int i1, unsigned int j1,
        matpoly_ptr src, unsigned int i0, unsigned int j0,
        unsigned int n)
{
    ASSERT_ALWAYS(src->size <= dst->alloc);
    ASSERT_ALWAYS(dst->size == src->size);
    for(unsigned int k = 0 ; k < n ; k++)
        abvec_set(ab,
                matpoly_part(dst, i1, j1 + k, 0),
                matpoly_part(src, i0, j0 + k, 0), dst->size);
}/*}}}*/

void matpoly_rshift(abdst_field ab, matpoly_ptr dst, matpoly_ptr src, unsigned int k)/*{{{*/
{
    ASSERT_ALWAYS(k <= src->size);
    unsigned int newsize = src->size - k;
    if (matpoly_check_pre_init(dst)) {
        matpoly_init(ab, dst, src->m, src->n, newsize);
    }
    ASSERT_ALWAYS(dst->m == src->m);
    ASSERT_ALWAYS(dst->n == src->n);
    ASSERT_ALWAYS(newsize <= dst->alloc);
    dst->size = newsize;
    for(unsigned int i = 0 ; i < src->m ; i++) {
        for(unsigned int j = 0 ; j < src->n ; j++) {
            abvec_set(ab,
                    matpoly_part(dst, i, j, 0),
                    matpoly_part(src, i, j, k),
                    newsize);
        }
    }
    if (dst == src)
        matpoly_realloc(ab, dst, newsize);
}/*}}}*/


void matpoly_addmul(abdst_field ab, matpoly c, matpoly a, matpoly b)/*{{{*/
{
    ASSERT_ALWAYS(a->n == b->m);
    if (matpoly_check_pre_init(c)) {
        matpoly_init(ab, c, a->m, b->n, a->size + b->size - 1);
    }
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    ASSERT_ALWAYS(c->alloc >= a->size + b->size - 1);
    c->size = a->size + b->size - 1;
    abvec_ur tmp[2];
    abvec_ur_init(ab, &tmp[0], c->size);
    abvec_ur_init(ab, &tmp[1], c->size);
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < b->n ; j++) {
            abvec_ur_set_vec(ab, tmp[1], matpoly_part(c, i, j, 0), c->size);
            for(unsigned int k = 0 ; k < a->n ; k++) {
                abvec_conv_ur(ab, tmp[0],
                        matpoly_part(a, i, k, 0), a->size,
                        matpoly_part(b, k, j, 0), b->size);
                abvec_ur_add(ab, tmp[1], tmp[1], tmp[0], c->size);
            }
            abvec_reduce(ab, matpoly_part(c, i, j, 0), tmp[1], c->size);
        }
    }
    abvec_ur_clear(ab, &tmp[0], c->size);
    abvec_ur_clear(ab, &tmp[1], c->size);
}/*}}}*/

void matpoly_mul(abdst_field ab, matpoly c, matpoly a, matpoly b)/*{{{*/
{
    ASSERT_ALWAYS(a->n == b->m);
    if (matpoly_check_pre_init(c)) {
        matpoly_init(ab, c, a->m, b->n, a->size + b->size - 1);
    }
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    ASSERT_ALWAYS(c->alloc >= a->size + b->size - 1);
    c->size = a->size + b->size - 1;
    abvec_set_zero(ab, c->x, c->m * c->n * c->size);
    matpoly_addmul(ab, c, a, b);
}/*}}}*/

void matpoly_addmp(abdst_field ab, matpoly b, matpoly a, matpoly c)/*{{{*/
{
    unsigned int nb = MAX(a->size, c->size) - MIN(a->size, c->size) + 1;
    ASSERT_ALWAYS(a->n == c->m);
    if (matpoly_check_pre_init(b)) {
        matpoly_init(ab, b, a->m, c->n, nb);
    }
    ASSERT_ALWAYS(b->m == a->m);
    ASSERT_ALWAYS(b->n == c->n);
    ASSERT_ALWAYS(b->alloc >= nb);
    b->size = nb;

    /* We are going to make it completely stupid for a beginning. */
    /* XXX XXX XXX FIXME !!! */
    abvec_ur tmp[2];
    abvec_ur_init(ab, &tmp[0], a->size + c->size - 1);
    abvec_ur_init(ab, &tmp[1], b->size);
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < c->n ; j++) {
            abvec_ur_set_vec(ab, tmp[1], matpoly_part(b, i, j, 0), b->size);
            for(unsigned int k = 0 ; k < a->n ; k++) {
                abvec_conv_ur(ab, tmp[0],
                        matpoly_part(a, i, k, 0), a->size,
                        matpoly_part(c, k, j, 0), c->size);
                abvec_ur_add(ab, tmp[1], tmp[1],
                        tmp[0] + MIN(a->size, c->size) - 1,
                        nb);
            }
            abvec_reduce(ab, matpoly_part(b, i, j, 0), tmp[1], b->size);
        }
    }
    abvec_ur_clear(ab, &tmp[0], a->size + c->size - 1);
    abvec_ur_clear(ab, &tmp[1], b->size);
}/*}}}*/

void matpoly_mp(abdst_field ab, matpoly b, matpoly a, matpoly c)/*{{{*/
{
    unsigned int nb = MAX(a->size, c->size) - MIN(a->size, c->size) + 1;
    ASSERT_ALWAYS(a->n == c->m);
    if (matpoly_check_pre_init(b)) {
        matpoly_init(ab, b, a->m, c->n, nb);
    }
    ASSERT_ALWAYS(b->m == a->m);
    ASSERT_ALWAYS(b->n == c->n);
    ASSERT_ALWAYS(b->alloc >= nb);
    b->size = nb;
    abvec_set_zero(ab, b->x, b->m * b->n * b->size);
    matpoly_addmp(ab, b, a, c);
}/*}}}*/



