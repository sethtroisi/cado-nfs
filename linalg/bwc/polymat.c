#include "cado.h"
#include "abase.h"
#include <stdlib.h>
#include "macros.h"
#include "polymat.h"

/* {{{ init/zero/clear interface for polymat */
void polymat_init(polymat_ptr p, unsigned int m, unsigned int n, int len) {
    memset(p, 0, sizeof(polymat));
    /* As a special case, we allow a pre-init state with m==n==len==0 */
    ASSERT(!m == !n);
    ASSERT(!m == !len);
    if (!m) return;
    p->m = m;
    p->n = n;
    p->alloc = len;
    p->size = 0;
    p->x = malloc(m*n*p->alloc*sizeof(abelt));
    memset(p->x, 0, m*n*p->alloc*sizeof(abelt));
}
void polymat_realloc(polymat_ptr p, int newalloc) {
    p->x = realloc(p->x, p->m*p->n*newalloc*sizeof(abelt));
    /* We should be copying only from p->size onwards. However our
     * tracking of the ->size field is not accurate. */
    memset(p->x + p->m*p->n*p->alloc, 0, p->m*p->n*(newalloc - p->alloc)*sizeof(abelt));
    p->alloc = newalloc;
}
void polymat_zero(polymat_ptr p) {
    p->size = 0;
    memset(p->x, 0, p->m*p->n*p->alloc*sizeof(abelt));
}
void polymat_clear(polymat_ptr p) {
    free(p->x);
    memset(p, 0, sizeof(polymat));
}
void polymat_swap(polymat_ptr a, polymat_ptr b)
{
    polymat x;
    memcpy(x, a, sizeof(polymat));
    memcpy(a, b, sizeof(polymat));
    memcpy(b, x, sizeof(polymat));
}

/* }}} */

/* {{{ init/zero/clear interface for polymat_ur */
void polymat_ur_init(polymat_ur_ptr p, unsigned int m, unsigned int n, int len) {
    memset(p, 0, sizeof(polymat_ur));
    p->m = m;
    p->n = n;
    p->alloc = len;
    p->size = 0;
    p->x = malloc(m*n*p->alloc*sizeof(abelt_ur));
    memset(p->x, 0, m*n*p->alloc*sizeof(abelt_ur));
}
void polymat_ur_realloc(polymat_ur_ptr p, int newalloc) {
    p->x = realloc(p->x, p->m*p->n*newalloc*sizeof(abelt_ur));
    memset(p->x + p->m*p->n*p->alloc, 0, p->m*p->n*(newalloc - p->alloc)*sizeof(abelt_ur));
    p->alloc = newalloc;
}
void polymat_ur_zero(polymat_ur_ptr p) {
    p->size = 0;
    memset(p->x, 0, p->m*p->n*p->alloc*sizeof(abelt_ur));
}
void polymat_ur_clear(polymat_ur_ptr p) {
    free(p->x);
    memset(p, 0, sizeof(polymat_ur));
}
void polymat_ur_swap(polymat_ur_ptr a, polymat_ur_ptr b)
{
    polymat_ur x;
    memcpy(x, a, sizeof(polymat_ur));
    memcpy(a, b, sizeof(polymat_ur));
    memcpy(b, x, sizeof(polymat_ur));
}

/* }}} */


void bwmat_copy_coeffs(abdst_field ab MAYBE_UNUSED, abelt * x0, int stride0, abelt * x1, int stride1, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        abset(ab, *x0, *x1);
        x0+=stride0;
        x1+=stride1;
    }
}
void bwmat_move_coeffs(abdst_field ab MAYBE_UNUSED, abelt * x0, int stride0, abelt * x1, int stride1, unsigned int n)
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

void polymat_mul(abdst_field ab, polymat c, polymat a, polymat b)
{
    polymat_ur cx;
    ASSERT_ALWAYS(a->n == b->m);
    if (!c->m) {
        ASSERT_ALWAYS(!c->m && !c->n && !c->size);
        polymat_init(c, a->m, b->n, a->size + b->size - 1);
    }
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    ASSERT_ALWAYS(c->alloc >= a->size + b->size - 1);

    polymat_ur_init(cx, c->m, c->n, 1);
    abelt_ur tmp_ur;
    abelt_ur_init(ab, &tmp_ur);
    c->size = a->size + b->size - 1;
    /* asize = dega + 1
     * bsize = degb + 1
     */
    for(unsigned int s = 0 ; s < c->size ; s++) {
        /* s < degc + 1, s<=degc = dega + degb */
        polymat_ur_zero(cx);
        unsigned int t0 = 0, t1 = s + 1;
        if (s - t0 >= b->size)          /* s-t0 > degb */
            t0 = s - (b->size - 1);     /* s-t0 = degb */
                                        /* NB: s-degb <= dega */
        if (t1 >= a->size)
            t1 = a->size;               /* t1 = dega + 1 */
                                        /* NB: s-dega <= degb */
        for(unsigned int t = t0 ; t < t1 ; t++) {
            unsigned int u = s - t;
            for(unsigned int i = 0 ; i < a->m ; i++) {
                for(unsigned int j = 0 ; j < b->n ; j++) {
                    for(unsigned int k = 0 ; k < a->n ; k++) {
                        abmul_ur(ab, tmp_ur,
                                polymat_coeff(a, i, k, t),
                                polymat_coeff(b, k, j, u));
                        abelt_ur_add(ab,
                                polymat_ur_coeff(cx, i, j, 0),
                                polymat_ur_coeff(cx, i, j, 0),
                                tmp_ur);
                    }
                }
            }
        }
        for(unsigned int i = 0 ; i < a->m ; i++) {
            for(unsigned int j = 0 ; j < b->n ; j++) {
                abreduce(ab, 
                        polymat_coeff(c, i, j, s),
                        polymat_ur_coeff(cx, i, j, 0));
            }
        }

    }
    polymat_ur_clear(cx);
    abelt_ur_clear(ab, &tmp_ur);
}

/* Compute the middle product of a and b. This is the part of the product
 * where there are the largest number of summands used to acccumulate the
 * terms of the result.
 *
 * Note that the first coefficient computed is stored as the coefficient
 * of degree 0 in c, even though it corresponds to something of larger
 * degree in a*b.
 *
 * MAX(s-(bsize-1), 0) <= t < MIN(a->size, s+1)
 * MIN(bsize, s+1) > s-t >= MAX(s-asize+1,0)
 * 
 * We want exactly MIN(asize, bsize) summands to reach degree s. These
 * correspond to indices:
 *
 * (s,0)....(s-MIN(asize,bsize)+1,MIN(asize,bsize)-1)
 *
 * We must have s>=MIN(asize,bsize)-1, otherwise the number of summands
 * will not be large enough. And this number of summands is possible only
 * as long as no summand goes off bounds for a or b. This implies that
 * the coefficient of degree 0 of one of a or b has to appear in all
 * products, and that the condition is thus s <= MAX(asize,bsize) - 1
 *
 * IOW:
 *      MIN(asize, bsize)-1 <= s < MAX(asize, bsize)
 *
 * which gives the following examples:
 *
 *      asize==bsize==n+1 --> s==n
 *      asize=n+1, bsize=2n --> n<=s<2n.
 */
void polymat_mp(abdst_field ab, polymat c, polymat a, polymat b)
{
    unsigned int s0 = MIN(a->size, b->size) - 1;
    unsigned int s1 = MAX(a->size, b->size);
    ASSERT_ALWAYS(a->n == b->m);
    if (!c->m) {
        ASSERT_ALWAYS(!c->m && !c->n && !c->size);
        polymat_init(c, c->m, c->n, s1 - s0);
    }
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    ASSERT_ALWAYS(c->alloc >= s1 - s0);
    c->size = s1 - s0;

    polymat_mp_raw(ab, c, 0, a, 0, a->size, b, 0, b->size);
}

/* This one is more generic, and allows arbitrary middle products. It
 * disregards a->size and b->size, and looks at the number of
 * coefficients na and nb, embodied within a and b at offsets xa and xb,
 * respectively. The coefficients of the result (their number is
 * MAX(na,nb)-MIN(na,nb) + 1) are put at offset xc within c
 */
void polymat_mp_raw(abdst_field ab,
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb)
{
    unsigned int s0 = MIN(na, nb) - 1;
    unsigned int s1 = MAX(na, nb);
    polymat_ur cx;
    ASSERT_ALWAYS(a->n == b->m);
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    ASSERT_ALWAYS(c->alloc >= s1 - s0);

    polymat_ur_init(cx, c->m, c->n, 1);
    abelt_ur tmp_ur;
    abelt_ur_init(ab, &tmp_ur);
    ASSERT_ALWAYS(s1 - s0 == MAX(na, nb) - MIN(na, nb) + 1);
    for(unsigned int s = s0 ; s < s1 ; s++) {
        polymat_ur_zero(cx);
        cx->size = 1;
        unsigned int t0 = 0;
        unsigned int t1 = MIN(na, nb);
        if (nb < na) {
            t0 += s - (nb - 1);
            t1 += s - (nb - 1);
        }
        ASSERT_ALWAYS(t0 < na);
        ASSERT_ALWAYS((t1-1) < na);
        ASSERT_ALWAYS(s-t0 < nb);
        ASSERT_ALWAYS(s-(t1-1) < nb);
        ASSERT_ALWAYS(s >= t0);
        ASSERT_ALWAYS(s >= (t1-1));
        // printf("Coefficient of degree %u in product\n", xa + xb + s);
        for(unsigned int t = t0 ; t < t1 ; t++) {
            unsigned int u = s - t;
            // printf("a%u * b%u\n", xa + t, xb + u);
            for(unsigned int i = 0 ; i < a->m ; i++) {
                for(unsigned int j = 0 ; j < b->n ; j++) {
                    for(unsigned int k = 0 ; k < a->n ; k++) {
                        abmul_ur(ab, tmp_ur,
                                polymat_coeff(a, i, k, xa + t),
                                polymat_coeff(b, k, j, xb + u));
                        abelt_ur_add(ab,
                                polymat_ur_coeff(cx, i, j, 0),
                                polymat_ur_coeff(cx, i, j, 0),
                                tmp_ur);
                    }
                }
            }
        }
        for(unsigned int i = 0 ; i < a->m ; i++) {
            for(unsigned int j = 0 ; j < b->n ; j++) {
                abreduce(ab, 
                        polymat_coeff(c, i, j, xc + s - s0),
                        polymat_ur_coeff(cx, i, j, 0));
            }
        }
    }
    polymat_ur_clear(cx);
    abelt_ur_clear(ab, &tmp_ur);
}
