#include "cado.h"
#include "abase.h"
#include <stdlib.h>
#include "portability.h"
#include "macros.h"
#include "polymat.h"

#define MIDDLE_PRODUCT_KARA_THRESHOLD   10
#define MIDDLE_PRODUCT_SUBDIVIDE_THRESHOLD      MIDDLE_PRODUCT_KARA_THRESHOLD
#define PRODUCT_KARA_THRESHOLD   10
#define PRODUCT_SUBDIVIDE_THRESHOLD      PRODUCT_KARA_THRESHOLD


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
void bwmat_zero_coeffs(abdst_field ab MAYBE_UNUSED, abelt * x0, int stride0, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        abset_zero(ab, *x0);
        x0+=stride0;
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

/* {{{ Elementary operations on the matrices, which are the coefficients
 * of our polynomials */
void polymat_addmat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb)
{
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < a->n ; j++) {
            abadd(ab,
                    polymat_coeff(c, i, j, kc),
                    polymat_coeff(a, i, j, ka),
                    polymat_coeff(b, i, j, kb));
        }
    }
}

void polymat_submat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb)
{
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < a->n ; j++) {
            absub(ab,
                    polymat_coeff(c, i, j, kc),
                    polymat_coeff(a, i, j, ka),
                    polymat_coeff(b, i, j, kb));
        }
    }
}

void polymat_addmulmat_ur(abdst_field ab,
        polymat_ur c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb)
{
    abelt_ur tmp_ur;
    abelt_ur_init(ab, &tmp_ur);
    ASSERT_ALWAYS(a->n == b->m);
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < b->n ; j++) {
            for(unsigned int k = 0 ; k < a->n ; k++) {
                absrc_elt ea = polymat_coeff(a, i, k, ka);
                absrc_elt eb = polymat_coeff(b, k, j, kb);
                abmul_ur(ab, tmp_ur, ea, eb);
                abelt_ur_add(ab,
                        polymat_ur_coeff(c, i, j, kc),
                        polymat_ur_coeff(c, i, j, kc),
                        tmp_ur);
            }
        }
    }
    abelt_ur_clear(ab, &tmp_ur);
}

void polymat_reducemat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat_ur a, unsigned int ka)
{
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < a->n ; j++) {
            abreduce(ab, 
                    polymat_coeff(c, i, j, kc),
                    polymat_ur_coeff(a, i, j, ka));
        }
    }
}

void polymat_mulmat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb)
{
    polymat_ur cx;
    polymat_ur_init(cx, c->m, c->n, 1);
    polymat_addmulmat_ur(ab, cx, 0, a, ka, b, kb);
    polymat_reducemat(ab, c, kc, cx, 0);
    polymat_ur_clear(cx);
}

void polymat_addmulmat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb)
{
    polymat cc;
    polymat_init(cc, c->m, c->n, 1);
    polymat_mulmat(ab,cc,0,a,ka,b,kb);
    polymat_addmat(ab,c,kc,c,kc,cc,0);
    polymat_clear(cc);
}

/* }}} */

void polymat_mul_raw_basecase(abdst_field ab,/*{{{*/
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb,
        int transpose, int add)
{
    unsigned int nc = na + nb - 1;
    ASSERT_ALWAYS(c->alloc >= xc + nc);

    polymat_ur tmat_ur;
    polymat_ur_init(tmat_ur, c->m, c->n, 1);

    for(unsigned int k = 0 ; k < nc ; k++) {
        unsigned int i0 = k >= nb ? k + 1 - nb : 0;
        unsigned int i1 = k + 1 < na ? k + 1 : na;
        polymat_ur_zero(tmat_ur);
        for(unsigned int i = i0 ; i < i1 ; i++) {
            unsigned int j = k - i;
            if (!transpose) {
                polymat_addmulmat_ur(ab, tmat_ur, 0, a, xa + i, b, xb + j);
            } else {
                polymat_addmulmat_ur(ab, tmat_ur, 0, b, xb + j, a, xa + i);
            }
        }
        if (!add) {
            polymat_reducemat(ab, c, xc + k, tmat_ur, 0);
        } else {
            polymat tmat;
            polymat_init(tmat, c->m, c->n, 1);
            polymat_reducemat(ab, tmat, 0, tmat_ur, 0);
            polymat_addmat(ab, c, xc + k, c, xc + k, tmat, 0);
            polymat_clear(tmat);
        }
    }
    polymat_ur_clear(tmat_ur);
}/*}}}*/

void polymat_mul_raw_kara(abdst_field ab,/*{{{*/
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb,
        int transpose, int add)
{
    if (MIN(na, nb) < PRODUCT_KARA_THRESHOLD) {
        polymat_mul_raw_basecase(ab, c, xc, a, xa, na, b, xb, nb, transpose, add);
        return;
    }
    unsigned int nc = na + nb - 1;
    ASSERT_ALWAYS(c->alloc >= xc + nc);
    if (!add) {
        bwmat_zero_coeffs(ab, polymat_part(c, 0, 0, xc), 1, c->m*c->n*nc);
    }
    /* This code works __only__ for kara-sized products */
    ASSERT_ALWAYS(nb == na);
    unsigned int m0 = na / 2;
    unsigned int m1 = na - m0;

    polymat q0, q2, s, t;
    polymat_init(s, a->m, a->n, m1);
    polymat_init(t, b->m, b->n, m1);
    for(unsigned int k = 0 ; k < m0 ; k++) {
        polymat_addmat(ab, s, k, a, xa + k, a, xa + m0 + k);
        polymat_addmat(ab, t, k, b, xb + k, b, xb + m0 + k);
    }
    if (m1 > m0) {
        polymat_addmat(ab, s, m0, s, m0, a, xa + m0 + m0);
        polymat_addmat(ab, t, m0, t, m0, b, xb + m0 + m0);
    }
    polymat_mul_raw_kara(ab, c, xc+m0, s, 0, m1, t, 0, m1, transpose, 1);
    polymat_clear(s);
    polymat_clear(t);

    polymat_init(q0, c->m, c->n, 2*m0-1);
    polymat_init(q2, c->m, c->n, 2*m1-1);
    polymat_mul_raw_kara(ab, q0, 0, a, xa, m0, b, xb, m0, transpose, 0);
    polymat_mul_raw_kara(ab, q2, 0, a, xa+m0, m1, b, xb+m0, m1, transpose, 0);
    for(unsigned int k = 0 ; k < 2*m0 - 1 ; k++) {
        polymat_addmat(ab, c, xc + k, c, xc + k, q0, k);
        polymat_submat(ab, c, xc + m0 + k, c, xc + m0 + k, q0, k);
    }
    for(unsigned int k = 0 ; k < 2*m1 - 1 ; k++) {
        polymat_addmat(ab, c, xc + 2 * m0 + k, c, xc + 2 * m0 + k, q2, k);
        polymat_submat(ab, c, xc + m0 + k, c, xc + m0 + k, q2, k);
    }
    polymat_clear(q0);
    polymat_clear(q2);
} /* }}} */

void polymat_mul_raw_subdivide(abdst_field ab,/*{{{*/
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb,
        int transpose, int add)
{
    if (!na || !nb) {
        return;
    }
    if (MIN(na, nb) < PRODUCT_SUBDIVIDE_THRESHOLD) {
        polymat_mul_raw_basecase(ab, c, xc, a, xa, na, b, xb, nb, transpose, add);
        return;
    }
    unsigned int nc = na + nb - 1;
    if (!add) {
        bwmat_zero_coeffs(ab, polymat_part(c, 0, 0, xc), 1, c->m*c->n*nc);
    }

    for( ; nb >= na ; ) {
        polymat_mul_raw_kara(ab, c, xc,
                a, xa, na, b, xb, na, transpose, 1);
        xb += na;
        xc += na;
        nb -= na;
    }
    if (!nb) {
        /* This condition also catches the Karatsuba case n, n */
        return;
    }
    ASSERT_ALWAYS(nb < na);
    /* Fixup needed, now. Treat the largest possible subblock. We want
     * the intervals [0,na-k[ and [0, nb[ to match a kara scheme, i.e.
     * have k=nb-na.
     */
    unsigned int chop = na - nb;
    polymat_mul_raw_kara(ab, c, xc,
            a, xa, nb, b, xb, nb, transpose, 1);
    polymat_mul_raw_subdivide(ab, c, xc + na - chop,
            a, xa + nb, chop, b, xb, nb, transpose, 1);
}/*}}}*/

void polymat_mul_raw(abdst_field ab,/*{{{*/
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb,
        int transpose, int add)
{
    polymat_mul_raw_subdivide(ab, c, xc, a, xa, na, b, xb, nb, transpose, add);
}/*}}}*/

void polymat_mul(abdst_field ab, polymat c, polymat a, polymat b)/*{{{*/
{
    ASSERT_ALWAYS(a->n == b->m);
    if (!c->m) {
        ASSERT_ALWAYS(!c->m && !c->n && !c->size);
        polymat_init(c, a->m, b->n, a->size + b->size - 1);
    }
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    ASSERT_ALWAYS(c->alloc >= a->size + b->size - 1);
    c->size = a->size + b->size - 1;
    polymat_mul_raw(ab, c, 0, a, 0, a->size, b, 0, b->size, 0, 0);
}/*}}}*/

/* Middle product b of a and c. This is the part of the product where
 * there are the largest number of summands used to acccumulate the terms
 * of the result.
 *
 * Note that the first coefficient computed is stored as the coefficient
 * of degree 0 in b, even though it corresponds to something of larger
 * degree in a*c.
 *
 * We require a to be the small operand, and c the larger one (we support
 * the other case, though, by swapping arguments).
 *
 * Exactly nc-na+1 coefficients are stored in b (so that nc = na + nb -
 * 1). Those are the coefficients of degrees [na-1..nc[ in the product
 * a*c.
 *
 * Some examples:
 *
 *      na==nc==n+1 --> nb==1: [n]
 *      na=n+1, nc=2n --> nb==n: [n,2n[
 */
void polymat_mp_raw_basecase(abdst_field ab,/*{{{*/
        polymat b, unsigned int xb,
        polymat a, unsigned int xa, unsigned int na,
        polymat c, unsigned int xc, unsigned int nc,
        int transpose, int add)
{
    if (nc < na) {
        polymat_mp_raw_basecase(ab, b, xb, c, xc, nc, a, xa, na, !transpose, add);
        return;
    }
    unsigned int nb = nc - na + 1;
    ASSERT_ALWAYS(b->alloc >= xb + nb);

    polymat_ur tmat_ur;
    polymat_ur_init(tmat_ur, b->m, b->n, 1);

    for(unsigned int j = 0 ; j < nb ; j++) {
        polymat_ur_zero(tmat_ur);
        for(unsigned int i = 0 ; i < na ; i++) {
            unsigned int k = j + na - 1 - i;
            if (!transpose) {
                polymat_addmulmat_ur(ab, tmat_ur, 0, a, xa + i, c, xc + k);
            } else {
                polymat_addmulmat_ur(ab, tmat_ur, 0, c, xc + k, a, xa + i);
            }
        }
        if (!add) {
            polymat_reducemat(ab, b, xb + j, tmat_ur, 0);
        } else {
            polymat tmat;
            polymat_init(tmat, b->m, b->n, 1);
            polymat_reducemat(ab, tmat, 0, tmat_ur, 0);
            polymat_addmat(ab, b, xb + j, b, xb + j, tmat, 0);
            polymat_clear(tmat);
        }
    }
    polymat_ur_clear(tmat_ur);
}/*}}}*/

void polymat_mp_raw_kara(abdst_field ab,/*{{{*/
        polymat b, unsigned int xb,
        polymat a, unsigned int xa, unsigned int na,
        polymat c, unsigned int xc, unsigned int nc,
        int transpose, int add)
{
    // printf("mp_kara(%d,%d)%b\n",na,nc,transpose?'T':' ');
    if (MIN(na, nc) < MIDDLE_PRODUCT_KARA_THRESHOLD) {
        polymat_mp_raw_basecase(ab, b, xb, a, xa, na, c, xc, nc, transpose, add);
        return;
    }
    if (nc < na) {
        polymat_mp_raw_kara(ab, b, xb, c, xc, nc, a, xa, na, !transpose, add);
        return;
    }
    /* This code works __only__ for kara-sized middle products */
    ASSERT_ALWAYS(nc == 2*na - 1);
    unsigned int m0 = na / 2;
    unsigned int m1 = na - m0;
    /* Spans of different chunks, in (offset, length) format */
    unsigned int span_a0[2] =  {xa,        m0};
    unsigned int span_a1[2] =  {xa+m0,     m1};
    unsigned int span_c0[2] =  {xc,      2*m1-1};
    unsigned int span_c10[2] = {xc+m1,   2*m0-1};
    unsigned int span_c11[2] = {xc+m1,   2*m1-1};
    unsigned int span_c2[2] =  {xc+2*m1, 2*m0-1};
    /* Dammit, we need some temp storage, of course */
    polymat q1, s, t;
    /* The polymats s and t are used to build temporaries, respectively
     * from a and c, of respective maximum lengths m1 and 2*m1-1 */
    /* q0 is MP(a1, b0+b11), i.e. MP(m1, 2*m1-1), which
     * produces m1 coefficients. This goes to offset 0 (well, xb) in b. */
    polymat_init(t, c->m, c->n, 2*m1-1);
    for(unsigned int k = 0 ; k < 2*m1 - 1 ; k++) {
        polymat_addmat(ab, t, k, c, span_c0[0] + k, c, span_c11[0] + k);
    }
    polymat_mp_raw_kara(ab, b, xb,
            a, span_a1[0], span_a1[1],
            t, 0, 2*m1-1,
            transpose, add);
    polymat_clear(t);
    /* q1 is MP(a1-a0, b11), i.e. MP(m1, 2*m1-1) again */
    polymat_init(q1, b->m, b->n, m1);
    polymat_init(s, a->m, a->n, m1);
    bwmat_copy_coeffs(ab,
            polymat_part(s,0,0,0),1,
            polymat_part(a,0,0,span_a1[0]),1,
            m1*a->m*a->n);
    for(unsigned int k = 0 ; k < m0 ; k++) {
        /* We've already copied the a1 coefficients above */
        polymat_submat(ab,
                s, k + m1 - m0,
                s, k + m1 - m0,
                a, span_a0[0] + k);
    }
    polymat_mp_raw_kara(ab, q1, 0,
            s, 0, m1,
            c, span_c11[0], span_c11[1],
            transpose, 0);
    polymat_clear(s);
    /* q2 is MP(a0, b10+b2), i.e. MP(m0, 2*m0-1) */
    /* This goes to offset m1 in b */
    polymat_init(t, c->m, c->n, 2*m0-1);
    for(unsigned int k = 0 ; k < 2*m0 - 1 ; k++) {
        polymat_addmat(ab, t, k, c, span_c10[0] + k, c, span_c2[0] + k);
    }
    polymat_mp_raw_kara(ab, b, xb + m1,
            a, span_a0[0], span_a0[1],
            t, 0, 2*m0-1,
            transpose, add);
    polymat_clear(t);

    /* We now have to append the coefficients to form the result. Not all
     * coefficients in q1 are read for the high part of the middle
     * product. */
    {
        for(unsigned int k = 0 ; k < m1 ; k++) {
            polymat_submat(ab, b, xb + k, b, xb + k, q1, k);
        }
        for(unsigned int k = 0 ; k < m0 ; k++) {
            polymat_addmat(ab, b, xb + m1 + k, b, xb + m1 + k, q1, k);
        }
    }
    polymat_clear(q1);
}/*}}}*/

void polymat_mp_raw_subdivide(abdst_field ab,/*{{{*/
        polymat b, unsigned int xb,
        polymat a, unsigned int xa, unsigned int na,
        polymat c, unsigned int xc, unsigned int nc,
        int transpose, int add)
{
    // printf("mp_subdivide(%d,%d)%b\n",na,nc,transpose?'T':' ');
    /* TODO. Presently the strategy below uses only areas of balanced
     * size (karatsuba-friendly). This is possibly inefficient, if for
     * example nc/na on input is classically close to a ratio which is
     * much different. In block Wiedemann for example, nc/na is typically
     * 1+2n/m, which may well be 3 for instance if m=n=1. In such a case,
     * the present code will reduce our mp to 2*3 = 6 balanced middle
     * products.
     *
     * Other options could be to rely on e.g. a polymat_mp_raw_toom24,
     * which would use the following operations to fall back on only 5
     * balanced middle products:
     *
     * Formulae for MP(2k, 6k)-> 3k
     *  t0 = 1/2*(2*x0 - x1 - 2*x2 + x3)*y1
     *  t1 = 1/6*(y0 - y1)*(2*x1 - 3*x2 + x3)
     *  t2 = (2*x1 - x2 - 2*x3 + x4)*y0
     *  t3 = 1/2*(y0 + y1)*(2*x1 + x2 - x3)
     *  t4 = -1/6*(x1 - x3)*(2*y0 + y1)
     *  z0 = t0 + t1 + t3 + t4
     *  z1 = -t1 + t3 + 2*t4
     *  z2 = t1 + t3 + 4*t4
     *  z3 = -t1 + t2 + t3 + 8*t4
     *
     * Another option, slightly less efficient, uses a
     * polymat_mp_raw_toom23 to fall back on balanced middle products in
     * 2 recursion steps. The used formulae, which absorb a ratio nc/na
     * close to 2.5, are as follows:
     *
     * Formulae for MP(2k, 5k)-> 3k
     *  t0 (x0 - x2)*y1
     *  t1 1/2*(x1 - x2)*(y0 - y1)
     *  t2 -(x1 - x3)*y0
     *  t3 1/2*(x1 + x2)*(y0 + y1)
     *  z0 t0 + t1 + t3
     *  z1 -t1 + t3
     *  z2 t1 + t2 + t3
     *
     * The different options, for reaching balanced products of size
     * k/4,2*(k/4) from something of size k,3k are:
     *  - as we do here: Karatsuba twice, two steps: 2*3*3=18M
     *  - using toom24, followed by Karatsuba: 5*3=15M
     *  - using toom23 twice: 4*4=16M.
     */
    if (!na || !nc) {
        return;
    }
    if (MIN(na, nc) < MIDDLE_PRODUCT_SUBDIVIDE_THRESHOLD) {
        polymat_mp_raw_basecase(ab, b, xb, a, xa, na, c, xc, nc, transpose, add);
        return;
    }
    if (nc < na) {
        polymat_mp_raw_subdivide(ab, b, xb, c, xc, nc, a, xa, na, !transpose, add);
        return;
    }
    ASSERT(nc >= na);
    for( ; nc >= 2 * na - 1 ; ) {
        polymat_mp_raw_kara(ab, b, xb,
                a, xa, na, c, xc, 2*na-1, transpose, add);
        xc += na;
        xb += na;
        nc -= na;
    }
    if (nc < na) {
        /* This condition also catches the Karatsuba case n, 2n-1 */
        return;
    }
    // printf("mp_subdivide_tail(%d,%d)%b\n",na,nc,transpose?'T':' ');
    ASSERT_ALWAYS(nc >= na);
    /* Fixup needed, now. Treat the largest possible subblock. We want
     * the intervals [0,na-k[ and [k, nc[ to match a kara scheme, i.e.
     * have the constraints 0<=k<=na, and nc-k = 2*(na-k)-1, i.e. k =
     * 2*na-1-nc.  By assumption, the expression above, with nc>=na,
     * guarantees that the first condition 0<=k<=na is satisfied.
     */
    unsigned int chop = 2*na-1-nc;
    polymat_mp_raw_kara(ab, b, xb,
            a, xa, na - chop, c, xc + chop, nc - chop, transpose, add);
    polymat_mp_raw_subdivide(ab, b, xb,
            a, xa + na - chop, chop, c, xc, na - 1, transpose, 1);
}/*}}}*/

void polymat_mp_raw(abdst_field ab,/*{{{*/
        polymat b, unsigned int xb,
        polymat a, unsigned int xa, unsigned int na,
        polymat c, unsigned int xc, unsigned int nc,
        int transpose, int add)
{
    polymat_mp_raw_subdivide(ab, b, xb, a, xa, na, c, xc, nc, transpose, add);
}/*}}}*/

void polymat_mp(abdst_field ab, polymat b, polymat a, polymat c)/*{{{*/
{
    unsigned int nb = MAX(a->size, c->size) - MIN(a->size, c->size) + 1;
    ASSERT_ALWAYS(a->n == c->m);
    if (!b->m) {
        ASSERT_ALWAYS(!b->m && !b->n && !b->size);
        polymat_init(b, b->m, b->n, nb);
    }
    ASSERT_ALWAYS(b->m == a->m);
    ASSERT_ALWAYS(b->n == c->n);
    ASSERT_ALWAYS(b->alloc >= nb);
    b->size = nb;

    polymat_mp_raw(ab, b, 0, a, 0, a->size, c, 0, c->size, 0, 0);
}/*}}}*/

