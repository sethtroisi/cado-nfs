#include "cado.h"
#include "abase.h"
#include <stdlib.h>
#include "macros.h"
#include "polymat.h"

#define MIDDLE_PRODUCT_KARA_THRESHOLD   8
#define MIDDLE_PRODUCT_SUBDIVIDE_THRESHOLD      MIDDLE_PRODUCT_KARA_THRESHOLD


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
void polymat_mp_raw_basecase(abdst_field ab,
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb,
        int transpose, int add)
{
    // printf("mp_basecase(%d,%d)%c\n",na,nb,transpose?'T':' ');
    unsigned int s0 = MIN(na, nb) - 1;
    unsigned int s1 = MAX(na, nb);
    polymat_ur cx;
    if (!transpose) {
        ASSERT_ALWAYS(a->n == b->m);
        ASSERT_ALWAYS(c->m == a->m);
        ASSERT_ALWAYS(c->n == b->n);
    } else {
        ASSERT_ALWAYS(b->n == a->m);
        ASSERT_ALWAYS(c->m == b->m);
        ASSERT_ALWAYS(c->n == a->n);
    }
    ASSERT_ALWAYS(c->alloc >= s1 - s0);

    polymat_ur_init(cx, c->m, c->n, 1);
    abelt_ur tmp_ur;
    abelt_ur_init(ab, &tmp_ur);
    ASSERT_ALWAYS(s1 - s0 == MAX(na, nb) - MIN(na, nb) + 1);
    for(unsigned int s = s0 ; s < s1 ; s++) {
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
        polymat_ur_zero(cx);
        cx->size = 1;
        for(unsigned int t = t0 ; t < t1 ; t++) {
            unsigned int u = s - t;
            // printf("a%u * b%u\n", xa + t, xb + u);
            if (!transpose) {
                polymat_addmulmat_ur(ab, cx, 0, a, xa + t, b, xb + u);
            } else {
                polymat_addmulmat_ur(ab, cx, 0, b, xb + u, a, xa + t);
            }
        }
        if (!add) {
            polymat_reducemat(ab, c, xc + s - s0, cx, 0);
        } else {
            polymat cc;
            polymat_init(cc, c->m, c->n, 1);
            polymat_reducemat(ab, cc, 0, cx, 0);
            polymat_addmat(ab, c, xc + s - s0, c, xc + s - s0, cc, 0);
            polymat_clear(cc);
        }
    }
    polymat_ur_clear(cx);
    abelt_ur_clear(ab, &tmp_ur);
}

void polymat_mp_raw_kara(abdst_field ab,
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb,
        int transpose, int add)
{
    // printf("mp_kara(%d,%d)%c\n",na,nb,transpose?'T':' ');
    if (MIN(na, nb) < MIDDLE_PRODUCT_KARA_THRESHOLD) {
        polymat_mp_raw_basecase(ab, c, xc, a, xa, na, b, xb, nb, transpose, add);
        return;
    }
    if (nb < na) {
        polymat_mp_raw_kara(ab, c, xc, b, xb, nb, a, xa, na, !transpose, add);
    }
    /* This code works __only__ for kara-sized middle products */
    ASSERT_ALWAYS(nb == 2*na - 1);
    unsigned int m0 = na / 2;
    unsigned int m1 = na - m0;
    /* Spans of different chunks, in (offset, length) format */
    unsigned int span_a0[2] =  {xa,        m0};
    unsigned int span_a1[2] =  {xa+m0,     m1};
    unsigned int span_b0[2] =  {xb,      2*m1-1};
    unsigned int span_b10[2] = {xb+m1,   2*m0-1};
    unsigned int span_b11[2] = {xb+m1,   2*m1-1};
    unsigned int span_b2[2] =  {xb+2*m1, 2*m0-1};
    /* Dammit, we need some temp storage, of course */
    polymat q0, q1, q2, s, t;
    /* The polymats s and t are used to build temporaries, respectively
     * from a and b, of respective maximum lengths m1 and 2*m1-1 */
    /* q0 is MP(a1, b0+b11), i.e. MP(m1, 2*m1-1), which
     * produces m1 coefficients */
    /* TODO: one can afford _not_ allocating q0, and update c directly */
    polymat_init(q0, c->m, c->n, m1);
    polymat_init(t, b->m, b->n, 2*m1-1);
    for(unsigned int k = 0 ; k < 2*m1 - 1 ; k++) {
        polymat_addmat(ab, t, k, b, span_b0[0] + k, b, span_b11[0] + k);
    }
    polymat_mp_raw_kara(ab, q0, 0,
            a, span_a1[0], span_a1[1],
            t, 0, 2*m1-1,
            transpose, 0);
    polymat_clear(t);
    /* q1 is MP(a1-a0, b11), i.e. MP(m1, 2*m1-1) again */
    polymat_init(q1, c->m, c->n, m1);
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
            b, span_b11[0], span_b11[1],
            transpose, 0);
    polymat_clear(s);
    /* q2 is MP(a0, b10+b2), i.e. MP(m0, 2*m0-1) */
    /* TODO: one can afford _not_ allocating q2, and update c directly */
    polymat_init(q2, c->m, c->n, m1);
    polymat_init(t, b->m, b->n, 2*m0-1);
    for(unsigned int k = 0 ; k < 2*m0 - 1 ; k++) {
        polymat_addmat(ab, t, k, b, span_b10[0] + k, b, span_b2[0] + k);
    }
    polymat_mp_raw_kara(ab, q2, 0,
            a, span_a0[0], span_a0[1],
            t, 0, 2*m0-1,
            transpose, 0);
    polymat_clear(t);

    /* We now have to append the coefficients to form the result. Not all
     * coefficients in q1 are read for the high part of the middle
     * product. */
    if (!add) {
        for(unsigned int k = 0 ; k < m1 ; k++) {
            polymat_submat(ab, c, xc + k, q0, k, q1, k);
        }
        for(unsigned int k = 0 ; k < m0 ; k++) {
            polymat_addmat(ab, c, xc + m1 + k, q2, k, q1, k);
        }
    } else {
        for(unsigned int k = 0 ; k < m1 ; k++) {
            polymat_addmat(ab, c, xc + k, c, xc + k, q0, k);
            polymat_submat(ab, c, xc + k, c, xc + k, q1, k);
        }
        for(unsigned int k = 0 ; k < m0 ; k++) {
            polymat_addmat(ab, c, xc + m1 + k, c, xc + m1 + k, q2, k);
            polymat_addmat(ab, c, xc + m1 + k, c, xc + m1 + k, q1, k);
        }
    }
    polymat_clear(q0);
    polymat_clear(q1);
    polymat_clear(q2);
}

void polymat_mp_raw_subdivide(abdst_field ab,
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb,
        int transpose, int add)
{
    // printf("mp_subdivide(%d,%d)%c\n",na,nb,transpose?'T':' ');
    /* TODO. Presently the strategy below uses only areas of balanced
     * size (karatsuba-friendly). This is possibly inefficient, if for
     * example nb/na on input is classically close to a ratio which is
     * much different. In block Wiedemann for example, nb/na is typically
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
     * 2 recursion steps. The used formulae, which absorb a ratio nb/na
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
    if (!na || !nb) {
        return;
    }
    if (MIN(na, nb) < MIDDLE_PRODUCT_SUBDIVIDE_THRESHOLD) {
        polymat_mp_raw_basecase(ab, c, xc, a, xa, na, b, xb, nb, transpose, add);
        return;
    }
    if (nb < na) {
        polymat_mp_raw_subdivide(ab, c, xc, b, xb, nb, a, xa, na, !transpose, add);
        return;
    }
    ASSERT(nb >= na);
    for( ; nb >= 2 * na - 1 ; ) {
        polymat_mp_raw_kara(ab, c, xc,
                a, xa, na, b, xb, 2*na-1, transpose, add);
        xb += na;
        xc += na;
        nb -= na;
    }
    if (nb < na) {
        /* This condition also catches the Karatsuba case n, 2n-1 */
        return;
    }
    // printf("mp_subdivide_tail(%d,%d)%c\n",na,nb,transpose?'T':' ');
    ASSERT_ALWAYS(nb >= na);
    /* Fixup needed, now. Treat the largest possible subblock. We want
     * the intervals [0,na-k[ and [k, nb[ to match a kara scheme, i.e.
     * have the constraints 0<=k<=na, and nb-k = 2*(na-k)-1, i.e. k =
     * 2*na-1-nb.  By assumption, the expression above, with nb>=na,
     * guarantees that the first condition 0<=k<=na is satisfied.
     */
    unsigned int chop = 2*na-1-nb;
    polymat_mp_raw_kara(ab, c, xc,
            a, xa, na - chop, b, xb + chop, nb - chop, transpose, add);
    polymat_mp_raw_subdivide(ab, c, xc,
            a, xa + na - chop, chop, b, xb, na - 1, transpose, 1);
}

void polymat_mp_raw(abdst_field ab,
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb)
{
    polymat_mp_raw_subdivide(ab, c, xc, a, xa, na, b, xb, nb, 0, 0);
}
