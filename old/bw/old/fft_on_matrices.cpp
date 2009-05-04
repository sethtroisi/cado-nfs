#include <sys/time.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <gmp.h>
#include "lingen_params.h"
#include "params.h"
#include "types.h"
#include "macros.h"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "variables.h"
#include "structure.h"
#include "gmp-hacks.h"
#include "modulus_hacks.h"
#include "twisting_polynomials.h"
#include "e_polynomial.h"
#include "ops_poly.hpp"
#include "fft_on_matrices.hpp"
#include "field_def.h"
#include "field_prime.h"
#include "field_quad.h"
#include "field_usage.h"
#include "timer.h"
#include "manu.h"

/**********************************************************************/
/**********************************************************************/


/**********************************************************************/

/* XXX : locality is certainly important, unless multiplications take a
 * huge time anyway. So it might be better to keep polynomials gathered
 * in a single area.
 */
struct dft_mb *fft_ec_dft(struct e_coeff *ec, ft_order_t const& order, double *tm)
{
    struct dft_mb *res;
    int i, j;
    struct timeval tv;

    timer_r(&tv, TIMER_SET);

    res = (dft_mb *) malloc(sizeof(struct dft_mb));
    if (res == NULL)
	return res;
    res->degree = ec->degree;
    res->order = order;
    if (!(order.mat_mb_alloc(res->p))) {
	free(res);
	return NULL;
    }

    /*
     * We want to compute the order n DFT, n==res->order. ec->degree
     * might actually be bigger than 1<<n, so we reduce e(X)
     * transparently modulo X^{1<<n} in order to get the proper DFT
     */
    for (i = 0; i < m_param; i++) {
	for (j = 0; j < bigdim; j++) {
            mp_limb_t * src0, * src1;
            src0 = mbmat_scal(mbpoly_coeff(ec->p, 0), i, j);
            src1 = mbmat_scal(mbpoly_coeff(ec->p, 1), i, j);
            order.zero(order.mat_mb_get(res->p, i, j));
            order.transform(order.mat_mb_get(res->p, i, j),
                    src0, src1 - src0, ec->degree);
	    printf("."); fflush(stdout);
	}
    }
    printf("\n");

    *tm = timer_r(&tv, TIMER_ASK);

    return res;
}


struct dft_bb *fft_tp_dft(struct t_poly *tp, ft_order_t const& order, double *tm)
{
    struct dft_bb *res;
    int i, j, jr;
    struct timeval tv;

    timer_r(&tv, TIMER_SET);

    res = (dft_bb *) malloc(sizeof(struct dft_bb));
    if (res == NULL)
	return res;
    res->degree = tp->degree;
    res->order = order;
    if (!(order.mat_bb_alloc(res->p))) {
	free(res);
	return NULL;
    }

    /* Same thing as above */
    for (i = 0; i < bigdim; i++) {
	for (j = 0; j < bigdim; j++) {
            mp_limb_t * src0, * src1;

	    jr = tp->clist[j];

            src0 = bbmat_scal(bbpoly_coeff(tp->p, 0), i, jr);
            src1 = bbmat_scal(bbpoly_coeff(tp->p, 1), i, jr);

            order.zero(order.mat_bb_get(res->p, i, j));
            order.transform(order.mat_bb_get(res->p, i, j),
                    src0, src1-src0, tp->degnom[jr]);

	    printf("."); fflush(stdout);
	}
    }
    printf("\n");

    *tm = timer_r(&tv, TIMER_ASK);

    return res;
}

#ifdef  HAS_CONVOLUTION_SPECIAL
struct dft_mb *fft_mbb_conv_sp(struct dft_mb *p,
			       struct dft_bb *q,
			       unsigned int dg_kill, double *tm)
{
    struct dft_mb *res;
    int i, j, l;
    struct timeval tv;

    timer_r(&tv, TIMER_SET);

    res = (dft_mb *) malloc(sizeof(struct dft_mb));
    if (res == NULL)
	return res;
    res->degree = p->degree + q->degree - dg_kill;
    ASSERT(p->order == q->order);;
    res->order = p->order;
    if (!(res->order.mat_mb_alloc(res->p))) {
	free(res);
	return NULL;
    }

    /* TODO: Multiply using Strassen's algorithm instead. */
    for (i = 0; i < m_param; i++) {
	for (j = 0; j < bigdim; j++) {
            res->order.zero(res->order.mat_mb_get(res->p, i, j));
	    for (l = 0; l < bigdim; l++) {
                res->order.convolution_special(
                        res->order.mat_mb_get(res->p, i, j),
                        res->order.mat_mb_get(p->p, i, l),
                        res->order.mat_bb_get(q->p, l, j),
                        dg_kill);
	    }
	    printf("."); fflush(stdout);
	}
    }
    printf("\n");

    *tm = timer_r(&tv, TIMER_ASK);

    return res;
}
#endif

/* This one is the same, without the feature that we divide by X^d
 * directly in the dft. That way, the dft is longer. This function is
 * only relevant for testing.
 */
struct dft_mb *fft_mbb_conv(struct dft_mb *p,
			       struct dft_bb *q,
			       double *tm)
{
    struct dft_mb *res;
    int i, j, l;
    struct timeval tv;

    timer_r(&tv, TIMER_SET);

    res = (dft_mb *) malloc(sizeof(struct dft_mb));
    if (res == NULL)
	return res;
    res->degree = p->degree + q->degree;
    ASSERT(p->order == q->order);;
    res->order = p->order;
    if (!(res->order.mat_mb_alloc(res->p))) {
	free(res);
	return NULL;
    }

    /* TODO: Multiply using Strassen's algorithm instead. */
    for (i = 0; i < m_param; i++) {
	for (j = 0; j < bigdim; j++) {
            res->order.zero(res->order.mat_mb_get(res->p, i, j));
	    for (l = 0; l < bigdim; l++) {
                res->order.convolution(
                        res->order.mat_mb_get(res->p, i, j),
                        res->order.mat_mb_get(p->p, i, l),
                        res->order.mat_bb_get(q->p, l, j));
	    }
	    printf("."); fflush(stdout);
	}
    }
    printf("\n");

    *tm = timer_r(&tv, TIMER_ASK);

    return res;
}

struct dft_bb *fft_bbb_conv(struct dft_bb *p, struct dft_bb *q, double *tm)
{
    struct dft_bb *res;
    int i, j, l;
    struct timeval tv;

    timer_r(&tv, TIMER_SET);

    res = (dft_bb *) malloc(sizeof(struct dft_bb));
    if (res == NULL)
	return res;
    res->degree = 0;		/* we don't care, this is meaningless anyway */
    ASSERT(p->order == q->order);;
    res->order = p->order;
    if (!(res->order.mat_bb_alloc(res->p))) {
	free(res);
	return NULL;
    }

    /* TODO: Multiply using Strassen's algorithm instead. */
    for (i = 0; i < bigdim; i++) {
	for (j = 0; j < bigdim; j++) {
	    res->order.zero(res->order.mat_bb_get(res->p, i, j));
	    for (l = 0; l < bigdim; l++) {
		res->order.convolution(
			       res->order.mat_bb_get(res->p, i, j),
			       res->order.mat_bb_get(p->p, i, l),
			       res->order.mat_bb_get(q->p, l, j));
	    }
	    printf("."); fflush(stdout);
	}
    }
    printf("\n");

    *tm = timer_r(&tv, TIMER_ASK);

    return res;
}

void fft_mb_invdft(bw_mbpoly dest,
		   struct dft_mb *p, unsigned int deg, double *tm)
{
    int i, j;
    struct timeval tv;

    timer_r(&tv, TIMER_SET);

    ASSERT(p->order.fits(deg + 1));

    for (i = 0; i < m_param; i++) {
	for (j = 0; j < bigdim; j++) {
            mp_limb_t * dst0;
            mp_limb_t * dst1;
            dst0 = mbmat_scal(mbpoly_coeff(dest, 0), i, j);
            dst1 = mbmat_scal(mbpoly_coeff(dest, 1), i, j);
            p->order.itransform(
                    dst0, dst1-dst0, deg, p->order.mat_mb_get(p->p, i, j));

	    printf("."); fflush(stdout);
	}
    }
    printf("\n");

    *tm = timer_r(&tv, TIMER_ASK);
}

void fft_tp_invdft(struct t_poly *tp, struct dft_bb *p, double *tm)
{
    int i, j;
    struct timeval tv;

    timer_r(&tv, TIMER_SET);

    /*
     * Les degrés de destinations sont déjà dans tp.
     *
     * Les colonnes ont déjà été réordonnées par tp_comp_alloc, donc
     * la permutation associée à tp est triviale : clist[j]==j.
     *
     * Par ailleurs, tp_comp_alloc a déjà tout mis à zéro.
     */

    for (i = 0; i < bigdim; i++) {
	for (j = 0; j < bigdim; j++) {
            mp_limb_t * dst0;
            mp_limb_t * dst1;
            dst0 = bbmat_scal(bbpoly_coeff(tp->p, 0), i, j);
            dst1 = bbmat_scal(bbpoly_coeff(tp->p, 1), i, j);
            p->order.itransform(dst0, dst1-dst0, tp->degnom[j],
                    p->order.mat_bb_get(p->p, i, j));
	    printf("."); fflush(stdout);
	}
    }
    printf("\n");

    *tm = timer_r(&tv, TIMER_ASK);
}

struct dft_bb *fft_bb_dft_init_one(unsigned int deg)
{
    struct dft_bb *res;
    int i, j;

    res = (dft_bb *) malloc(sizeof(struct dft_bb));
    if (res == NULL)
	return res;
    res->degree = 0;
    res->order.set(deg + 1);
    if (!(res->order.mat_bb_alloc(res->p))) {
	free(res);
	return NULL;
    }

    for (i = 0; i < bigdim; i++) {
	for (j = 0; j < bigdim; j++) {
            res->order.zero(res->order.mat_bb_get(res->p, i, j));
        }
        res->order.one(res->order.mat_bb_get(res->p, i, i));
    }

    return res;
}

void dft_mb_free(struct dft_mb *p)
{
    p->order.mat_mb_free(p->p);
    free(p);
}

void dft_bb_free(struct dft_bb *p)
{
    p->order.mat_bb_free(p->p);
    free(p);
}

/* vim:set sw=4 sta et: */
