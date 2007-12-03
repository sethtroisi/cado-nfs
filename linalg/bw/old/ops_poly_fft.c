#include <gmp.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "field_def.h"
#include "field_prime.h"
#include "field_quad.h"
#include "field_usage.h"
#include "fft_core.h"
#include "fft_core.h"
#include "ops_poly.h"
#include "auxfuncs.h"

#define	ft_coeff_stride	l_size

void ft_order(ft_order_t * s, int ncoeffs)
{
	*s = ceil_log2(ncoeffs);
}

void ft_zero(ft_order_t order, ft_t p)
{
    memset(p, 0, ft_coeff_stride * sizeof(mp_limb_t) << order);
}

void ft_one(ft_order_t order, ft_t p)
{
    int k;
    for (k = 0; k < (1 << order); k++) {
        l_set_one(p);
        p += ft_coeff_stride;
    }
}

/* This really is an ``addmul'' */
void ft_convolution(ft_order_t order, ft_t r, ft_t p, ft_t q)
{
    int i;
    for (i = 0; i < (1 << order); i++) {
	l_addmul(r, p, q);
	r += ft_coeff_stride;
	p += ft_coeff_stride;
	q += ft_coeff_stride;
    }
}

/* Computes the DFT of p(X)q(X)/X^n (where n = dg_kill)*/
void ft_convolution_special(ft_order_t order, ft_t r, ft_t p, ft_t q,
        unsigned int dg_kill)
{
    int k;
    mp_limb_t * tmp;
    tmp = malloc(l_size * sizeof(mp_limb_t));
    for (k = 0; k < (1 << order); k++) {
        int kr;
        kr = phi_tab[(phi_tab[k] *
                (-dg_kill)) & ((1 << tabulated_order) - 1)];
        l_mul(tmp, p, q);
        l_addmul(r, tmp, l_roots + kr * ft_coeff_stride);
        p += ft_coeff_stride;
        q += ft_coeff_stride;
        r += ft_coeff_stride;
    }
    free(tmp);
}

void ft_itransform(ft_order_t order, mp_limb_t * dst, ptrdiff_t stride, int deg,
        ft_t p)
{
    int k;
    mp_limb_t *one_over_n;

    one_over_n = malloc(k_size * sizeof(mp_limb_t));
    k_set_int(one_over_n, 1 << (order));
    k_inv(one_over_n, one_over_n);

    (*inverse_fft) (p, ft_coeff_stride, l_roots, order);
    for (k = 0; k <= deg; k++) {
        /*
         * The proper function to use is
         * restrict_scalars. But there's some more to do.
         */
#ifndef HAS_NATIVE_FFT
        assert(k_is_zero(p + k * ft_coeff_stride + k_size));
#endif
        k_mul(dst + k * stride, p + k * ft_coeff_stride, one_over_n);
    }
    /* Watch out: The other coefficients might very well be
     * non-zero. That's not incorrect, but these are
     * meaningless */
    free(one_over_n);
}

void ft_transform(ft_order_t order, ft_t p, mp_limb_t * src, ptrdiff_t stride, int deg)
{
    int k,l;
    for (k = 0; k < (1 << order); k++) {
        for (l = k; l <= deg; l += (1 << order)) {
            /* The proper function to use is extend_scalars, but it
             * doesn't do the addup thing... */
            k_addto(p + k * ft_coeff_stride, src + l * stride);
        }
    }
    (*direct_fft) (p, ft_coeff_stride, l_roots, order);
}

void ft_ops_init(unsigned int nmax, struct charptr_list_t * args)
{
	int fermat_prime = 0;
	int enable_cplx = 0;

	printf("Preparing FFT engine\n");

	for( ; args != NULL ; args = args->next) {
		if (strcmp(args->p, "--enable-complex-field") == 0) {
			enable_cplx = 1;
		} else if (strcmp(args->p, "--fermat-prime") == 0) {
			fermat_prime = 1;
		}
	}


	if (fermat_prime) {
#ifdef	HAS_NATIVE_FFT
		enable_cplx=-1;
#else
		die("To use the --fermat-prime option, "
			"recompile with -DHAS_NATIVE_FFT\n",1);
#endif
	}

	prepare_fft_engine(ceil_log2(nmax), enable_cplx);

}
void ft_ops_cleanup()
{
	cleanup_fft_engine();
}
