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
#include "auxfuncs.h"
#include "ops_poly_fft.hpp"

int ops_poly_fft::coeff_stride;

void ops_poly_fft::set(int ncoeffs)
{
	s = ceil_log2(ncoeffs);
}
void ops_poly_fft::zero(transform_t p) const
{
    memset(p, 0, coeff_stride * sizeof(mp_limb_t) << s);
}

void ops_poly_fft::one(transform_t p) const
{
    int k;
    for (k = 0; k < (1 << s); k++) {
        l_set_one(p);
        p += coeff_stride;
    }
}

/* This really is an ``addmul'' */
void ops_poly_fft::convolution(transform_t r, transform_t p, transform_t q) const
{
    int i;
    for (i = 0; i < 1 << s; i++) {
	l_addmul(r, p, q);
	r += coeff_stride;
	p += coeff_stride;
	q += coeff_stride;
    }
}

/* Computes the DFT of p(X)q(X)/X^n (where n = dg_kill)*/
void ops_poly_fft::convolution_special(transform_t r, transform_t p, transform_t q,
        unsigned int dg_kill) const
{
    int k;
    mp_limb_t * tmp;
    tmp = (mp_limb_t *) malloc(l_size * sizeof(mp_limb_t));
    for (k = 0; k < (1 << s); k++) {
        int kr;
        kr = phi_tab[(phi_tab[k] *
                (-dg_kill)) & ((1 << tabulated_order) - 1)];
        l_mul(tmp, p, q);
        l_addmul(r, tmp, l_roots + kr * coeff_stride);
        p += coeff_stride;
        q += coeff_stride;
        r += coeff_stride;
    }
    free(tmp);
}

void ops_poly_fft::itransform(mp_limb_t * dst, ptrdiff_t stride, int deg,
        transform_t p) const
{
    int k;
    mp_limb_t *one_over_n;

    one_over_n = (mp_limb_t *) malloc(k_size * sizeof(mp_limb_t));
    k_set_int(one_over_n, 1 << s);
    k_inv(one_over_n, one_over_n);

    (*inverse_fft) (p, coeff_stride, l_roots, s);
    for (k = 0; k <= deg; k++) {
        /*
         * The proper function to use is
         * restrict_scalars. But there's some more to do.
         */
#ifndef HAS_NATIVE_FFT
        assert(k_is_zero(p + k * coeff_stride + k_size));
#endif
        k_mul(dst + k * stride, p + k * coeff_stride, one_over_n);
    }
    /* Watch out: The other coefficients might very well be
     * non-zero. That's not incorrect, but these are
     * meaningless */
    free(one_over_n);
}

void ops_poly_fft::transform(transform_t p, mp_limb_t * src, ptrdiff_t stride, int deg) const
{
    int k,l;
    for (k = 0; k < (1 << s); k++) {
        for (l = k; l <= deg; l += (1 << s)) {
            /* The proper function to use is extend_scalars, but it
             * doesn't do the addup thing... */
            k_addto(p + k * coeff_stride, src + l * stride);
        }
    }
    (*direct_fft) (p, coeff_stride, l_roots, s);
}

void ops_poly_fft::init(unsigned int nmax, std::list<char *> const& args)
{
	int fermat_prime = 0;
	int enable_cplx = 0;
	std::list<char *>::const_iterator foo;

	printf("Preparing FFT engine\n");

	for(foo = args.begin() ; foo != args.end() ; foo++) {
		if (strcmp(*foo, "--enable-complex-field") == 0) {
			enable_cplx = 1;
		} else if (strcmp(*foo, "--fermat-prime") == 0) {
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

	coeff_stride = l_size;
}

void ops_poly_fft::cleanup()
{
	cleanup_fft_engine();
}

bool ops_poly_fft::operator==(ops_poly_fft const& b) const {
	return s == b.s;
}
bool ops_poly_fft::fits(int n) const {
	return n<=(1 << s);
}
ops_poly_fft::operator int() const {
	return s;
}
ops_poly_fft& ops_poly_fft::operator=(ops_poly_fft const& o) {
	s = o.s;
	return *this;
}
void * ops_poly_fft::mat_mb_alloc(mat_mb& x) const {
	return mbdft_alloc(s,x);
}
void ops_poly_fft::mat_mb_free(mat_mb x) const {
	mbdft_free(s,x);
}
ops_poly_fft::transform_t ops_poly_fft::mat_mb_get(mat_mb x, int i, int j) const {
	return mbdft_poly(s,x,i,j);
}
void * ops_poly_fft::mat_bb_alloc(mat_bb& x) const {
	return bbdft_alloc(s,x);
}
void ops_poly_fft::mat_bb_free(mat_bb x) const {
	bbdft_free(s,x);
}
ops_poly_fft::transform_t ops_poly_fft::mat_bb_get(mat_bb x, int i, int j) const {
	return bbdft_poly(s,x,i,j);
}
