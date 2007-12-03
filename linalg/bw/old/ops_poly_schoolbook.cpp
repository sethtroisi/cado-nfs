#include <string.h>
#include "field_def.h"
#include "field_prime.h"
#include "field_quad.h"
#include "field_usage.h"
#include "fft_core.h"
#include "ops_poly_schoolbook.hpp"
#include "auxfuncs.h"
#include "bw_scalar.h"
#include "variables.h"

/* This is a simplistic approach where the ``transform'' is actually the
 * identity ; all the work is done in the convolution steps using
 * schoolbook multiplication. Therefore this code is SLOW, and only of
 * expository value.
 */

int ops_poly_scbk::coeff_stride;
void ops_poly_scbk::set(int ncoeffs)
{
	/* A security margin seems wise. This one is really a rule of
	 * thumb for one small example.
	 */
	// *s = ncoeffs * 1.2;
	
	/* And this one is here so that we fail if and only if the _fft
	 * version fails */
	s = 1 << ceil_log2(ncoeffs);
}

void ops_poly_scbk::zero(transform_t p) const
{                       
    memset(p, 0, coeff_stride * s * sizeof(mp_limb_t));
}       

void ops_poly_scbk::one(transform_t p) const
{
	zero(p);
	k_set_one(p);
}

/* This really is an ``addmul'' */
void ops_poly_scbk::convolution(transform_t r, transform_t p, transform_t q) const
{
    int i,k;
    for (k = 0; k < s; k++) {
        for (i = 0; i < s; i++) {
            k_addmul(r + i * coeff_stride,
			    p + k * coeff_stride,
			    q + ((i+s-k) % s) * coeff_stride);
        }
    }
}

/* Computes the DFT of p(X)q(X)/X^n */
/* This really is an ``addmul'' */
void ops_poly_scbk::convolution_special(transform_t r,
		transform_t p, transform_t q,
		unsigned int dg_kill) const
{
    int i, k;
    for (i = 0 ; i < s; i++) {
        for(k = 0 ; k < s ; k++) {
            k_addmul(r + i * coeff_stride,
                    p + k * coeff_stride,
                    q + ((dg_kill + s + i - k) % s) * coeff_stride);
        }
    }
}


void ops_poly_scbk::itransform(mp_limb_t * dst, ptrdiff_t stride, int deg,
        transform_t p) const
{
    int k;
    for (k = 0; k <= deg; k++) {
        k_set(dst, p);
        dst += stride;
        p += coeff_stride;
    }
}

void ops_poly_scbk::transform(transform_t p,
		mp_limb_t * src, ptrdiff_t stride, int deg) const
{
    int k, l;
    for (k = 0; k < s; k++) {
        for(l = k ; l <= deg ; l += s) {
            k_addto(p + k * coeff_stride, src + l * stride);
        }
    }
}

void ops_poly_scbk::init(unsigned int nmax, std::list<char *> const& args)
{
	field_k=new_prime_field(modulus_plain,bw_allocsize);
	coeff_stride = bw_allocsize;
}

void ops_poly_scbk::cleanup()
{
}

bool ops_poly_scbk::operator==(ops_poly_scbk const& b) const {
        return s == b.s;
}
bool ops_poly_scbk::fits(int n) const {
        return n<=s;
}
ops_poly_scbk::operator int() const {
        return s;
}
ops_poly_scbk& ops_poly_scbk::operator=(ops_poly_scbk const& o) {
        s = o.s;
        return *this;
}

void * ops_poly_scbk::mat_mb_alloc(mat_mb& x) const {
        return x = (mat_mb) mymalloc(m_param * bigdim * (coeff_stride * s) * sizeof(mp_limb_t));
}
void ops_poly_scbk::mat_mb_free(mat_mb x) const {
	free(x);
}
ops_poly_scbk::transform_t ops_poly_scbk::mat_mb_get(mat_mb x, int i, int j) const {
	return ((x) + ((i) * bigdim + (j)) * (coeff_stride * s));
}

void * ops_poly_scbk::mat_bb_alloc(mat_bb& x) const {
        return x = (mat_bb) mymalloc(bigdim * bigdim * (coeff_stride * s) * sizeof(mp_limb_t));
}
void ops_poly_scbk::mat_bb_free(mat_bb x) const {
	free(x);
}
ops_poly_scbk::transform_t ops_poly_scbk::mat_bb_get(mat_bb x, int i, int j) const {
	return ((x) + ((i) * bigdim + (j)) * (coeff_stride * s));
}


