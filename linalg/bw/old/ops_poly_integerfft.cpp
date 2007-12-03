#include <gmp.h>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>

#include "ops_poly_integerfft.hpp"
#include "bw_scalar.h"
#include "gmp-hacks.h"
#include "field_def.h"
#include "field_common.h"
#include "field_prime.h"
#include "field_usage.h"
#include "variables.h"

using namespace std;

/* This file requires a patched version of GMP, using the patch at the
 * following URL:
 *
 * http://www.loria.fr/~zimmerma/bignum/fft_patch.gmp-4.2.1
 *
 */
extern "C" {
extern mp_size_t mpn_fft_init (mp_fft_t, mp_size_t, unsigned long);
extern void mpn_fft_transform (mp_fft_t, mp_ptr, mp_size_t);
extern void mpn_fft_mul (mp_fft_t, mp_fft_t);
extern void mpn_fft_add (mp_fft_t, mp_fft_t);
extern void mpn_fft_sub (mp_fft_t, mp_fft_t);
extern void mpn_fft_itransform (mp_ptr, mp_fft_t);
extern void mpn_fft_clear (mp_fft_t);
extern int __gmpn_fft_best_k (mp_size_t n, int sqr);
extern mp_size_t __gmpn_fft_next_size (mp_size_t pl, int k);
void mpn_fft_copy (mp_fft_t r, mp_fft_t s)
{
   MPN_COPY (r->d[0], s->d[0], (s->nprime + 1) << s->k);
}
void mpn_fft_zero (mp_fft_t s)
{
   MPN_ZERO (s->d[0], (s->nprime + 1) << s->k);
}
void mpn_fft_one (mp_fft_t s)
{
	mpn_fft_zero(s);
	for(int i = 0 ; i < (1 << s->k) ; i++) {
		s->d[i][0] = 1UL;
	}
}
}

// XXX debug.
mpz_t bound;

mp_size_t ops_poly_ifft::limbs_per_coeff;

void ops_poly_ifft::set(int nc)
{
    // unsigned long k, K, M, Nprime, nprime, maxLK;
    // int n;
    ncoeffs = limbs_per_coeff * nc;
    /* These are used only for comparison between fft blocks. They have
     * to match the parameters chosen by fft_init, but fortunately it is
     * not too messy.
     */
    k = __gmpn_fft_best_k (ncoeffs, 0);
    n = __gmpn_fft_next_size (ncoeffs, k);
    // K = 1 << k;
    // M = (n * GMP_NUMB_BITS) / K;
    // maxLK = (K > GMP_NUMB_BITS) ? K : GMP_NUMB_BITS;
    // Nprime = ((2 * M + k + 2 + maxLK) / maxLK) * maxLK;
    // nprime = Nprime / GMP_NUMB_BITS;
}

static void _ft_mat_alloc_init_many(ops_poly_ifft::transform_t & x, int nb, int ncoeffs)
{
	unsigned long l;
	x = (mp_fft_t *) malloc(nb * sizeof(mp_fft_t));
	if (x == NULL) {
		return;
	}
	for(l = 1 ; bigdim > (1 << l) ; l++);
	for(int i = 0 ; i < nb ; i++) {
		mpn_fft_init(x[i], ncoeffs, l);
	}
}

static void _ft_mat_clear_free_many(ops_poly_ifft::transform_t p, int nb)
{
	for(int i = 0 ; i < nb ; i++) {
		mpn_fft_clear(p[i]);
	}
	free(p);
}


void ops_poly_ifft::zero(transform_t p) const { mpn_fft_zero(p[0]); }

void ops_poly_ifft::one(transform_t p) const { mpn_fft_one(p[0]); }

void ops_poly_ifft::convolution(transform_t r, transform_t p, transform_t q) const
{
	/* ugly ; find another way, perhaps directly in the order field
	 */
	mp_fft_t tmp;
	unsigned long l;
	for(l = 1 ; bigdim > (1 << l) ; l++);
	mpn_fft_init(tmp, ncoeffs, bigdim);
	mpn_fft_copy(tmp, p[0]);
	mpn_fft_mul(tmp, q[0]);
	mpn_fft_add(r[0], tmp);
	mpn_fft_clear(tmp);
}

#if 0
void ops_poly_ifft::convolution_special(ft_order_t order,
		transform_t r, transform_t p, transform_t q,
		unsigned int dg_kill)
{
	BUG_ON(dg_kill);
	/* probably doable, not necessarily trivial though */
	ft_convolution(order, r, p, q);
}
#endif

void ops_poly_ifft::itransform(
		mp_limb_t * dst, ptrdiff_t stride, int deg,
		transform_t p) const
{
	mp_limb_t * tmp;
	int i;
	assert((deg+1) <= ncoeffs);
	tmp = (mp_limb_t *) malloc((p[0]->n+1) * sizeof(mp_limb_t));
	memset(tmp, 0, (p[0]->n+1) * sizeof(mp_limb_t));
	mpn_fft_itransform(tmp, p[0]);
	for(i = 0 ; i <= deg ; i++) {
		if (mpn_cmp(tmp + i * limbs_per_coeff, PTR(bound), limbs_per_coeff) >= 0) {
			fprintf(stderr, "Argh, ifft bound failed\n");
			abort();
		}
		k_set_mpn(dst + i * stride, 
				tmp + i * limbs_per_coeff,
				limbs_per_coeff);
	}
	free(tmp);
}

void ops_poly_ifft::transform(
		transform_t p,
		mp_limb_t * src, ptrdiff_t stride, int deg) const
{
	mp_limb_t * tmp;
	int i;
	mp_size_t nlimbs = (deg + 1) * limbs_per_coeff;
	assert((deg+1) <= ncoeffs);
	tmp = (mp_limb_t *) malloc(nlimbs * sizeof(mp_limb_t));
	memset(tmp, 0, nlimbs * sizeof(mp_limb_t));
	for(i = 0 ; i <= deg ; i++) {
		k_reduce(src + i * stride);
		memcpy(tmp + i * limbs_per_coeff, src + i * stride, 
				bw_allocsize * sizeof(mp_limb_t));
	}
	mpn_fft_transform(p[0], tmp, nlimbs);
	free(tmp);
}

void ops_poly_ifft::init(unsigned int nmax, std::list<char *> const& args)
{
	mpz_t z;
	mp_size_t nbits, nlimbs;

	field_k=new_prime_field(modulus_plain,bw_allocsize);
	
	mpz_init(z);
	mpz_set(z, modulus);
	mpz_mul(z, z, z);
	mpz_mul_ui(z, z, nmax);
	mpz_mul_ui(z, z, bigdim);

	nbits = mpz_sizeinbase(z, 2);
	nlimbs = SIZ(z);

	mpz_init_set(bound, z);

	limbs_per_coeff = nlimbs;

	mpz_clear(z);

	printf("// Using %ld limbs per coefficient\n", limbs_per_coeff);
}
void ops_poly_ifft::cleanup()
{
	mpz_clear(bound);
}

bool ops_poly_ifft::operator==(ops_poly_ifft const& b) const {
	return (n==b.n) && (k == b.k);
}
bool ops_poly_ifft::fits(int nc) const {
	return (nc * limbs_per_coeff) <= n;
}
ops_poly_ifft::operator int() const {
	return k;
}
ops_poly_ifft& ops_poly_ifft::operator=(ops_poly_ifft const& o) {
	memcpy(this, &o, sizeof(self));
	return *this;
}
void * ops_poly_ifft::mat_mb_alloc(mat_mb& x) const {
	_ft_mat_alloc_init_many(x, m_param*bigdim, ncoeffs);
	return x;
}
void ops_poly_ifft::mat_mb_free(mat_mb x) const {
	_ft_mat_clear_free_many(x,m_param*bigdim);
}
ops_poly_ifft::transform_t ops_poly_ifft::mat_mb_get(mat_mb x, int i, int j) const {
	return ((x) + ((i) * bigdim + (j)));
}
void * ops_poly_ifft::mat_bb_alloc(mat_bb& x) const {
	_ft_mat_alloc_init_many(x, bigdim*bigdim, ncoeffs);
	return x;
}
void ops_poly_ifft::mat_bb_free(mat_bb x) const {
	_ft_mat_clear_free_many(x, bigdim*bigdim);
}
ops_poly_ifft::transform_t ops_poly_ifft::mat_bb_get(mat_bb x, int i, int j) const {
	return ((x) + ((i) * bigdim + (j)));
}
