/* mul_fft -- radix 2 fft routines for MPIR.
 * 
 * Copyright 2009, 2011 William Hart. All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY William Hart ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN 
 * NO EVENT SHALL William Hart OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are
 * those of the authors and should not be interpreted as representing
 * official policies, either expressed or implied, of William Hart.
 * 
 */

/******************************************************************************

    Copyright (C) 2009, 2011 William Hart
 
******************************************************************************/

#ifndef FFT_H
#define FFT_H

#undef ulong
#define ulong ulongxx		/* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t
#include "flint.h"
#include "mpn_extras.h"

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__MPIR_VERSION)

#if !defined(__MPIR_RELEASE ) || __MPIR_RELEASE < 20600
#define mpn_sumdiff_n __MPN(sumdiff_n)
extern
mp_limb_t mpn_sumdiff_n(mp_ptr, mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
#endif

#else

static __inline__ mp_limb_t
mpn_sumdiff_n(mp_ptr s, mp_ptr d, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    mp_limb_t ret;
    mp_ptr t;

    if (n == 0)
	return 0;

    if ((s == x && d == y) || (s == y && d == x)) {
	t = flint_malloc(n * sizeof(mp_limb_t));
	ret = mpn_sub_n(t, x, y, n);
	ret += 2 * mpn_add_n(s, x, y, n);
	flint_mpn_copyi(d, t, n);
	flint_free(t);
	return ret;
    }

    if (s == x || s == y) {
	ret = mpn_sub_n(d, x, y, n);
	ret += 2 * mpn_add_n(s, x, y, n);
	return ret;
    }

    ret = 2 * mpn_add_n(s, x, y, n);
    ret += mpn_sub_n(d, x, y, n);
    return ret;
}

#endif

#define fft_sumdiff(t, u, r, s, n) \
   (n == 0 ? 0 : mpn_sumdiff_n(t, u, r, s, n))


#define SWAP_PTRS(xx, yy) \
   do { \
      mp_limb_t * __ptr = xx; \
      xx = yy; \
      yy = __ptr; \
   } while (0)

/* used for generating random values mod p in test code */
#define random_fermat(nn, state, limbs) \
   do { \
      if (n_randint(state, 10) == 0) { \
         flint_mpn_zero(nn, limbs); \
         nn[limbs] = 1; \
      } else { \
         if (n_randint(state, 2) == 0) \
            flint_mpn_rrandom(nn, state->gmp_state, limbs); \
         else \
            flint_mpn_urandomb(nn, state->gmp_state, limbs*FLINT_BITS); \
         nn[limbs] = n_randint(state, 1024); \
      } \
      if (n_randint(state, 2)) \
         nn[limbs] = -nn[limbs]; \
   } while (0)

static __inline__
/*
 * Adds the signed limb \code{c} to the generalised fermat number \code{r}
 * modulo \code{B^limbs + 1}. The compiler should be able to inline
 * this for the case that there is no overflow from the first limb.
 * 
 */
    void mpn_addmod_2expp1_1(mp_limb_t * r, mp_size_t limbs,
			     mp_limb_signed_t c)
{
    mp_limb_t sum = r[0] + c;

    /* check if adding c would cause a carry to propagate */
    if ((mp_limb_signed_t) (sum ^ r[0]) >= 0)
	r[0] = sum;
    else {
	if (c >= 0)
	    mpn_add_1(r, r, limbs + 1, c);
	else
	    mpn_sub_1(r, r, limbs + 1, -c);
    }
}

void fft_combine_limbs(mp_limb_t * res, mp_limb_t ** poly, slong length,
		       mp_size_t coeff_limbs, mp_size_t output_limbs,
		       mp_size_t total_limbs);

void fft_combine_bits(mp_limb_t * res, mp_limb_t ** poly, slong length,
		      mp_bitcnt_t bits, mp_size_t output_limbs,
		      mp_size_t total_limbs);

mp_size_t fft_split_limbs(mp_limb_t ** poly, mp_srcptr limbs,
			  mp_size_t total_limbs, mp_size_t coeff_limbs,
			  mp_size_t output_limbs);

mp_size_t fft_split_bits(mp_limb_t ** poly, mp_srcptr limbs,
			 mp_size_t total_limbs, mp_bitcnt_t bits,
			 mp_size_t output_limbs);

void fermat_to_mpz(mpz_t m, mp_limb_t * i, mp_size_t limbs);

void mpn_normmod_2expp1(mp_limb_t * t, mp_size_t limbs);

void butterfly_lshB(mp_limb_t * t, mp_limb_t * u, mp_limb_t * i1,
		    mp_limb_t * i2, mp_size_t limbs, mp_size_t x,
		    mp_size_t y);

void butterfly_rshB(mp_limb_t * t, mp_limb_t * u, mp_limb_t * i1,
		    mp_limb_t * i2, mp_size_t limbs, mp_size_t x,
		    mp_size_t y);

void mpn_mul_2expmod_2expp1(mp_limb_t * t,
			    mp_limb_t * i1, mp_size_t limbs, mp_bitcnt_t d);

void mpn_div_2expmod_2expp1(mp_limb_t * t,
			    mp_limb_t * i1, mp_size_t limbs, mp_bitcnt_t d);

void fft_adjust(mp_limb_t * r, mp_limb_t * i1,
		mp_size_t i, mp_size_t limbs, mp_bitcnt_t w);

void fft_butterfly(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1,
		   mp_limb_t * i2, mp_size_t i, mp_size_t limbs,
		   mp_bitcnt_t w);

void ifft_butterfly(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1,
		    mp_limb_t * i2, mp_size_t i, mp_size_t limbs,
		    mp_bitcnt_t w);

void fft_radix2(mp_limb_t ** ii,
		mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2);

void fft_truncate1(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w,
		   mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc);

void fft_truncate(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w,
		  mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc);

void ifft_radix2(mp_limb_t ** ii, mp_size_t n,
		 mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2);

void ifft_truncate1(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w,
		    mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc);

void ifft_truncate(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w,
		   mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc);

void fft_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t,
			 mp_limb_t * i1, mp_limb_t * i2, mp_size_t i,
			 mp_size_t limbs, mp_bitcnt_t w, mp_limb_t * temp);

void ifft_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1,
			  mp_limb_t * i2, mp_size_t i, mp_size_t limbs,
			  mp_bitcnt_t w, mp_limb_t * temp);

void fft_adjust_sqrt2(mp_limb_t * r, mp_limb_t * i1,
		      mp_size_t i, mp_size_t limbs, mp_bitcnt_t w,
		      mp_limb_t * temp);

void fft_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w,
			mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
			mp_size_t trunc);

void ifft_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w,
			 mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp,
			 mp_size_t trunc);

void mul_truncate_sqrt2(mp_ptr r1, mp_srcptr i1, mp_size_t n1,
			mp_srcptr i2, mp_size_t n2, mp_bitcnt_t depth,
			mp_bitcnt_t w);

void fft_butterfly_twiddle(mp_limb_t * u, mp_limb_t * v,
			   mp_limb_t * s, mp_limb_t * t, mp_size_t limbs,
			   mp_bitcnt_t b1, mp_bitcnt_t b2);

void ifft_butterfly_twiddle(mp_limb_t * u, mp_limb_t * v,
			    mp_limb_t * s, mp_limb_t * t, mp_size_t limbs,
			    mp_bitcnt_t b1, mp_bitcnt_t b2);

void fft_radix2_twiddle(mp_limb_t ** ii, mp_size_t is,
			mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1,
			mp_limb_t ** t2, mp_size_t ws, mp_size_t r,
			mp_size_t c, mp_size_t rs);

void ifft_radix2_twiddle(mp_limb_t ** ii, mp_size_t is,
			 mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1,
			 mp_limb_t ** t2, mp_size_t ws, mp_size_t r,
			 mp_size_t c, mp_size_t rs);

void fft_truncate1_twiddle(mp_limb_t ** ii, mp_size_t is,
			   mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1,
			   mp_limb_t ** t2, mp_size_t ws, mp_size_t r,
			   mp_size_t c, mp_size_t rs, mp_size_t trunc);

void ifft_truncate1_twiddle(mp_limb_t ** ii, mp_size_t is,
			    mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1,
			    mp_limb_t ** t2, mp_size_t ws, mp_size_t r,
			    mp_size_t c, mp_size_t rs, mp_size_t trunc);

void fft_mfa_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n,
			    mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
			    mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc);

void ifft_mfa_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n,
			     mp_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2,
			     mp_limb_t ** temp, mp_size_t n1,
			     mp_size_t trunc);

void mul_mfa_truncate_sqrt2(mp_ptr r1, mp_srcptr i1, mp_size_t n1,
			    mp_srcptr i2, mp_size_t n2, mp_bitcnt_t depth,
			    mp_bitcnt_t w);

void fft_mfa_truncate_sqrt2_outer(mp_limb_t ** ii, mp_size_t n,
				  mp_bitcnt_t w, mp_limb_t ** t1,
				  mp_limb_t ** t2, mp_limb_t ** temp,
				  mp_size_t n1, mp_size_t trunc);

void fft_mfa_truncate_sqrt2_inner(mp_limb_t ** ii, mp_limb_t ** jj,
				  mp_size_t n, mp_bitcnt_t w, mp_limb_t ** t1,
				  mp_limb_t ** t2, mp_limb_t ** temp,
				  mp_size_t n1, mp_size_t trunc,
				  mp_limb_t * tt);

void ifft_mfa_truncate_sqrt2_outer(mp_limb_t ** ii, mp_size_t n,
				   mp_bitcnt_t w, mp_limb_t ** t1,
				   mp_limb_t ** t2, mp_limb_t ** temp,
				   mp_size_t n1, mp_size_t trunc);

void fft_negacyclic(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w,
		    mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp);

void ifft_negacyclic(mp_limb_t ** ii, mp_size_t n, mp_bitcnt_t w,
		     mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp);

void fft_naive_convolution_1(mp_limb_t * r, mp_limb_t * ii,
			     mp_limb_t * jj, mp_size_t m);

void _fft_mulmod_2expp1(mp_limb_t * r1, mp_limb_t * i1, mp_limb_t * i2,
			mp_size_t r_limbs, mp_bitcnt_t depth, mp_bitcnt_t w);

slong fft_adjust_limbs(mp_size_t limbs);

void fft_mulmod_2expp1(mp_limb_t * r, mp_limb_t * i1, mp_limb_t * i2,
		       mp_size_t n, mp_size_t w, mp_limb_t * tt);

void flint_mpn_mul_fft_main(mp_ptr r1, mp_srcptr i1, mp_size_t n1,
			    mp_srcptr i2, mp_size_t n2);

void fft_convolution(mp_limb_t ** ii, mp_limb_t ** jj, slong depth,
		     slong limbs, slong trunc, mp_limb_t ** t1,
		     mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t * tt);



struct fft_transform_info {
    mp_bitcnt_t bits1;
    mp_bitcnt_t bits2;
    unsigned int nacc;
    mp_size_t w;        /* Use sqrt(2)^w as a root of unity */
    mp_size_t depth;    /* Let n=2^depth. We work modulo 2^(wn)+1. Do a
                           transform length 4n. */
    mp_bitcnt_t bits;   /* Chunk sizes in bits */
    mp_size_t trunc0;   /* Number of coeffs of the transform computed.
                           This is not exactly the fourier transform
                           truncation point, because the truncation point
                           is also subject to a few extra criteria. */
    mp_bitcnt_t ks_coeff_bits;  /* This is used only for kronecker substitution */
    mp_bitcnt_t minwrap;        /* zero when no wraparound wanted */
    int alg;            /* alg==1: use matrix fourier algorithm */
};

/* The transform data is provided as an opaque void* pointer ; in reality
 * its structure may be written as the following pseudo-C code.
 *  struct transform_data {
 *      union {
 *          mp_limb_t * x[4 << depth];
 *          mp_limb_t * xy[2 << depth2][1 << depth1];
 *          // 1 + depth = depth1 + depth2
 *      } p;
 *      mp_limb_t * extra[2];
 *      mp_limb_t data[0];
 *  };
 * with pointers in x[], or xy[][], or extra[], all pointing to areas
 * within the zone beginning at data[]. Those areas are coefficients in
 * R=Z/(2^(nw)+1), each occupying fti_rsize0(fti)+1 limbs. All these
 * coefficient areas are disjoint. xy[] and x[] are of course to ways of
 * accessing the same set of pointers.

 * A transform data object represents a sequence of 4n coefficients in R
 * (pointed to by x[]). Before a DFT, this may be intepreted as a
 * polynomial P modulo x^(4n)-1. After, it's a sequence of pointwise
 * evaluations of P.  The layout is to be interpreted as follows.

 * with fti->alg == 0:
 *     x[i] == P( sqrt(2)^(bitrev(i, depth+2)) )
 * with fti->alg == 1 (matrix algorithm), The meaning of xy[] differs in
 * the two halves of the array.  We use the notation n1 == 1<<depth1,
 * n2==1<<depth2. The contents of xy[][] are given as follows.  Let
 * i<n2. 
 
 * xy[bitrev(i,depth2)][j]             == P(sqrt(2)^(    2*(j + i<<depth1))
 *                                        P(             2^(j + i<<depth1))
 * xy[bitrev(i,depth2) + 1<<depth2][j] == P(sqrt(2)^(1 + 2*(j + i<<depth1)))
 *                                        P(   sqrt(2) * 2^(j + i<<depth1))

 * In case of truncation, only some of the entries are considered. For
 * fti->alg==0, these are the entries with i < trunc. For fti->alg==1, the
 * rule is instead (i+b<<depth2) < trunc, with b being the most significant
 * bit of the row index.
 */

void fft_get_transform_info(struct fft_transform_info * fti, mp_bitcnt_t bits1, mp_bitcnt_t bits2, unsigned int nacc);
void fft_get_transform_info_mulmod(struct fft_transform_info * fti, mp_bitcnt_t xbits, mp_bitcnt_t ybits, unsigned int nacc, mp_bitcnt_t minwrap);
void fft_get_transform_allocs(size_t sizes[3], struct fft_transform_info * fti);
void fft_transform_prepare(void * x, struct fft_transform_info * fti);
void fft_do_dft(void * y, mp_limb_t * x, mp_size_t nx, void * temp, struct fft_transform_info * fti);
void fft_do_ift(mp_limb_t * x, mp_size_t nx, void * y, void * temp, struct fft_transform_info * fti);
void fft_mul(void * z, void * y0, void * y1, void * temp, struct fft_transform_info * fti);
void fft_add(void * z, void * y0, void * y1, struct fft_transform_info * fti);
void fft_get_transform_info_fppol(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t n1, mp_size_t n2, unsigned int nacc);
void fft_get_transform_info_fppol_mp(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t nmin, mp_size_t nmax, unsigned int nacc);
void fft_do_dft_fppol(void * y, mp_limb_t * x, mp_size_t nx, void * temp, struct fft_transform_info * fti, mpz_srcptr p);
void fft_do_ift_fppol(mp_limb_t * x, mp_size_t nx, void * y, void * temp, struct fft_transform_info * fti, mpz_srcptr p);

/* indicates whether the integer returned is actually reduced modulo some
 * B^n-a, with a=\pm1. Returns n, and sets a. If the result is known to
 * be valid in Z, then n is returned as 0.
 */
static inline mp_bitcnt_t fft_get_mulmod(struct fft_transform_info * fti, int * a)
{
    *a=1;
    return fti->minwrap ? (4<<fti->depth)*fti->bits : 0;
}

static inline mp_size_t fft_get_mulmod_output_minlimbs(struct fft_transform_info * fti)
{
    if (!fti->minwrap) return 0;
    mp_size_t w = fti->w;
    mp_size_t n = 1 << fti->depth;
    mp_bitcnt_t need = (4*n-1)*fti->bits+n*w;
    return (need + FLINT_BITS - 1) / FLINT_BITS;
}


#ifdef __cplusplus
}
#endif

#endif
