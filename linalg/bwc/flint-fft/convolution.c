/* 
 * Copyright (C) 2008-2011 William Hart
 * 
 * This file is part of FLINT.
 * 
 * FLINT is free software: you can redistribute it and/or modify it under the 
 * terms of the GNU Lesser General Public License (LGPL) as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.  See <http://www.gnu.org/licenses/>. */

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
// #include "fmpz.h"
// #include "fmpz_vec.h"
// #include "fmpz_poly.h"
#include "fft.h"

/*
 * Perform an FFT convolution of ``ii`` with ``jj``, both of length 
 * ``4*n`` where ``n = 2^depth``. Assume that all but the first 
 * ``trunc`` coefficients of the output (placed in ``ii``) are zero.
 * Each coefficient is taken modulo ``B^limbs + 1``. The temporary 
 * spaces ``t1``, ``t2`` and ``s1`` must have ``limbs + 1`` 
 * limbs of space and ``tt`` must have ``2*(limbs + 1)`` of free 
 * space.
 * 
 */
void fft_convolution(mp_limb_t ** ii, mp_limb_t ** jj, slong depth,
		     slong limbs, slong trunc, mp_limb_t ** t1,
		     mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt)
{
    slong n = (WORD(1) << depth), j;
    slong w = (limbs * FLINT_BITS) / n;
    slong sqrt = (WORD(1) << (depth / 2));

    if (depth <= 6) {
	trunc = 2 * ((trunc + 1) / 2);

	fft_truncate_sqrt2(ii, n, w, t1, t2, s1, trunc);

	if (ii != jj)
	    fft_truncate_sqrt2(jj, n, w, t1, t2, s1, trunc);

	for (j = 0; j < trunc; j++) {
	    mpn_normmod_2expp1(ii[j], limbs);
	    if (ii != jj)
		mpn_normmod_2expp1(jj[j], limbs);

	    fft_mulmod_2expp1(ii[j], ii[j], jj[j], n, w, *tt);
	}

	ifft_truncate_sqrt2(ii, n, w, t1, t2, s1, trunc);

	for (j = 0; j < trunc; j++) {
	    mpn_div_2expmod_2expp1(ii[j], ii[j], limbs, depth + 2);
	    mpn_normmod_2expp1(ii[j], limbs);
	}
    } else {
	trunc = 2 * sqrt * ((trunc + 2 * sqrt - 1) / (2 * sqrt));

	fft_mfa_truncate_sqrt2_outer(ii, n, w, t1, t2, s1, sqrt, trunc);

	if (ii != jj)
	    fft_mfa_truncate_sqrt2_outer(jj, n, w, t1, t2, s1, sqrt, trunc);

	fft_mfa_truncate_sqrt2_inner(ii, jj, n, w, t1, t2, s1, sqrt, trunc,
				     tt);

	ifft_mfa_truncate_sqrt2_outer(ii, n, w, t1, t2, s1, sqrt, trunc);
    }
}
