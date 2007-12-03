#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "field_def.h"
#include "field_common.h"
#include "field_prime.h"
#include "field_usage.h"
#include "variables.h"

#include "ops_poly_wrap.hpp"

void ops_poly_wrap::itransform  (
		mp_limb_t * dst, ptrdiff_t stride, int deg,
		transform_t p) const
{
	int i;

	scbk.itransform(dst, stride, deg, p.scbk);

	{
		mp_limb_t * res_fft;
		res_fft = (mp_limb_t *) malloc((deg+1)*bw_allocsize*sizeof(mp_limb_t));
		fft.itransform(res_fft, bw_allocsize, deg, p.fft);
		for(i = 0 ; i <= deg ; i++) {
			mp_limb_t * ref = dst + i * stride;
			mp_limb_t * me = res_fft + i * bw_allocsize;
			if (!k_eq(me, ref)) {
				fprintf(stderr, "fft result bad, coeff %d\n",
						i);
				abort();
			}
		}
		free(res_fft);
	}

	{
		mp_limb_t * res_ifft;
		res_ifft = (mp_limb_t *) malloc((deg+1)*bw_allocsize*sizeof(mp_limb_t));
		ifft.itransform(res_ifft, bw_allocsize, deg, p.ifft);
		for(i = 0 ; i <= deg ; i++) {
			mp_limb_t * ref = dst + i * stride;
			mp_limb_t * me = res_ifft + i * bw_allocsize;
			if (!k_eq(me, ref)) {
				fprintf(stderr, "ifft result bad, coeff %d\n",
						i);
				abort();
			}
		}
		free(res_ifft);
	}
	fprintf(stderr, "checked degree %d ok\n", deg);
}
