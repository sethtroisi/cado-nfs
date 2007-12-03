#ifndef ADDMUL_HPP_
#define ADDMUL_HPP_

#include <gmp.h>
#include <gmpxx.h>
#include "gmp-hacks.h"

#include <cstdlib>

namespace globals {
	extern mpz_class modulus;
}


namespace core_ops {
/* CORE OPS : reduce & addmul */
/* reduce only a scalar */
template<mp_size_t width>
mpz_class reduce(mp_limb_t tmp[width + 2])
{
	mpz_class p;
	mpz_class r;

	if (tmp[width + 1] >> mp_bits_per_limb) {
		mp_size_t i;
		for(i = 0 ; i < width + 2 && !(tmp[i] = -tmp[i]) ; i++);
		for( i++  ; i < width + 2 ; i++) tmp[i] = ~ tmp[i];
		MPZ_SET_MPN(p.get_mpz_t(), tmp, width + 2);
		SIZ(p.get_mpz_t()) *= -1;
	} else {
		MPZ_SET_MPN(p.get_mpz_t(), tmp, width + 2);
	}

	/* FIXME: improve this */
	mpz_fdiv_r(r.get_mpz_t(), p.get_mpz_t(), globals::modulus.get_mpz_t());

	return r;
}

template<mp_size_t width>
void reduce(mp_limb_t dst[width], mp_limb_t src[width + 2])
{
	mpz_class r;
	r = reduce<width>(src);
	BUG_ON(r < 0);
	BUG_ON(r >= globals::modulus);
	MPN_SET_MPZ(dst, width, r.get_mpz_t());
}

/* reduce a vector */
template<mp_size_t width>
void reduce(mp_limb_t dst[][width], mp_limb_t src[][width + 2],
		unsigned int i0, unsigned int i1)
{
	for(unsigned int i = i0 ; i < i1 ; i++) {
		reduce<width>(dst[i], src[i - i0]);
	}
}

template<mp_size_t width>  
void addmul(mp_limb_t dst[width + 2], const mp_limb_t src[width], int32_t x)
{
	mp_limb_t c;
	if (x < 0) {
		c = mpn_submul_1(dst, src, width, -x);
		mpn_sub_1(dst + width, dst + width, 2, c);
	} else {
		c = mpn_addmul_1(dst, src, width, x);
		mpn_add_1(dst + width, dst + width, 2, c);
	}
}
template<mp_size_t width>  
inline void zero(mp_limb_t dst[width])
{
	memset(dst, 0, sizeof(mp_limb_t[width]));
}

template<mp_size_t width>  
inline void assign(mp_limb_t dst[width], const mpz_class& z)
{
	if (z < 0 || z >= globals::modulus) {
		mpz_class r;
		mpz_fdiv_r(r.get_mpz_t(), z.get_mpz_t(),
				globals::modulus.get_mpz_t());
		MPN_SET_MPZ(dst, width, r.get_mpz_t());
	} else {
		MPN_SET_MPZ(dst, width, z.get_mpz_t());
	}
}

#if 0
/* It'd be nice if default arguments for function templates were allowed */
template<mp_size_t width, mp_size_t nw, bool small = (nw < width)>
inline void assign(mp_limb_t dst[width], const mp_limb_t z[nw]);

template<mp_size_t width, mp_size_t nw>
inline void assign<width, nw, true>(mp_limb_t dst[width], const mp_limb_t z[nw])
{
	memcpy(dst, z, sizeof(mp_limb_t[nw]));
	memset(dst + nw, 0, sizeof(mp_limb_t[width - nw]));
}

template<mp_size_t width, mp_size_t nw>
inline void assign<width, nw, false>(mp_limb_t dst[width],
		const mp_limb_t z[nw])
{
	mpz_class r, s;
	MPZ_SET_MPN(s.get_mpz_t(), z, nw);
	mpz_fdiv_r(r.get_mpz_t(), s.get_mpz_t(),
			globals::modulus.get_mpz_t());
	MPN_SET_MPZ(dst, width, r.get_mpz_t());
}
#else
template<mp_size_t width, mp_size_t nw>
inline void assign(mp_limb_t dst[width], const mp_limb_t z[nw])
{
	if (nw < width) {
		memcpy(dst, z, sizeof(mp_limb_t[nw]));
		/* As the code is always generated here, can't use
		 * mp_limb_t[width-nw] as a type */
		memset(dst + nw, 0, (width - nw) * sizeof(mp_limb_t));
	} else {
		mpz_class r, s;
		MPZ_SET_MPN(s.get_mpz_t(), z, nw);
		mpz_fdiv_r(r.get_mpz_t(), s.get_mpz_t(),
				globals::modulus.get_mpz_t());
		MPN_SET_MPZ(dst, width, r.get_mpz_t());
	}
}
#endif

}

#endif	/* ADDMUL_HPP_ */
