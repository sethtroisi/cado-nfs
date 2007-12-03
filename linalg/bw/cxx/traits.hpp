#ifndef TRAITS_HPP_
#define TRAITS_HPP_

#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include "addmul.hpp"

template<mp_size_t width>
struct typical_scalar_traits {
	struct scalar_t { mp_limb_t p[width]; };
	struct wide_scalar_t { mp_limb_t p[width+2]; };
	static inline mpz_class reduce(wide_scalar_t & t) {
		return core_ops::reduce<width>(t.p);
	}
	static inline void reduce(scalar_t & d, wide_scalar_t & s) {
		core_ops::reduce<width>(d.p,s.p);
	}
	static inline void reduce(scalar_t * dst, wide_scalar_t * src,
		unsigned int i0, unsigned int i1)
	{
		for(unsigned int i = i0 ; i < i1 ; i++) {
			reduce(dst[i], src[i - i0]);
		}
	}
	static inline void addmul(wide_scalar_t & dst, scalar_t const & src, int32_t x) {
		core_ops::addmul<width>(dst.p, src.p, x);
	}

	static inline void zero(scalar_t & x) {
		core_ops::zero<width>(x.p);
	}
	static inline void zero(wide_scalar_t & x) {
		core_ops::zero<width+2>(x.p);
	}
	static inline bool is_zero(scalar_t const& x) {
		for(unsigned int i = 0 ; i < width ; i++) {
			if (x.p[i])
				return false;
		}
		return true;
	}
	/*
	static inline void wzero(wide_scalar_t & x) {
		core_ops::zero<width+2>(x.p);
	}
	*/
	static inline void assign(scalar_t & x, const mpz_class& z) {
		core_ops::assign<width>(x.p, z);
	}
	static inline void assign(mpz_class& z, scalar_t const & x) {
		MPZ_SET_MPN(z.get_mpz_t(), x.p, width);
	}
	static inline void addmul(scalar_t & dst,
			scalar_t const & a,
			scalar_t const & b)
	{
		mp_limb_t t[2 * width + 1];
		mp_limb_t c;
		mpn_mul_n(t, a.p, b.p, width);
		c = mpn_add_n(t, t, dst.p, width);
		if (c)
			c = mpn_add_1(t + width, t + width, width, c);
		t[2 * width] = c;
		core_ops::assign<width, 2 * width + 1>(dst.p, t);
	}
	static inline void
	assign(scalar_t & x, std::vector<mpz_class> const& z,  unsigned int i)
	{
		/* This one affects from a vector. We provide the
		 * position, in case we're doing SIMD.
		 */
		assign(x,z[i]);
	}

	/*
	template<mp_size_t nw>
	inline void assign(scalar_t dst, const mp_limb_t z[nw]) {
		return core_ops::assign<width, nw>(dst, z);
	}
	*/
#if 0
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
/*
		std::ostream& operator<<(std::ostream& o, scalar_t const& x)
		{
			mpz_class z;
			MPZ_SET_MPN(z.get_mpz_t(), x.p, width);
			o << z;
			return o;
		}
	*/
	static std::ostream& print(std::ostream& o, scalar_t const& x) {
		mpz_class z;
		MPZ_SET_MPN(z.get_mpz_t(), x.p, width);
		o << z;
		return o;
	}
};
#endif	/* TRAITS_HPP_ */
