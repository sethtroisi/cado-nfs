#ifndef SSE2_8WORDS_TRAITS_HPP_
#define SSE2_8WORDS_TRAITS_HPP_

#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include "addmul.hpp"

namespace globals {
	extern uint8_t modulus_u8;
	extern uint16_t modulus_u16;
	extern uint32_t modulus_u32;
}

struct sse2_8words_traits {
	static const int max_accumulate = 250;
	typedef unsigned short sse2_scalars __attribute__((vector_size(16)));

	struct scalar_t { sse2_scalars p; };
	typedef scalar_t wide_scalar_t;

	static inline mpz_class get_y(scalar_t const & x, int i) {
		unsigned short foo[8] __attribute__((aligned(16)));
		memcpy(foo, &x.p, sizeof(sse2_scalars));
		BUG_ON(i >= 8 || i < 0);
		return mpz_class(foo[i]);
	}
	static inline void reduce(scalar_t & d, wide_scalar_t & s) {
		unsigned short foo[8] __attribute__((aligned(16)));
		memcpy(foo, &s.p, sizeof(sse2_scalars));
		for(unsigned int i = 0 ; i < 8 ; i++) {
			foo[i] %= globals::modulus_u8;
		}
		d.p = (sse2_scalars) {
			foo[0], foo[1], foo[2], foo[3],
			foo[4], foo[5], foo[6], foo[7],
		};
	}
	static inline void reduce(scalar_t * dst, wide_scalar_t * src,
		unsigned int i0, unsigned int i1)
	{
		for(unsigned int i = i0 ; i < i1 ; i++) {
			reduce(dst[i], src[i - i0]);
		}
	}
	static inline void addmul(wide_scalar_t & dst, scalar_t const & src, int32_t x) {
		sse2_scalars w = (sse2_scalars) { x,x,x,x,x,x,x,x, };
		dst.p += src.p * w;
	}

	static inline void zero(scalar_t * p, unsigned int n) {
		memset(p, 0, n * sizeof(scalar_t));
	}
	static inline void zero(scalar_t & x) { zero(&x, 1); }
	static inline void copy(scalar_t * q, const scalar_t * p, unsigned int n) {
		memcpy(q, p, n * sizeof(scalar_t));
	}

	static inline bool is_zero(scalar_t const& x) {
		unsigned short foo[8] __attribute__((aligned(16)));
		memcpy(foo, &x.p, sizeof(sse2_scalars));
		for(unsigned int i = 0 ; i < 8 ; i++) {
			if (foo[i])
				return false;
		}
		return true;
	}
	/*
	static inline void assign(mpz_class& z, scalar_t const & x) {
		MPZ_SET_MPN(z.get_mpz_t(), x.p, width);
	}
	*/
	static inline void addmul(scalar_t & dst,
			scalar_t const & a,
			scalar_t const & b)
	{
		dst.p += a.p * b.p;
	}
	static inline void
	assign(scalar_t & x, std::vector<mpz_class> const& z,  unsigned int i)
	{
		/* This one affects from a vector. We provide the
		 * position, in case we're doing SIMD.
		 */
		x.p = (sse2_scalars) {
			z[i].get_si(),
			z[i+1].get_si(),
			z[i+2].get_si(),
			z[i+3].get_si(),
			z[i+4].get_si(),
			z[i+5].get_si(),
			z[i+6].get_si(),
			z[i+7].get_si(),
		};
	}

	static std::ostream& print(std::ostream& o, scalar_t const& x) {
		unsigned short foo[8] __attribute__((aligned(16)));
		memcpy(foo, &x.p, sizeof(sse2_scalars));
		for(int i = 0 ; i < 8 ; i++) {
			if (i) { o << " "; }
			o << foo[i];
		}
		return o;
	}
};

#endif	/* SSE2_8WORDS_TRAITS_HPP_ */
