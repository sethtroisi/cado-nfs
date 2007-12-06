#ifndef SSE2_2WORDS_TRAITS_HPP_
#define SSE2_2WORDS_TRAITS_HPP_

#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <stdint.h>
#include "traits_globals.hpp"
#include "matrix_repr_prime.hpp"

struct sse2_2words_traits {
	typedef matrix_repr_prime representation;

#if defined(__OpenBSD__) && defined(__x86_64)
	/* I know it's outright stupid */
	typedef long inner_type;
#else
	typedef int64_t inner_type;
#endif
	static const int sse2_vectoring = 2;

	static const int max_accumulate = 50;
	static const int max_accumulate_wide = 50;
	typedef inner_type sse2_scalars __attribute__((vector_size(16)));
	typedef inner_type vec_t[sse2_vectoring] __attribute__((aligned(16)));

	typedef prime_field_any coeff_field;

	struct scalar_t { sse2_scalars p; };
	typedef scalar_t wide_scalar_t;

	struct name {
		const char * operator()() const {
			return "SSE-2 code [ 2 words ]";
		}
	};

	static int can() {
		return globals::nbys == 2 && MODBITS() < 20;
	}

	static inline mpz_class get_y(scalar_t const & x, int i) {
		vec_t foo;
		memcpy(foo, &x.p, sizeof(sse2_scalars));
		BUG_ON(i >= sse2_vectoring || i < 0);
		return mpz_class(foo[i]);
	}
	static inline void reduce(scalar_t & d, wide_scalar_t & s) {
		vec_t foo;
		memcpy(foo, &s.p, sizeof(sse2_scalars));
		for(int i = 0 ; i < sse2_vectoring ; i++) {
			/* FIXME : make this generic here */
			foo[i] %= globals::modulus_u32;
		}
		/* FIXME : make this generic here */
		d.p = (sse2_scalars) {
			foo[0], foo[1],
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
		/* FIXME : make this generic here */
		sse2_scalars w = (sse2_scalars) { x,x, };
		dst.p += src.p * w;
	}

	static inline void zero(scalar_t * p, unsigned int n) {
		memset(p, 0, n * sizeof(scalar_t));
	}
	static inline void zero(scalar_t & x) { zero(&x, 1); }
	static inline void copy(scalar_t * q, const scalar_t * p, unsigned int n) {
		memcpy(q, p, n * sizeof(scalar_t));
	}
        /*      [duplicate]
	static inline void copy(wide_scalar_t * q, const wide_scalar_t * p, unsigned int n) {
		memcpy(q, p, n * sizeof(wide_scalar_t));
	}
        */


	static inline bool is_zero(scalar_t const& x) {
		vec_t foo;
		memcpy(foo, &x.p, sizeof(sse2_scalars));
		for(int i = 0 ; i < sse2_vectoring ; i++) {
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
	static inline void addmul_wide(scalar_t & dst,
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
		/* FIXME : make this generic here */
		x.p = (sse2_scalars) {
			z[i].get_si(),
			z[i+1].get_si(),
		};
	}
	static inline void assign(std::vector<mpz_class>& z, scalar_t const & x) {
		vec_t foo;
		memcpy(foo, &x.p, sizeof(sse2_scalars));
		BUG_ON(z.size() != (unsigned int) sse2_vectoring);
		for(int i = 0 ; i < sse2_vectoring ; i++) {
			z[i] = foo[i];
		}
	}

	static std::ostream& print(std::ostream& o, scalar_t const& x) {
		vec_t foo;
		memcpy(foo, &x.p, sizeof(sse2_scalars));
		for(int i = 0 ; i < sse2_vectoring ; i++) {
			if (i) { o << " "; }
			o << foo[i];
		}
		return o;
	}
};

#endif	/* SSE2_2WORDS_TRAITS_HPP_ */
