#ifndef ULONG_TRAITS_HPP_
#define ULONG_TRAITS_HPP_

#include <gmp.h>
#include <gmpxx.h>
#include <vector>
#include <limits.h>
#include "traits_globals.hpp"
#include "matrix_repr_prime.hpp"

template <typename T, int acc>
struct pod_traits {
	typedef matrix_repr_prime representation;

	static const int max_accumulate = acc;
	static const int max_accumulate_wide = acc;
	static const int nbits = sizeof(T)*CHAR_BIT;

	typedef prime_field_any coeff_field;
	struct scalar_t { T p; };
	typedef scalar_t wide_scalar_t;

	struct name {
		char x[40];
		name() {
			snprintf(x,sizeof(x),"pod-type (%d bits) code", nbits);
		}
		const char * operator()() const {
			return x;
		}
	};

	static int can() {
		return globals::nbys == 1 &&
			nbits >= 2*(int) MODBITS() &&
			acc < 1UL << (nbits - 2*MODBITS()); 
	}

	static inline mpz_class get_y(scalar_t const & x, int i) {
		BUG_ON(i);
		return mpz_class(x.p);
	}
	static inline void reduce(scalar_t & d, wide_scalar_t & s) {
		d.p = s.p % (T) globals::modulus_ulong;
	}
	static inline void reduce(scalar_t * dst, wide_scalar_t * src,
		unsigned int i0, unsigned int i1)
	{
		for(unsigned int i = i0 ; i < i1 ; i++) {
			reduce(dst[i], src[i - i0]);
		}
	}
	static inline void addmul(wide_scalar_t & dst, scalar_t const & src, int32_t x) {
		dst.p += src.p * x;
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
		return x.p == 0;
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
		x.p = z[i].get_si();
	}
	static inline void assign(std::vector<mpz_class>& z, scalar_t const & x) {
		BUG_ON(z.size() != 1);
		z[0] = x.p;
	}

	static std::ostream& print(std::ostream& o, scalar_t const& x) {
		return o << x.p;
	}
};

typedef pod_traits<long, 250> ulong_traits;
typedef pod_traits<int64_t, 250> int64_traits;
typedef pod_traits<int32_t, 50> int32_traits;
typedef pod_traits<int16_t, 10> int16_traits;
typedef pod_traits<int8_t, 6> int8_traits;

#endif	/* ULONG_TRAITS_HPP_ */
